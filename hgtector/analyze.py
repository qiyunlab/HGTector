#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re
from os import makedirs
from os.path import join, isdir, isfile, dirname
from io import StringIO

import numpy as np
import pandas as pd

from scipy.signal import find_peaks
from scipy.stats import zscore

from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KernelDensity, NearestCentroid
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import silhouette_samples

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hgtector.util import (
    timestamp, load_configs, get_config, file2id, id2file_map, arg2bool,
    read_taxdump, read_file, list_from_param, dict_from_param, get_lineage,
    sort_by_hierarchy, refine_taxdump, add_children, describe_taxon, find_lca,
    taxid_at_rank, get_descendants, save_figure)


description = """predict HGT events based on search results"""

arguments = [
    'basic',
    ['-i|--input',   'input search result file, or directory where one or '
                     'more input files are located'],
    ['-o|--output',  'directory where analysis results will be saved',
                     {'required': True}],
    ['-t|--taxdump', 'directory of taxonomy database files (nodes.dmp and '
                     'names.dmp); required if they are not found in input '
                     'directory'],

    'hit filtering',
    ['-k|--maxhits', 'maximum number of sequence similarity search hits per '
                     'gene (protein) to preserve', {'type': int}],
    ['--evalue',     'maximum E-value cutoff', {'type': float}],
    ['--identity',   'minimum percent identity cutoff', {'type': int}],
    ['--coverage',   'minimum percent query coverage cutoff', {'type': int}],

    'taxonomic assignment',
    ['--input-tax',  'taxIds of input samples, format: sample1:taxId1,sample2:'
                     'taxId2..., will auto-infer if omitted'],
    ['--input-cov',  'for auto-inference: an input sample will be assigned '
                     'to a taxon if it represents this percentage or more '
                     'best hits', {'type': int}],

    'taxonomic grouping',
    ['--self-tax',   'taxIds of "self" group, delimited by comma, will auto-'
                     'infer if omitted'],
    ['--close-tax',  'taxIds of "close" group, delimited by comma, will auto-'
                     'infer if omitted'],
    ['--self-rank',  'for auto-inference: "self" group must be at or above '
                     'this rank'],
    ['--close-size', 'for auto-inference: "close" group must have at least '
                     'this number of taxa', {'type': int}],
    ['--distal-top', 'find match in "distal" group which is LCA of hits with '
                     'bit score at most this percentage lower than best hit',
                     {'type': int}],

    'scoring',
    ['--weighted',   'score is sum of weighted bit scores; otherwise simple '
                     'counts', {'choices': ['yes', 'no']}],

    'filtering',
    ['--outliers',   'detect and remove outliers using selected method',
                     {'choices': ['none', 'zscore', 'boxplot']}],
    ['--orphans',    'keep orphans (proteins without non-self hits) in '
                     'statistical analysis',
                     {'choices': ['yes', 'no']}],

    'prediction',
    ['--bandwidth',  'bandwidth for Gaussian KDE (auto, grid, silverman, or a '
                     'number between 0.1 and 1.0)'],
    ['--bw-steps',   'steps for auto and grid kernel bandwidth optimization',
                     {'type': int}],
    ['--low-part',   'maximum percentage below threshold for automatic '
                     'bandwidth optimization', {'type': float}],
    ['--noise',      'percent valley-to-peak distance to exclude from cluster',
                     {'type': float}],
    ['--fixed',      'use this percentage as threshold if KDE clustering '
                     'fails', {'type': float}],
    ['--silhouette', 'silhouette score threshold for cluster refinement',
                     {'type': float}],
    ['--self-low',   'HGT has low "self" score (an optional criterion)',
                     {'choices': ['yes', 'no']}],

    'program behavior',
    ['--from-scores', 'if score table already exists, use it and skip search '
                      'result parsing and taxonomy inference; otherwise '
                      'overwrite it', {'action': 'store_true'}],
]


class Analyze(object):

    def __init__(self):
        self.arguments = arguments
        self.description = description

    def __call__(self, args):
        print(f'Analysis started at {timestamp()}.')

        # load configurations
        self.cfg = load_configs()

        # read and validate arguments
        self.set_parameters(args)

        # use existing score table
        score_file = join(self.output, 'scores.tsv')
        if self.from_scores and isfile(score_file):
            self.df = pd.read_csv(score_file, sep='\t')
        else:

            # read input taxonomy and search results
            self.read_input()

            # assign taxonomy to input genomes
            self.assign_taxonomy()

            # define the three groups for search hits
            self.define_groups()

            # calculate scores for the three groups per protein
            self.calc_scores()

            # generate a table of calculated scores
            self.make_score_table()

        # no distal hits
        if not self.df['distal'].any():
            print('WARNING: No hit is assigned to distal group. Cannot '
                  'predict HGTs.')
            return

        # remove orphans
        if not self.orphans:
            self.remove_orphans()

        # remove outliers
        if self.outliers != 'none':
            self.remove_outliers()

        # predict HGTs
        self.predict_hgt()

        print(f'Analysis finished at {timestamp()}.')

    def set_parameters(self, args):
        """Validate and set parameters.

        Parameters
        ----------
        args : dict
            command-line arguments

        Raises
        ------
        ValueError
            Found invalid parameter(s).
        """
        # load arguments
        for key, val in vars(args).items():
            setattr(self, key, val)

        # check input directory and data
        if self.input:
            if isfile(self.input):
                self.input_map = {file2id(self.input): self.input}
            elif isdir(self.input):
                self.input_map = {k: join(self.input, v) for k, v in
                                  id2file_map(self.input, ext='tsv').items()}
            else:
                raise ValueError(
                    f'Invalid input data file or directory: {self.input}.')
            if len(self.input_map) == 0:
                raise ValueError(
                    f'No input data are found under: {self.input}.')

        # check / create output directory
        makedirs(self.output, exist_ok=True)
        self.prev_map = id2file_map(self.output, 'tsv')

        # load configurations
        get_config(self, 'evalue', 'analyze.evalue', float)
        for key in ('maxhits', 'identity', 'coverage'):
            get_config(self, key, f'analyze.{key}')
        for key in ('input_cov', 'self_rank', 'close_size', 'distal_top'):
            get_config(self, key, f'grouping.{key.replace("_", "")}')
        for key in ('weighted', 'outliers', 'orphans', 'bandwidth', 'bw_steps',
                    'low_part', 'noise', 'fixed', 'silhouette', 'self_low'):
            get_config(self, key, f'predict.{key.replace("_", "")}')

        # convert boolean values
        for key in ('weighted', 'orphans', 'self_low'):
            setattr(self, key, arg2bool(getattr(self, key, None)))

        # convert fractions to percentages
        for metric in ('input_cov', 'noise', 'fixed', 'distal_top'):
            val = getattr(self, metric)
            if val and val < 1:
                setattr(self, metric, val * 100)

        # convert distal top to a factor to save compute
        self.match_th = 1 - self.distal_top / 100

        # force coverage >= 50 to ensure that candidates are sequential
        if (self.input_cov or 0) < 50:
            raise ValueError('Taxonomy coverage for auto-interence must be at '
                             'least 50%.')

    def read_input(self):
        """Workflow for reading input data.

        Notes
        -----
        1. Read taxonomy database into `taxdump`.
        2. Read homology search results into `data`.
        """
        # read taxonomy database
        if self.taxdump is not None:
            print('Reading local taxonomy database...')
            self.taxdump = read_taxdump(self.taxdump)
        elif (isfile(join(self.input, 'names.dmp')) and
                isfile(join(self.input, 'nodes.dmp'))):
            print('Reading custom taxonomy database...')
            self.taxdump = read_taxdump(self.input)
        elif (isfile(join(dirname(self.input), 'names.dmp')) and
                isfile(join(dirname(self.input), 'nodes.dmp'))):
            print('Reading custom taxonomy database...')
            self.taxdump = read_taxdump(dirname(self.input))
        else:
            raise ValueError('Missing taxonomy database.')
        print(f'Done. Read {len(self.taxdump)} taxa.')

        # read search results
        print('Reading homology search results...')
        self.data = {}
        for sid, fname in self.input_map.items():
            self.data[sid] = self.read_search_results(fname)
            print(f'  {sid}: {len(self.data[sid])} proteins.')
        print(f'Done. Read search results of {len(self.data)} samples.')

    @staticmethod
    def read_search_results(file, maxhits=None):
        """Read homology search results of one sample.

        Parameters
        ----------
        file : str
            input filepath
        maxhits : int
            maximum number of hits per protein to preserve

        Returns
        -------
        list of dict
            search results
        """
        p = re.compile(r'# (\S+): (.*)')
        data = []
        with read_file(file) as f:
            for line in f:
                line = line.rstrip('\r\n')
                m = p.match(line)
                if m:
                    if m.group(1) == 'ID':
                        data.append({'id': m.group(2), 'hits': []})
                    elif m.group(1) == 'Length':
                        data[-1]['length'] = int(m.group(2))
                    elif m.group(1) == 'Product':
                        data[-1]['product'] = m.group(2)
                    elif m.group(1) == 'Score':
                        data[-1]['score'] = float(m.group(2))
                else:
                    data[-1]['hits'].append(line)
                    if len(data[-1]['hits']) == (maxhits or 0):
                        break

        # convert hit table to DataFrame
        for i in range(len(data)):
            data[i]['hits'] = pd.read_csv(
                StringIO('\n'.join(data[i]['hits'])), sep='\t', na_values='*',
                names=['id', 'identity', 'evalue', 'score', 'coverage',
                       'taxid'],
                dtype={'id': str,
                       'identity': np.float32,
                       'evalue': np.float64,
                       'score': np.float32,
                       'coverage': np.float32,
                       'taxid': str}).set_index('id')
        return data

    def assign_taxonomy(self):
        """Assign taxonomy to genomes.
        """
        # take user-defined taxIds of input genomes
        if self.input_tax:
            try:
                self.input_tax = dict_from_param(self.input_tax)
            except ValueError:
                if len(self.data) > 1:
                    raise ValueError('Invalid input taxonomy format.')
                # for single-sample analysis, one can simply enter a taxId
                self.input_tax = {max(self.data.keys()): self.input_tax}
            print('User-specified TaxIDs of input genomes:')
            for sid, tid in sorted(self.input_tax.items()):
                if tid not in self.taxdump:
                    # TODO: read from both temp and master taxdump
                    raise ValueError(
                        f'TaxID {tid} is not present in taxonomy database.')
                print(f'  {sid}: {tid} ({self.taxdump[tid]["name"]}).')
        else:
            self.input_tax = {}

        # auto-infer taxIds of remaining genomes
        sids = sorted([x for x in self.data if x not in self.input_tax])
        if sids:
            print('Auto-inferring plausible taxIds for input genomes based on '
                  'taxonomy of search results...')
            for sid in sids:
                try:
                    tid, cov = self.infer_genome_tax(
                        self.data[sid], self.taxdump, self.input_cov)
                    self.input_tax[sid] = tid
                except ValueError:
                    raise ValueError(f'Cannot auto-infer taxonomy for {sid}.'
                                     ' Please specify manually.')
                print(f'  {sid}: {tid} ({self.taxdump[tid]["name"]}) '
                      f'(covering {cov:2g}% best hits).')

        # refine taxonomy database
        print('Refining taxonomy database...')
        self.taxdump = refine_taxdump(self.sum_taxids(), self.taxdump)
        add_children(self.taxdump)
        print(f'Done. Retained {len(self.taxdump)} taxa.')

        # find lowest common ancestor (LCA) of all genomes
        self.lca = find_lca(self.input_tax.values(), self.taxdump)
        print(f'All input genomes belong to {self.lca} '
              f'({describe_taxon(self.lca, self.taxdump)}).')

    @staticmethod
    def infer_genome_tax(prots, taxdump, coverage):
        """Automatically infer taxId of a genome based on best hits.

        Parameters
        ----------
        prots : list of dicts
            proteins of a genome, including hit tables
        taxdump : dict
            taxonomy database
        coverage : int
            minimum percentage coverage of best hits

        Returns
        -------
        str, float
            taxId, and percentage of best hits that belong to this taxId

        Raises
        ------
        ValueError
            Cannot infer taxon (e.g., none of taxIds cover more best hits than
            threshold).
        """
        # collect taxIds of best hits of all proteins
        bestids = []
        for prot in prots:
            if prot['hits'].shape[0]:
                bestids.append(prot['hits']['taxid'].iloc[0])

        # count proteins that have at least one hit
        n = len(bestids)

        # get frequencies of all taxIds at or above taxId of each best hit
        freqs = {}
        for bestid in bestids:
            for tid in get_lineage(bestid, taxdump):
                freqs[tid] = freqs.get(tid, 0) + 1

        # loop frequencies from low to high, and select taxId(s) at the lowest
        # frequency
        lowfreq = 0
        candidates = []
        for tid, freq in sorted(freqs.items(), key=lambda x: x[1]):
            if freq / n * 100 >= coverage:
                if lowfreq == 0:
                    lowfreq = freq
                elif freq > lowfreq:
                    break
                candidates.append(tid)

        # return the lowest-frequency taxId
        m = len(candidates)
        if not m:
            raise ValueError('Cannot auto-infer taxonomy.')
        elif m == 1:
            res = candidates[0]
            return res, freqs[res] / n * 100

        # if there are multiple, sort by hierarchy, and take the lowest one
        else:
            res = sort_by_hierarchy(candidates, taxdump)[0]
            return res, freqs[res] / n * 100

    def sum_taxids(self):
        """Get all taxIds mentioned in all sets and hit tables.

        Returns
        -------
        set of str
            all taxIds
        """
        res = set(self.input_tax.values())
        for sid, prots in self.data.items():
            for prot in prots:
                res.update(prot['hits']['taxid'].tolist())
        return res

    def define_groups(self):
        """Define the three (actually two) groups: "self" and "close".

        Notes
        -----
        Assign these attributes:
        1. `self_tax`: top-level taxId(s) of the self group.
        2. `close_tax`: top-level taxId(s) of the close group.
        3. `groups` (keys: self, close, distal): all taxIds under each group.
        """
        self.groups = {}
        for key in ('self', 'close'):
            tids = getattr(self, f'{key}_tax')

            # user-defined group
            if tids:
                setattr(self, f'{key}_tax', list_from_param(tids))
                print(f'User-defined {key} group:')

            # auto-infer group
            else:
                getattr(self, f'infer_{key}_group')()
                print(f'Auto-inferred {key} group:')

            # collect taxIds that belong to group
            tids = getattr(self, f'{key}_tax')
            if key not in self.groups:
                self.groups[key] = set().union(*[[x] + get_descendants(
                    x, self.taxdump) for x in tids])

                # subtract self group from close group
                if key == 'close':
                    self.groups['close'] = self.groups['close'].difference(
                        self.groups['self'])

            # report group content
            for tid in tids:
                print(f'  {tid} ({describe_taxon(tid, self.taxdump)})')
            print(f'{key.capitalize()} group has {len(self.groups[key])} '
                  'taxa.')

    def infer_self_group(self):
        """Infer self group automatically.

        Notes
        -----
        Assign `self_tax` as top-level taxId(s) of the self group.
        """
        # just use LCA
        if not self.self_rank:
            self.self_tax = [self.lca]

        # try to raise LCA to given rank, but if LCA is already above that
        # rank, just use LCA
        else:
            tid_ = taxid_at_rank(self.lca, self.self_rank, self.taxdump)
            self.self_tax = [tid_ or self.lca]

    def infer_close_group(self):
        """Infer close group automatically.

        Notes
        -----
        1. Assign `close_tax` as top-level taxId(s) of the close group.
        2. Assign `groups['close']` as all taxIds under the close group.
        """
        mems = []

        # start from the LCA of self group
        cid = find_lca(self.self_tax, self.taxdump)
        while True:

            # close group should exclude self group
            mems = set([cid] + get_descendants(
                cid, self.taxdump)).difference(self.groups['self'])

            # stop when size limit is reached
            if mems and (not self.close_size or len(mems) >= self.close_size):
                break

            # move up one level
            pid = self.taxdump[cid]['parent']
            if pid == cid or pid == '0':
                break
            cid = pid
        self.close_tax = [cid]
        self.groups['close'] = mems

    def calc_scores(self):
        """Summarize search scores for proteins.

        Notes
        -----
        1. Append column `group` to the hit table of each protein based on the
          `taxid` column.
        2. Calculate scores of `self`, `close` and `distal` for each protein.
        3. Infer `match` for each protein based on its distal hits.
        """
        print('Calculating protein scores by group...', flush=True)
        for sid, prots in sorted(self.data.items()):
            for prot in prots:
                # assign each hit to one of the three groups
                prot['hits']['group'] = prot['hits']['taxid'].apply(
                    lambda x: 'self' if x in self.groups['self'] else 'close'
                    if x in self.groups['close'] else 'distal')

                # sum up scores per group
                gb = prot['hits'].groupby('group')
                scores = (gb.size() if not self.weighted else gb['score']
                          .sum() / prot['score']).to_dict()
                for group in ('self', 'close', 'distal'):
                    prot[group] = scores[group] if group in scores else 0.0

                # find best match taxId in distal group
                prot['match'] = self.find_match(prot['hits'].query(
                    'group == "distal"'))
            print(f'  {sid}')
        print('Done.')

    def find_match(self, df):
        """Find a taxId that best describes top hits.

        Parameters
        ----------
        df : pd.DataFrame
            hit table

        Returns
        -------
        str
            taxId of match, or '0' if not found

        Notes
        -----
        The best match taxId is the LCA of top hits. The "top hits" are
        defined as those whose bit scores are no less than a certain
        percentage of that of the best hit. This behavior is similar to
        DIAMOND's taxonomic classification function.
        """
        try:
            th = df.iloc[0]['score'] * self.match_th
        except IndexError:
            return '0'
        return find_lca(df[df['score'] >= th]['taxid'], self.taxdump)

    def make_score_table(self):
        """Make a data frame for the entire protein set.

        Notes
        -----
        1. Generate score table (pd.DataFrame) and assign to `df`.
        2. Write score table to file `scores.tsv` (tab-delimited).
        """
        print('Summarizing scores of all proteins...', end='', flush=True)
        self.df = {}
        data = []
        for sid, prots in self.data.items():
            for prot in prots:
                data.append([sid, prot['id'], prot['length'],
                             prot['hits'].shape[0], prot['self'],
                             prot['close'], prot['distal'], prot['match']])
        self.df = pd.DataFrame(data, columns=[
            'sample', 'protein', 'length', 'hits', 'self', 'close', 'distal',
            'match'])
        self.df.to_csv(join(self.output, 'scores.tsv'), sep='\t',
                       index=False, float_format='%g')
        print(' done.')
        print('Protein scores written to scores.tsv.')

    def remove_orphans(self):
        """Remove ORFans (genes without non-self hits).

        Notes
        -----
        Remove ORFan rows from the score table `df` in-place.
        """
        n = self.df.shape[0]
        self.df.query('close + distal > 0', inplace=True)
        print(f'Removed {n - self.df.shape[0]} ORFans.')

    def remove_outliers(self):
        """Remove outliers from selected groups of scores.

        Notes
        -----
        Remove outlier rows from the score table `df` in-place.
        Only outliers at the right (high) side will be removed.
        """
        # TODO: add other methods
        groups = self.relevant_groups()
        n = self.df.shape[0]

        if self.outliers == 'zscore':
            self.df = self.outliers_zscore(self.df, groups)
        elif self.outliers == 'boxplot':
            self.df = self.outliers_boxplot(self.df, groups)

        print(f'Removed {n - self.df.shape[0]} outliers.')

    def relevant_groups(self):
        """Get groups that are relevant in HGT prediction.

        Returns
        -------
        list of str
            relevant groups in order

        Notes
        -----
        The `self` group is relevant only when `self_low` is True, otherwise
        only `close` and `distal` groups are relevant.
        """
        return (['self'] if self.self_low else []) + ['close', 'distal']

    @staticmethod
    def outliers_zscore(df, keys):
        """Remove outliers using the Z-score method.

        Parameters
        ----------
        df : pd.DataFrame
            input dataframe
        keys : list of str
            relevant columns of dataframe

        Returns
        -------
        pd.DataFrame
            output dataframe with outliers removed

        Notes
        -----
        Z-score >= 3 is the criterion for outliers.
        """
        return df[(zscore(df[keys]) < 3).all(axis=1)]

    @staticmethod
    def outliers_boxplot(df, keys):
        """Remove outliers using the boxplot (IQR) method.

        Parameters
        ----------
        df : pd.DataFrame
            input dataframe
        keys : list of str
            relevant columns of dataframe

        Returns
        -------
        pd.DataFrame
            output dataframe with outliers removed

        Notes
        -----
        Criterion for outliers: < Q1 - 1.5 * IQR or > Q3 + 1.5 * IQR.
        """
        return df[(df[keys] <= df[keys].quantile(0.75) * 2.5 -
                   df[keys].quantile(0.25) * 1.5).all(axis=1)]

    def predict_hgt(self):
        """Predict HGTs.

        Returns
        -------
        int
            number of predicted HGTs

        Notes
        -----
        Central pipeline of this script.
        - Report threshold for each group.
        - Generate one scatter plot for close vs distal, and one density plot
          for each group.
        - Generate one text file for putative HGTs of each sample under `hgt/`.
        """
        print('Predicting HGTs...')

        # perform kernel density estimation (KDE), identify "atypical"
        # cluster, and determine threshold
        ths = {}
        groups = self.relevant_groups()
        print('Calculating thresholds for clustering...')
        for group in groups:
            print(f'{group.capitalize()} group:')
            self.plot_hist(self.df[group].tolist(),
                           join(self.output, f'{group}.hist.png'))

            # cannot cluster constant data
            if self.df[group].std() == 0.0:
                print(f'WARNING: {group.capitalize()} group is constant. '
                      'Cannot predict HGTs.')
                return 0

            # calculate threshold using KDE
            ths[group] = self.cluster_kde(group)

            # use a fixed global threshold if KDE fails
            if not ths[group] and self.fixed:
                print(f'WARNING: Cannot cluster {group} group using KDE. '
                      f'Use fixed threshold {self.fixed} instead.')
                ths[group] = self.df[group].quantile(self.fixed / 100)

            print(f'  Threshold: {ths[group]:g}.')
        print('Done.')

        # identify atypical cluster
        print('Labeling cluster...', end='')
        self.df['hgt'] = (self.df['close'] <= ths['close']) &\
                         (self.df['distal'] >= ths['distal']) &\
                         (not self.self_low or self.df['self'] <= ths['self'])
        print(' done.')
        n = self.df[self.df['hgt']].shape[0]
        print(f'  Total predicted HGTs: {n:g}.')
        if not n:
            return 0

        # calculate silhouette scores and centroid
        print('Calculating cluster properties...', end='')
        cent = self.calc_cluster_props()
        print(' done.')

        # refine cluster by silhouette score
        if self.silhouette:
            print('Refining cluster...', end='')
            self.refine_cluster(cent)
            print(' done.')
            n = self.df[self.df['hgt']].shape[0]
            print(f'  Total predicted HGTs after refinement: {n:g}.')
            if not n:
                return 0

        # summarize prediction results
        print('Predicted HGTs by sample:')
        makedirs(join(self.output, 'hgts'), exist_ok=True)
        for sample in self.df['sample'].unique():
            df_ = self.df[self.df['hgt'] & (self.df['sample'] == sample)]
            print(f'  {sample}: {df_.shape[0]}.')
            df_[['protein', 'silh']].to_csv(
                join(self.output, 'hgts', f'{sample}.txt'),
                sep='\t', index=False, header=False, float_format='%g')
        print('Prediction results saved to hgts/.')

        # plot prediction results
        self.plot_hgts()
        return self.df[self.df['hgt']].shape[0]

    def cluster_kde(self, group):
        """Cluster data by KDE.

        Parameters
        ----------
        group : str
            which group to cluster

        Returns
        -------
        float
            clustering threshold, or 0 if not determined
        """
        if self.bandwidth != 'auto':
            data = self.df[group].values

            # kernel density estimation
            x, y, bw = self.perform_kde(data)

            # find first peak and first valley
            try:
                peak, valley = self.first_hill(x, y)
            except ValueError:
                return 0.0

            # determine threshold
            th = valley - (valley - peak) * self.noise / 100

            # plot density function and thresholds
            self.plot_density(x, y, peak, valley, th,
                              join(self.output, f'{group}.kde.png'))

            return th
        else:
            return self.smart_kde(group)

    def perform_kde(self, data):
        """Perform kernel density estimation (KDE) on data.

        Parameters
        ----------
        data : np.array
            sample data (1D)

        Returns
        -------
        np.array, np.array, float
            x values, y values, bandwidth

        .. _scikit-learn tutorial 1:
            https://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d
            .html
        .. _scikit-learn tutorial 2:
            https://scikit-learn.org/stable/auto_examples/neighbors/plot_digits_
            kde_sampling.html
        .. _A useful article:
            https://jakevdp.github.io/blog/2013/12/01/kernel-density-estimation/
        """
        bw = self.bandwidth
        data_ = data[:, np.newaxis]
        scaler = StandardScaler()
        data_ = scaler.fit_transform(data[:, np.newaxis])
        estimator = KernelDensity(kernel='gaussian')

        # grid search optimization
        if bw == 'grid':
            kde = self.grid_kde(data_, estimator, self.bw_steps)
            bw = kde.bandwidth
            print(f'  Grid search-optimized bandwidth: {bw:g}.')

        # Silverman's rule-of-thumb
        elif bw == 'silverman':
            bw = self.silverman_bw(data)
            print(f'  Bandwidth by Silverman\'s rule-of-thumb: {bw:g}.')
            setattr(estimator, 'bandwidth', bw)
            kde = estimator.fit(data_)

        # fixed bandwidth value
        elif isinstance(bw, float) and 0.1 <= bw <= 1.0:
            setattr(estimator, 'bandwidth', bw)
            kde = estimator.fit(data_)
        else:
            raise ValueError(f'Invalid bandwidth: {bw}.')

        # get density function
        x, y = self.density_func(data_, kde)
        x = scaler.inverse_transform(x)
        y = scaler.inverse_transform(y)
        return x, y, bw

    @staticmethod
    def grid_kde(data, estimator, steps):
        """Perform kernel density estimation using grid search with cross
        validations.

        Parameters
        ----------
        data : np.array
            input data
        estimator : sklearn.neighbors.KernelDensity
            kernel density estimator
        steps : int
            number of bandwidths to test

        Returns
        -------
        sklearn.neighbors.KernelDensity
            estimator trained on data

        Raises
        ------
        ValueError
            If data size < 5 (number of splits).
        """
        n = data.size
        if n < 5:
            raise ValueError(
                f'Cannot perform grid search on {n} data point(s).')
        bwspace = np.logspace(-1, 0, steps)
        params = {'bandwidth': bwspace}

        # removed parameter "iid=False" because it is deprecated since scikit-
        # learn 0.22; however in older versions it was true by default so be
        # cautious
        grid = GridSearchCV(estimator, params, cv=5)
        grid.fit(data)
        return grid.best_estimator_

    @staticmethod
    def silverman_bw(data):
        """Calculate kernel bandwidth using Silverman's rule-of-thumb.

        Parameters
        ----------
        data : iterable of float
            input data

        Returns
        -------
        float
            bandwidth

        Raises
        ------
        ValueError
            If data size < 2.

        Notes
        -----
        bw = 0.9 * min(std, IQR / 1.34) * n ^ (-1/5)

        .. _Wikipedia:
            https://en.wikipedia.org/wiki/Kernel_density_estimation
        """
        n = len(data)
        if n < 2:
            raise ValueError(f'Cannot calculate bandwidth on {n} data point.')
        iqr = np.subtract(*np.percentile(data, [75, 25]))
        std = np.std(data, ddof=1)
        if not std and not iqr:
            bw = 1.0
        elif not std:
            bw = iqr / 1.34
        elif not iqr:
            bw = std
        elif std <= iqr / 1.34:
            bw = std
        else:
            bw = iqr / 1.34
        return bw * 0.9 * len(data) ** -0.2

    @staticmethod
    def density_func(data, kde, num=10000):
        """Get density function of KDE within given range.

        Parameters
        ----------
        data : np.array
            sample data (1D)
        kde : sklearn.neighbors.KernelDensity
            kernel density estimator
        num : int, optional
            number of data points in density function

        Returns
        -------
        np.array, np.array
            x and y values
        """
        x = np.linspace(data.min(), data.max(), num=num)
        y = np.exp(kde.score_samples(x[:, np.newaxis]))
        return x, y

    @staticmethod
    def first_hill(x, y):
        """Find first peak and first valley of density function.

        Parameters
        ----------
        x, y : np.array
            x and y values

        Returns
        -------
        float, float
            peak and valley

        Raises
        ------
        ValueError
            Cannot identify at least two peaks.
        ValueError
            Cannot identify at least one valley (unlikely).
        ValueError
            Peak is larger than valley.
        """
        # find peaks
        peaks = find_peaks(y)[0]
        if peaks.size < 2:
            raise ValueError('Cannot identify at least two peaks.')

        # find valleys
        valleys = find_peaks(np.negative(y))[0]
        if not valleys.size:
            raise ValueError('Cannot identify at least one valley.')

        # get first peak and first valley
        peak, valley = peaks[0], valleys[0]
        if peak > valley:
            raise ValueError('Peak is larger than valley.')

        # convert from index to actual value
        return x[peak], x[valley]

    @staticmethod
    def plot_hist(data, file):
        """Plot histogram.

        Parameters
        ----------
        data : iterable of float
            data to plot histogram
        file : str
            filename to save plot
        """
        fig = plt.figure(figsize=(5, 5))
        plt.hist(data)
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        save_figure(fig, file)
        plt.close()

    @staticmethod
    def plot_density(x, y, peak, valley, th, file):
        """Plot kernel density function.

        Parameters
        ----------
        x, y : np.array
            data points of density function
        peak, valley : float
            x-coordinate of first peak or valley
        th : float
            x-coordinate of threshold
        file : str
            filename to save plot

        Notes
        -----
        1. Generate a density plot, with 1st peak and 1st valley circled, and
           threshold indicated by a vertical dashed line.
        2. Export density plot to image file `file`.
        """
        fig = plt.figure(figsize=(5, 5))
        plt.plot(x, y)
        plt.plot([peak], [y[np.where(x == peak)]], marker='o')
        plt.plot([valley], [y[np.where(x == valley)]], marker='o')
        plt.axvline(th, color='grey', linestyle='--')
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        save_figure(fig, file)
        plt.close()

    def smart_kde(self, group):
        """Automatically determine kernel bandwidth for the goal of this
        analysis (HGT prediction).

        Parameters
        ----------
        group : str
            which group to cluster

        Returns
        -------
        float
            threshold, or 0 if unable to determine
        """
        data = self.df[group].values
        scaler = StandardScaler()
        data_ = scaler.fit_transform(data[:, np.newaxis])
        estimator = KernelDensity(kernel='gaussian')
        bwspace = np.logspace(0, -1, self.bw_steps)

        # move through bandwidth space from high to low, until an atypical hill
        # at low end is no larger than threshold
        for bw in bwspace:
            setattr(estimator, 'bandwidth', bw)
            kde = estimator.fit(data_)
            x, y = self.density_func(data_, kde)
            x = scaler.inverse_transform(x)
            y = scaler.inverse_transform(y)
            try:
                peak, valley = self.first_hill(x, y)
            except ValueError:
                print(f'  {bw:.3f}: n/a')
                continue
            th = valley - (valley - peak) * self.noise / 100
            ratio = data[data < th].size / data.size * 100
            print(f'  {bw:.3f}: {th:g} - {ratio:.2f}%')
            if not self.low_part or ratio <= self.low_part:
                print(f'  Auto-determined bandwidth: {bw:g}.')
                self.plot_density(x, y, peak, valley, th,
                                  join(self.output, f'{group}.kde.png'))
                return th
        return 0.0

    def calc_cluster_props(self):
        """Calculate cluster properties.

        Returns
        -------
        tuple of float
            centroid of prediction results

        Notes
        -----
        Add column `silh` to score table as silhouette scores.
        """
        data = self.df[self.relevant_groups()]
        scaler = StandardScaler()
        data_ = scaler.fit_transform(data)
        labels = self.df['hgt'].tolist()

        # calculate silhouette scores for all samples
        self.df['silh'] = silhouette_samples(data_, labels)

        # calculate centroid
        clf = NearestCentroid()
        clf.fit(data_, self.df['hgt'])
        cent = scaler.inverse_transform(clf.centroids_[1])
        return cent

    def refine_cluster(self, cent):
        """Refine cluster.

        Parameters
        ----------
        cent : tuple of float
            centroid

        Notes
        -----
        1. Append boolean column `far` to score table indicating whether data
           points fall beyond cluster boundaries.
        2. Refine predictions based on 1) position relative to centroid and 2)
           whether beyond cluster boundaries.
        """
        # protect data on far side of centroid in all groups
        if self.self_low:
            self.df['far'] = (self.df['self'] < cent[0]) &\
                             (self.df['close'] < cent[1]) &\
                             (self.df['distal'] > cent[2])
        else:
            self.df['far'] = (self.df['close'] < cent[0]) &\
                             (self.df['distal'] > cent[1])

        # drop data that are low in silhouette score and on near side
        # of centroid
        self.df['hgt'] = self.df['hgt'] & (
            self.df['far'] | (self.df['silh'] >= self.silhouette))

    def plot_hgts(self):
        """Plot HGT prediction results.

        Notes
        -----
        1. Generate a scatter plot with putative HGTs colored.
        2. Export scatter plot to image file `scatter.png`.
        """
        fig = plt.figure(figsize=(5, 5))
        plt.scatter('close', 'distal', c='hgt', data=self.df)
        plt.xlabel('Close')
        plt.ylabel('Distal')
        save_figure(fig, join(self.output, 'scatter.png'))
        plt.close()

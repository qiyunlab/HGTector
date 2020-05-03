#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove, makedirs
from os.path import join, isdir, isfile, dirname, realpath
from shutil import rmtree, copy, move
from tempfile import mkdtemp

import numpy as np
import pandas as pd

from sklearn.neighbors import KernelDensity

from hgtector.analyze import Analyze
from hgtector.util import (
    load_configs, add_children, get_descendants, taxdump_from_text)


class AnalyzeTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

        np.random.seed(42)
        self.dist_norm1 = np.random.normal(5.0, 2.0, 1500)
        self.dist_norm2 = np.random.normal(1.0, 0.5, 500)
        self.dist_lognorm = np.random.lognormal(0.0, 1.0, 1000)
        self.dist_gamma = np.random.gamma(2, 2, 800)

    def tearDown(self):
        rmtree(self.tmpdir)

    def test___call__(self):
        # run Ecoli sample using the Silverman method
        me = Analyze()
        def args(): return None
        args.input = join(self.datadir, 'Ecoli', 'search')
        args.output = join(self.tmpdir, 'output')
        args.taxdump = join(self.datadir, 'Ecoli', 'taxdump')
        args.input_tax = None
        args.self_tax = None
        args.close_tax = None
        args.self_rank = None
        args.close_size = None
        args.distal_top = None
        args.bandwidth = 'silverman'
        args.from_scores = False
        me(args)
        self.assertEqual(me.df[me.df['hgt']].shape[0], 16)

        # use existing score table, run grid search
        args.input = None
        args.from_scores = True
        args.bandwidth = 'grid'
        me(args)
        self.assertEqual(me.df[me.df['hgt']].shape[0], 18)
        rmtree(args.output)

    def test_set_parameters(self):
        me = Analyze()
        me.cfg = load_configs()
        def args(): return None

        # input is file
        infile = join(self.datadir, 'DnaK', 'search', 'sample.tsv')
        outdir = join(self.tmpdir, 'output')
        args.input = infile
        args.output = outdir
        args.noise = 0.75
        me.set_parameters(args)
        self.assertEqual(me.input, infile)
        self.assertEqual(me.output, outdir)
        self.assertTrue(isdir(outdir))
        self.assertDictEqual(me.input_map, {'sample': infile})
        self.assertEqual(me.noise, 75)

        # coverage threshold too small
        args.input_cov = 25
        with self.assertRaises(ValueError) as ctx:
            me.set_parameters(args)
        msg = 'Taxonomy coverage for auto-interence must be at least 50%.'
        self.assertEqual(str(ctx.exception), msg)
        args.input_cov = 75

        # input is directory
        indir = join(self.datadir, 'DnaK', 'search')
        args.input = indir
        me.set_parameters(args)
        self.assertEqual(me.input, indir)
        self.assertDictEqual(me.input_map, {'sample': infile})
        rmtree(outdir)

        # input is invalid
        not_path = 'I am not a path'
        args.input = not_path
        with self.assertRaises(ValueError) as ctx:
            me.set_parameters(args)
        msg = f'Invalid input data file or directory: {not_path}.'
        self.assertEqual(str(ctx.exception), msg)

        # input has no search result
        args.input = self.tmpdir
        with self.assertRaises(ValueError) as ctx:
            me.set_parameters(args)
        msg = f'No input data are found under: {self.tmpdir}.'
        self.assertEqual(str(ctx.exception), msg)

        # no input (which is okay)
        delattr(me, 'input_map')
        args.input = None
        me.set_parameters(args)
        self.assertFalse(hasattr(me, 'input_map'))

    def test_read_input(self):
        me = Analyze()

        def batch_assert():
            self.assertEqual(len(me.taxdump), 76)
            self.assertEqual(me.data['sample'][0]['id'], 'WP_000516135.1')
            self.assertEqual(me.data['sample'][0]['hits'].shape, (12, 5))

        # DnaK - default mode
        me.taxdump = join(self.datadir, 'DnaK', 'taxdump')
        me.input_map = {'sample': join(
            self.datadir, 'DnaK', 'search', 'sample.tsv')}
        me.read_input()
        batch_assert()

        # missing taxonomy
        copy(join(self.datadir, 'DnaK', 'search', 'sample.tsv'),
             join(self.tmpdir, 'sample.tsv'))
        me.input = self.tmpdir
        me.taxdump = None
        with self.assertRaises(ValueError) as ctx:
            me.read_input()
        msg = 'Missing taxonomy database.'
        self.assertEqual(str(ctx.exception), msg)

        # taxonomy in same directory as search result
        copy(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
             join(self.tmpdir, 'nodes.dmp'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
             join(self.tmpdir, 'names.dmp'))
        me.input_map = {'sample': join(self.tmpdir, 'sample.tsv')}
        me.read_input()
        batch_assert()

        # taxonomy in parent directory as search result
        indir = join(self.tmpdir, 'search')
        makedirs(indir)
        move(join(self.tmpdir, 'sample.tsv'), join(indir, 'sample.tsv'))
        me.input = indir
        me.input_map = {'sample': join(indir, 'sample.tsv')}
        me.taxdump = None
        me.read_input()
        batch_assert()
        rmtree(indir)
        remove(join(self.tmpdir, 'nodes.dmp'))
        remove(join(self.tmpdir, 'names.dmp'))

    def test_read_search_results(self):
        file = join(self.datadir, 'DnaK', 'search', 'sample.tsv')
        obs = Analyze.read_search_results(file)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0]['id'], 'WP_000516135.1')
        self.assertAlmostEqual(obs[0]['score'], 1092.8)
        self.assertTupleEqual(obs[0]['hits'].shape, (12, 5))
        self.assertEqual(obs[0]['hits'].iloc[2].name, 'NP_454622.1')
        self.assertAlmostEqual(obs[0]['hits']['evalue']['NP_230502.1'],
                               5.9e-282)
        self.assertEqual(obs[0]['hits']['taxid']['NP_384288.1'], '266834')

        # maximum number of hits
        obs = Analyze.read_search_results(file, 5)
        self.assertEqual(len(obs[0]['hits']), 5)

    def test_assign_taxonomy(self):
        # input are two genomes with defined taxonomy
        me = Analyze()
        me.input_tax = 'S1:561,S2:620'  # Escherichia and Shigella
        me.data = {}
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.assign_taxonomy()
        # test input taxonomy extraction
        self.assertDictEqual(me.input_tax, {'S1': '561', 'S2': '620'})
        # test taxonomy refinement
        exp = {'1', '131567', '2', '1224', '1236', '91347', '543', '561',
               '620'}
        self.assertSetEqual(set(me.taxdump.keys()), exp)
        # test LCA discovery
        self.assertEqual(me.lca, '543')

        # helper for making hit table
        def _hits_df(d):
            return pd.Series(d, name='taxid', dtype=object).to_frame()

        # input is one genome with defined taxonomy
        me = Analyze()
        me.data = {'S1': [{'hits': pd.DataFrame(columns=['taxid'])}]}
        me.input_tax = '561'  # Escherichia
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.assign_taxonomy()
        self.assertDictEqual(me.input_tax, {'S1': '561'})

        # input taxonomy not found in database
        me.input_tax = '1234'
        me.taxdump = taxdump_from_text(taxdump_proteo)
        with self.assertRaises(ValueError) as ctx:
            me.assign_taxonomy()
        msg = 'TaxID 1234 is not present in taxonomy database.'
        self.assertEqual(str(ctx.exception), msg)

        # input are two genome whose taxonomies are to be inferred based on
        # search results
        me = Analyze()
        me.input_tax = None
        me.data = {'S1': [{'hits': _hits_df({'P1': '561', 'P2': '562'})},
                          {'hits': _hits_df({'P3': '543', 'P4': '561'})}],
                   'S2': [{'hits': _hits_df({'P5': '562', 'P6': '585056'})},
                          {'hits': _hits_df({'P7': '561', 'P8': '1038927'})},
                          {'hits': _hits_df({'P9': '2580236'})}]}
        me.input_cov = 75
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.assign_taxonomy()
        self.assertDictEqual(me.input_tax, {'S1': '543', 'S2': '561'})
        self.assertEqual(me.lca, '543')

        # cannot auto-infer taxonomy
        me.data['S3'] = [{'hits': _hits_df({})}]
        me.taxdump = taxdump_from_text(taxdump_proteo)
        with self.assertRaises(ValueError) as ctx:
            me.assign_taxonomy()
        msg = 'Cannot auto-infer taxonomy for S3. Please specify manually.'
        self.assertEqual(str(ctx.exception), msg)

        # invalid input taxonomy string
        me.input_tax = '561'
        with self.assertRaises(ValueError) as ctx:
            me.assign_taxonomy()
        msg = 'Invalid input taxonomy format.'
        self.assertEqual(str(ctx.exception), msg)

    def test_infer_genome_tax(self):
        taxdump = taxdump_from_text(taxdump_proteo)

        # five proteins, in which four have hits
        taxids = [['562', '620', '570'],  # E. coli
                  ['562', '585056', '1038927', '2'],  # E. coli
                  ['561', '543', '776'],  # Escherichia
                  ['548', '570', '1236'],  # K. aerogenes
                  []]
        prots = [{'hits': pd.DataFrame(x, columns=['taxid'])} for x in taxids]
        obs = Analyze.infer_genome_tax(prots, taxdump, 75)
        exp = ('561', 75.0)  # 3 / 4 best hits assigned to Escherichia
        self.assertTupleEqual(obs, exp)

        # reduce coverage threshold
        obs = Analyze.infer_genome_tax(prots, taxdump, 50)
        exp = ('562', 50.0)  # 2 / 4 best hits assigned to Escherichia
        self.assertTupleEqual(obs, exp)

        # remove one protein that best matches E. coli
        prots.pop(0)
        obs = Analyze.infer_genome_tax(prots, taxdump, 75)
        exp = ('543', 100.0)  # 3 / 3 best hits assigned to Enterobacteriaceae
        self.assertTupleEqual(obs, exp)

        # no input protein
        with self.assertRaises(ValueError) as ctx:
            Analyze.infer_genome_tax({}, taxdump, 75)
        msg = 'Cannot auto-infer taxonomy.'
        self.assertEqual(str(ctx.exception), msg)

    def test_sum_taxids(self):
        me = Analyze()
        me.input_tax = {'S1': '1', 'S2': '3'}

        def _hits_df(d):
            return pd.Series(d, name='taxid').to_frame()

        me.data = {'S1': [{'hits': _hits_df({'a': '4', 'b': '6'})},
                          {'hits': _hits_df({'a': '4', 'c': '8'})}],
                   'S2': [{'hits': _hits_df({'b': '6', 'd': '1'})}]}
        obs = me.sum_taxids()
        exp = {'1', '3', '4', '6', '8'}
        self.assertSetEqual(obs, exp)

    def test_define_groups(self):
        me = Analyze()
        me.taxdump = taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)
        me.groups = {}

        # user defined groups:
        # self: genera Escherichia and Shigella
        # close: family Enterobacteriaceae
        me.groups = {}
        me.self_tax = '561,620'
        me.close_tax = '543'
        me.define_groups()
        self.assertListEqual(me.self_tax, ['561', '620'])
        exp = {'561', '562', '585056', '1038927', '2580236', '620', '622'}
        self.assertSetEqual(me.groups['self'], exp)
        self.assertListEqual(me.close_tax, ['543'])
        exp = {'543', '548', '570'}
        self.assertSetEqual(me.groups['close'], exp)

        # auto-infer groups
        me.self_tax = {}
        me.close_tax = {}
        me.lca = '562'          # all inputs are E. coli
        me.self_rank = 'genus'  # but we want to raise self to genus
        me.close_size = 2       # close group must be this big or bigger
        me.define_groups()
        self.assertListEqual(me.self_tax, ['561'])
        exp = {'561', '562', '585056', '1038927', '2580236'}
        self.assertSetEqual(me.groups['self'], exp)
        self.assertListEqual(me.close_tax, ['543'])
        exp = {'543', '548', '570', '620', '622'}
        self.assertSetEqual(me.groups['close'], exp)

    def test_infer_self_group(self):
        me = Analyze()
        me.taxdump = taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)

        # assign to LCA of all genomes (E. coli)
        me.self_tax = None
        me.lca = '562'
        me.self_rank = None
        me.infer_self_group()
        self.assertListEqual(me.self_tax, ['562'])

        # raise LCA to genus level (Escherichia)
        me.self_tax = None
        me.lca = '562'
        me.self_rank = 'genus'
        me.infer_self_group()
        self.assertListEqual(me.self_tax, ['561'])

        # LCA (Enterobacteriaceae) is already above designated rank (genus)
        me.self_tax = None
        me.lca = '543'
        me.self_rank = 'genus'
        me.infer_self_group()
        self.assertListEqual(me.self_tax, ['543'])

    def test_infer_close_group(self):
        me = Analyze()
        me.taxdump = taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)
        me.groups = {}

        # close group is parent of LCA of self group
        me.self_tax = ['562']  # E. coli
        me.groups['self'] = set(['562'] + get_descendants('562', me.taxdump))
        me.close_tax = None
        me.close_size = None
        me.infer_close_group()
        self.assertListEqual(me.close_tax, ['561'])  # Escherichia
        self.assertSetEqual(me.groups['close'], {'561', '2580236'})

        # close group must have at least 5 taxa
        me.close_tax = None
        me.groups['close'] = None
        me.close_size = 5
        me.infer_close_group()
        self.assertListEqual(me.close_tax, ['543'])  # Enterobacteriaceae
        exp = {'543', '620', '622', '570', '548', '561', '2580236'}
        self.assertSetEqual(me.groups['close'], exp)

        # close group is LCA of multiple self groups
        me.self_tax = ['561', '620']  # Escherichia and Shigella
        me.groups['self'] = set().union(*[[x] + get_descendants(
            x, me.taxdump) for x in me.self_tax])
        me.close_tax = None
        me.groups['close'] = None
        me.close_size = None
        me.infer_close_group()
        self.assertListEqual(me.close_tax, ['543'])  # Enterobacteriaceae
        exp = {'543', '570', '548'}
        self.assertSetEqual(me.groups['close'], exp)

    def test_calc_scores(self):
        columns = ('id', 'taxid', 'score')

        # helper for making hit table
        def _hits_df(data):
            return pd.DataFrame(data, columns=columns).set_index('id')

        me = Analyze()
        me.taxdump = taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)
        me.groups = {'self':  {'561', '562', '585056'},
                     'close': {'543', '91347', '1236'}}
        me.data = {'S1': [
            {'score': 100, 'hits': _hits_df((('P1', '561', 100),
                                             ('P2', '562', 95)))},
            {'score': 90,  'hits': _hits_df((('P3', '561', 81),
                                             ('P4', '543', 72)))}],
                   'S2': [
            {'score': 96,  'hits': _hits_df((('P5', '561', 90),
                                             ('P6', '543', 84),
                                             ('P7', '620', 66)))}]}
        me.weighted = True
        me.match_th = 0.9
        me.calc_scores()

        # helper for get scores
        def _prot_scores(prot):
            return [prot[x] for x in ('self', 'close', 'distal')]

        s1_1 = me.data['S1'][0]
        self.assertListEqual(s1_1['hits']['group'].tolist(), ['self', 'self'])
        self.assertListEqual(_prot_scores(s1_1), [1.95, 0.0, 0.0])
        self.assertEqual(s1_1['match'], '0')
        s1_2 = me.data['S1'][1]
        self.assertListEqual(s1_2['hits']['group'].tolist(), ['self', 'close'])
        self.assertListEqual(_prot_scores(s1_2), [0.9, 0.8, 0.0])
        self.assertEqual(s1_2['match'], '0')
        s2_1 = me.data['S2'][0]
        self.assertListEqual(s2_1['hits']['group'].tolist(),
                             ['self', 'close', 'distal'])
        self.assertListEqual(_prot_scores(s2_1), [0.9375, 0.875, 0.6875])
        self.assertEqual(s2_1['match'], '620')

    def test_find_match(self):
        me = Analyze()
        me.taxdump = taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)
        df = pd.DataFrame(
            [[100, '585056'],  # E. coli UMN026
             [99, '1038927'],  # E. coli O104:H4
             [97, '562'],      # Escherichia coli
             [95, '622'],      # Shigella dysenteriae
             [92, '543'],      # Enterobacteriaceae
             [88, '548'],      # Klebsiella aerogenes
             [80, '766']],     # Rickettsiales
            columns=['score', 'taxid'])

        # keep top 1% hits
        me.match_th = 0.99
        self.assertEqual(me.find_match(df), '562')

        # keep top 10% hits
        me.match_th = 0.9
        self.assertEqual(me.find_match(df), '543')

        # keep top 20% hits
        me.match_th = 0.8
        self.assertEqual(me.find_match(df), '1224')

        # input DataFrame is empty
        self.assertEqual(me.find_match(pd.DataFrame()), '0')

    def test_make_score_table(self):
        me = Analyze()
        me.output = self.tmpdir
        me.data = {'S1': [{'id': 'P1', 'length': 100, 'match': '0',
                           'self': 1.5, 'close': 0.75, 'distal': 0.0,
                           'hits': pd.DataFrame([0] * 3)},
                          {'id': 'P2', 'length': 120, 'match': '1224',
                           'self': 1.625, 'close': 0.225, 'distal': 0.375,
                           'hits': pd.DataFrame([0] * 5)}],
                   'S2': [{'id': 'P1', 'length': 225, 'match': '620',
                           'self': 2.35, 'close': 1.05, 'distal': 0.75,
                           'hits': pd.DataFrame([0] * 6)}]}
        me.make_score_table()
        obs = me.df.values.tolist()
        exp = [['S1', 'P1', 100, 3, 1.5,   0.75,  0,        '0'],
               ['S1', 'P2', 120, 5, 1.625, 0.225, 0.375, '1224'],
               ['S2', 'P1', 225, 6, 2.35,  1.05,  0.75,   '620']]
        self.assertListEqual(obs, exp)
        fp = join(self.tmpdir, 'scores.tsv')
        with open(fp, 'r') as f:
            obs = [x.split('\t') for x in f.read().splitlines()[1:]]
        exp = [[str(y) for y in x] for x in exp]
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_remove_orphans(self):
        me = Analyze()
        me.df = pd.DataFrame([
            [1.0, 0.2], [0.5, 0.4], [0.0, 0.0], [0.8, 0.0], [0.0, 0.7]],
            columns=['close', 'distal'])
        me.remove_orphans()
        self.assertListEqual(me.df.values.tolist(), [
            [1.0, 0.2], [0.5, 0.4], [0.8, 0.0], [0.0, 0.7]])

    def test_remove_outliers(self):
        me = Analyze()
        me.self_low = False
        df = pd.DataFrame(np.array([self.dist_gamma,
                                    self.dist_lognorm[:800]]).T,
                          columns=['close', 'distal'])

        # Z-score
        me.df = df.copy()
        me.outliers = 'zscore'
        me.remove_outliers()
        self.assertEqual(me.df.shape[0], 781)

        # boxplot
        me.df = df.copy()
        me.outliers = 'boxplot'
        me.remove_outliers()
        self.assertEqual(me.df.shape[0], 710)

    def test_relevant_groups(self):
        me = Analyze()
        me.self_low = False
        self.assertListEqual(me.relevant_groups(), ['close', 'distal'])
        me.self_low = True
        self.assertListEqual(me.relevant_groups(), ['self', 'close', 'distal'])

    def test_outliers_zscore(self):
        df = pd.DataFrame(np.array([self.dist_gamma,
                                    self.dist_lognorm[:800]]).T,
                          columns=['close', 'distal'])
        obs = Analyze.outliers_zscore(df, ['close', 'distal'])
        self.assertEqual(obs.shape[0], 781)

    def test_outliers_boxplot(self):
        df = pd.DataFrame(np.array([self.dist_gamma,
                                    self.dist_lognorm[:800]]).T,
                          columns=['close', 'distal'])
        obs = Analyze.outliers_boxplot(df, ['close', 'distal'])
        self.assertEqual(obs.shape[0], 710)

    def test_predict_hgt(self):
        me = Analyze()

        # populate score table
        n = 1000
        data = {'sample': ['S1'] * n,
                'protein': [f'P{x}' for x in range(n)],
                'self': np.random.choice(self.dist_gamma, n),
                'close': np.concatenate((
                    np.random.choice(self.dist_norm1, int(n / 2)) / 3,
                    np.random.choice(self.dist_norm2, int(n / 2)))),
                'distal': np.concatenate((
                    np.random.choice(self.dist_lognorm, int(n * 3 / 4)),
                    np.random.choice(self.dist_gamma, int(n / 4)) / 2))}
        me.df = pd.DataFrame(data)

        # default setting
        me.output = self.tmpdir
        me.self_low = False
        me.bandwidth = 'auto'
        me.bw_steps = 20
        me.low_part = 75
        me.fixed = 25
        me.noise = 50
        me.silhouette = 0.5

        # run prediction
        self.assertEqual(me.predict_hgt(), 96)
        groups = ['self', 'close', 'distal']
        for group in groups[1:]:
            fp = join(self.tmpdir, f'{group}.hist.png')
            self.assertTrue(isfile(fp))
            remove(fp)
        fp = join(self.tmpdir, 'scatter.png')
        self.assertTrue(isfile(fp))
        remove(fp)
        fp = join(self.tmpdir, 'hgts')
        self.assertTrue(isfile(join(fp, 'S1.txt')))
        rmtree(fp)

        # constant values
        me.df['close'] = 1
        me.df.drop('hgt', axis=1, inplace=True)
        self.assertEqual(me.predict_hgt(), 0)
        self.assertNotIn('hgt', me.df.columns)
        remove(join(self.tmpdir, 'close.hist.png'))

    def test_cluster_kde(self):
        me = Analyze()
        data = np.concatenate([self.dist_norm1, self.dist_norm2])
        me.df = pd.Series(data, name='group').to_frame()
        me.bw_steps = 10
        me.noise = 50
        me.low_part = 75
        me.output = self.tmpdir

        # grid search
        me.bandwidth = 'grid'
        obs = me.cluster_kde('group')
        self.assertAlmostEqual(obs, 1.855525575742988)

        # Silverman's rule-of-thumb
        me.bandwidth = 'silverman'
        obs = me.cluster_kde('group')
        self.assertAlmostEqual(obs, 2.2279977615745703)

        # fixed value
        me.bandwidth = 0.5
        obs = me.cluster_kde('group')
        self.assertAlmostEqual(obs, 2.2507008281395433)

        # smart KDE
        me.bandwidth = 'auto'
        obs = me.cluster_kde('group')
        self.assertAlmostEqual(obs, 2.1903958075763343)

        # clean up
        remove(join(self.tmpdir, 'group.kde.png'))

        # cannot find threshold (unimodal distribution)
        me.df = pd.Series(self.dist_norm1, name='group').to_frame()
        me.bandwidth = 'silverman'
        obs = me.cluster_kde('group')
        self.assertEqual(obs, 0)

    def test_perform_kde(self):
        me = Analyze()
        me.bw_steps = 10
        data = np.concatenate([self.dist_norm1, self.dist_norm2])

        # grid search
        me.bandwidth = 'grid'
        obs = me.perform_kde(data)[2]
        self.assertAlmostEqual(obs, 0.21544346900318834)

        # Silverman's rule-of-thumb
        me.bandwidth = 'silverman'
        obs = me.perform_kde(data)[2]
        self.assertAlmostEqual(obs, 0.48713295460585126)

        # fixed value
        me.bandwidth = 0.5
        obs = me.perform_kde(data)[2]
        self.assertAlmostEqual(obs, 0.5)

        # invalid bandwidth
        me.bandwidth = 100
        with self.assertRaises(ValueError) as ctx:
            me.perform_kde(data)
        msg = 'Invalid bandwidth: 100.'
        self.assertEqual(str(ctx.exception), msg)

    def test_grid_kde(self):
        estimator = KernelDensity(kernel='gaussian')

        # unimodal
        data = self.dist_gamma[:, np.newaxis]
        obs = Analyze.grid_kde(data, estimator, 10).bandwidth
        self.assertAlmostEqual(obs, 0.774263682681127)

        # bimodal
        data = np.concatenate([
            self.dist_norm1, self.dist_norm2])[:, np.newaxis]
        obs = Analyze.grid_kde(data, estimator, 10).bandwidth
        self.assertAlmostEqual(obs, 0.46415888336127786)

        data = np.array([1, 2, 3, 4, 5])[:, np.newaxis]
        obs = Analyze.grid_kde(data, estimator, 5).bandwidth
        self.assertAlmostEqual(obs, 1.0)

        # very few data points (bw = high end)
        data = np.array([1, 2, 3, 4, 5])[:, np.newaxis]
        obs = Analyze.grid_kde(data, estimator, 5).bandwidth
        self.assertAlmostEqual(obs, 1.0)

        # constant values (bw = low end)
        data = np.array([1, 1, 1, 1, 1])[:, np.newaxis]
        obs = Analyze.grid_kde(data, estimator, 5).bandwidth
        self.assertAlmostEqual(obs, 0.1)

        # too few data points (less than splits)
        data = np.array([1, 2, 3])[:, np.newaxis]
        with self.assertRaises(ValueError) as ctx:
            Analyze.grid_kde(data, estimator, 5)
        msg = 'Cannot perform grid search on 3 data point(s).'
        self.assertEqual(str(ctx.exception), msg)

    def test_silverman_bw(self):
        # unimodal
        obs = Analyze.silverman_bw(self.dist_gamma)
        self.assertAlmostEqual(obs, 0.6148288686346546)
        obs = Analyze.silverman_bw(self.dist_lognorm)
        self.assertAlmostEqual(obs, 0.2384666552244172)

        # bimodal
        obs = Analyze.silverman_bw(np.concatenate([
            self.dist_norm1, self.dist_norm2]))
        self.assertAlmostEqual(obs, 0.48713295460585126)

        # constant values
        obs = Analyze.silverman_bw([1, 1, 1, 1, 1])
        self.assertAlmostEqual(obs, 0.652301697309926)

        # IQR = 0
        obs = Analyze.silverman_bw([1, 3, 3, 3, 5])
        self.assertAlmostEqual(obs, 0.9224939070946869)

        # one element
        with self.assertRaises(ValueError) as ctx:
            Analyze.silverman_bw([5])
        msg = 'Cannot calculate bandwidth on 1 data point.'
        self.assertEqual(str(ctx.exception), msg)

    def test_density_func(self):
        data = self.dist_norm1[:, np.newaxis]
        estimator = KernelDensity(kernel='gaussian', bandwidth=0.5)
        kde = estimator.fit(data)
        obs = Analyze.density_func(data, kde, 10)
        exp = (np.array([-1.48253468, 0.0939095, 1.67035369, 3.24679787,
                         4.82324206, 6.39968624, 7.97613043, 9.55257461,
                         11.1290188, 12.70546298]),
               np.array([0.00104342, 0.00788705, 0.0496806, 0.13173376,
                         0.19176352, 0.15754466, 0.06992292, 0.02140856,
                         0.00150463, 0.00053637]))
        np.testing.assert_array_almost_equal(obs, exp)

    def test_first_hill(self):
        # typical bimodal distribution
        data = np.concatenate([
            self.dist_norm1, self.dist_norm2])[:, np.newaxis]
        estimator = KernelDensity(kernel='gaussian', bandwidth=0.5)
        kde = estimator.fit(data)
        x, y = Analyze.density_func(data, kde, 100)
        obs_x, obs_y = Analyze.first_hill(x, y)
        exp_x, exp_y = 1.0971012583068704, 2.5302323352207674
        self.assertAlmostEqual(obs_x, exp_x)
        self.assertAlmostEqual(obs_y, exp_y)

        # peak larger than valley
        data = np.negative(data)
        kde = estimator.fit(data)
        x, y = Analyze.density_func(data, kde, 100)
        with self.assertRaises(ValueError) as ctx:
            Analyze.first_hill(x, y)
        msg = 'Peak is larger than valley.'
        self.assertEqual(str(ctx.exception), msg)

        # unimodal distribution
        data = self.dist_norm1[:, np.newaxis]
        kde = estimator.fit(data)
        x, y = Analyze.density_func(data, kde, 100)
        with self.assertRaises(ValueError) as ctx:
            Analyze.first_hill(x, y)
        msg = 'Cannot identify at least two peaks.'
        self.assertEqual(str(ctx.exception), msg)

    def test_plot_hist(self):
        fp = join(self.tmpdir, 'tmp.png')
        Analyze.plot_hist(self.dist_gamma, fp)
        self.assertTrue(isfile(fp))
        remove(fp)

    def test_plot_density(self):
        data = np.concatenate([
            self.dist_norm1, self.dist_norm2])[:, np.newaxis]
        estimator = KernelDensity(kernel='gaussian', bandwidth=0.5)
        kde = estimator.fit(data)
        x, y = Analyze.density_func(data, kde, 100)
        peak, valley = Analyze.first_hill(x, y)
        th = valley - (valley - peak) * 0.5 / 100
        fp = join(self.tmpdir, 'tmp.png')
        Analyze.plot_density(x, y, peak, valley, th, fp)
        self.assertTrue(isfile(fp))
        remove(fp)

    def test_smart_kde(self):
        me = Analyze()

        # typical case (bimodal distribution)
        me.df = pd.Series(np.concatenate([
            self.dist_norm1, self.dist_norm2]), name='group').to_frame()
        me.bw_steps = 10
        me.noise = 50
        me.low_part = 75
        me.output = self.tmpdir
        obs = me.smart_kde('group')
        self.assertAlmostEqual(obs, 2.1903958075763343)
        file = join(self.tmpdir, 'group.kde.png')
        self.assertTrue(isfile(file))
        remove(file)

        # unable to determine threshold
        me.low_part = 0.001
        me.df = pd.Series(self.dist_norm1, name='group').to_frame()
        self.assertEqual(me.smart_kde('group'), 0)

    def test_calc_cluster_props(self):
        me = Analyze()
        me.self_low = False
        me.df = pd.DataFrame(np.array(
            [self.dist_gamma, self.dist_lognorm[:800]]).T,
            columns=['close', 'distal'])
        me.df['hgt'] = (me.df['close'] < 2) & (me.df['distal'] > 2)
        obs = me.calc_cluster_props()
        self.assertAlmostEqual(obs[0], 1.094658052928843)
        self.assertAlmostEqual(obs[1], 4.30076698399293)
        obs = me.df['silh'].describe()
        self.assertAlmostEqual(obs['mean'], 0.312495082044277)
        self.assertAlmostEqual(obs['std'], 0.21945541659155993)
        self.assertEqual(me.df.query('hgt & silh < 0.5').shape[0], 35)

    def test_refine_cluster(self):
        me = Analyze()

        # only close and distal
        me.self_low = False
        me.silhouette = 0.5
        me.df = pd.DataFrame(np.array(
            [self.dist_gamma, self.dist_lognorm[:800]]).T,
            columns=['close', 'distal'])
        me.df['hgt'] = (me.df['close'] < 2) & (me.df['distal'] > 2)
        me.refine_cluster(me.calc_cluster_props())
        self.assertEqual(me.df[me.df['hgt']].shape[0], 11)

        # all three groups
        me.self_low = True
        me.df = pd.DataFrame(np.array([
            self.dist_norm1[:800], self.dist_gamma,
            self.dist_lognorm[:800]]).T,
            columns=['self', 'close', 'distal'])
        me.df['hgt'] = (me.df['close'] < 2) & (me.df['distal'] > 2)
        me.refine_cluster(me.calc_cluster_props())
        self.assertEqual(me.df[me.df['hgt']].shape[0], 4)

    def test_plot_hgts(self):
        me = Analyze()
        me.output = self.tmpdir
        me.df = pd.DataFrame(np.array(
            [self.dist_gamma, self.dist_lognorm[:800]]).T,
            columns=['close', 'distal'])
        me.df['hgt'] = (me.df['close'] < 2) & (me.df['distal'] > 2)
        me.plot_hgts()
        fp = join(self.tmpdir, 'scatter.png')
        self.assertTrue(isfile(fp))
        remove(fp)


"""Constants"""

taxdump_proteo = (
    '1,root,1,no rank',
    '131567,cellular organisms,1,no rank',
    '2,Bacteria,131567,superkingdom',
    '1224,Proteobacteria,2,phylum',
    '28211,Alphaproteobacteria,1224,class',
    '766,Rickettsiales,28211,order',
    '1236,Gammaproteobacteria,1224,class',
    '91347,Enterobacterales,1236,order',
    '543,Enterobacteriaceae,91347,family',
    '561,Escherichia,543,genus',
    '562,Escherichia coli,561,species',
    '585056,Escherichia coli UMN026,562,no rank',
    '1038927,Escherichia coli O104:H4,562,no rank',
    '2580236,synthetic Escherichia coli Syn61,561,species',
    '620,Shigella,543,genus',
    '622,Shigella dysenteriae,620,species',
    '570,Klebsiella,543,genus',
    '548,Klebsiella aerogenes,570,species',
    '118884,unclassified Gammaproteobacteria,1236,no rank',
    '126792,Plasmid pPY113,1,species')


if __name__ == '__main__':
    main()

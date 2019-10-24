#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp

import pandas as pd

from hgtector.analyze import Analyze
from hgtector.util import add_children, get_descendants


class AnalyzeTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test___call__(self):
        # TODO
        pass

    def test_read_input(self):
        # TODO
        pass

    def test_assign_taxonomy(self):
        # TODO
        pass

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

    def test_infer_genome_tax(self):
        taxdump = _taxdump_from_text(taxdump_proteo)

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

    def test_define_groups(self):
        me = Analyze()
        me.taxdump = _taxdump_from_text(taxdump_proteo)
        add_children(me.taxdump)
        me.groups = {}

        # user defined groups:
        # self: genera Escherichia and Shigella
        # close: family Enterobacteriaceae
        me.self_tax = '561,620'
        me.close_tax = '543'
        me.define_groups()
        self.assertListEqual(me.self_tax, ['561', '620'])
        exp = {'561', '562', '585056', '1038927', '2580236', '620', '622'}
        self.assertSetEqual(me.groups['self'], exp)
        self.assertListEqual(me.close_tax, ['543'])
        exp = {'543', '548', '570'}
        self.assertSetEqual(me.groups['close'], exp)

    def test_infer_self_group(self):
        me = Analyze()
        me.taxdump = _taxdump_from_text(taxdump_proteo)
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
        me.taxdump = _taxdump_from_text(taxdump_proteo)
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
        # TODO
        pass

    def test_find_match(self):
        me = Analyze()
        me.taxdump = _taxdump_from_text(taxdump_proteo)
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
        # TODO
        pass


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


"""Helpers"""


def _taxdump_from_text(text):
    """Read taxdump from text.

    Parameters
    ----------
    text : list of str
        multi-line, tab-delimited text

    Returns
    -------
    dict of dict
        taxonomy database
    """
    res = {}
    for line in text:
        x = line.split(',')
        res[x[0]] = {'name': x[1], 'parent': x[2], 'rank': x[3]}
    return res


if __name__ == '__main__':
    main()

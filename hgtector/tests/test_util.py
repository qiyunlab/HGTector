#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove, getcwd, chdir
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp, mkstemp
import gzip

from hgtector.util import (
    timestamp, load_configs, find_config_file, get_config, arg2bool,
    read_file, list_from_param, dict_from_param, run_command, file2id,
    id2file_map, read_taxdump, read_prot2taxid, is_ancestral, taxid_at_rank,
    read_fasta, read_input_prots, contain_words, get_product, seqid2accver,
    taxids_at_ranks, find_lca, get_lineage, sort_by_hierarchy, refine_taxdump,
    _get_taxon, add_children, get_descendants, is_latin, is_capital,
    taxdump_from_text, rank_plural)


class UtilTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_timestamp(self):
        self.assertIsNotNone(getattr(timestamp(), 'timestamp'))

    def test_load_configs(self):
        obs = load_configs()
        self.assertIn('database', obs)

    def test_find_config_file(self):
        self.assertTrue(find_config_file().endswith('hgtector/config.yml'))
        cwd = getcwd()
        fp = join(self.tmpdir, 'config.yml')
        open(fp, 'w').close()
        chdir(self.tmpdir)
        self.assertEqual(find_config_file(), 'config.yml')
        chdir(cwd)
        remove(fp)

    def test_get_config(self):
        def obj(): return None

        # no content
        obj.cfg = None
        get_config(obj, 'test', 'test')
        self.assertFalse(hasattr(obj, 'test'))

        # one-level entry
        key = 'threshold'
        obj.cfg = {key: 10}
        get_config(obj, key, key)
        self.assertEqual(getattr(obj, key), 10)

        # two-level entry
        obj.cfg['search'] = {'evalue': 1e-5}
        get_config(obj, 'evalue', 'search.evalue')
        self.assertEqual(getattr(obj, 'evalue'), 1e-5)

        # three-level entry
        obj.cfg['search']['taxonomy'] = {'rank': 'genus'}
        get_config(obj, 'taxrank', 'search.taxonomy.rank')
        self.assertEqual(getattr(obj, 'taxrank'), 'genus')

        # attribute already has value
        obj.cfg['evalue'] = 1e-10
        get_config(obj, 'evalue', 'evalue')
        self.assertEqual(getattr(obj, 'evalue'), 1e-5)

        # attribute present but is None:
        obj.evalue = None
        get_config(obj, 'evalue', 'evalue')
        self.assertEqual(getattr(obj, 'evalue'), 1e-10)

        # level that does not exist
        get_config(obj, 'key', 'level.key')
        self.assertFalse(hasattr(obj, 'key'))

        # entry that does not exist
        get_config(obj, 'identity', 'search.identity')
        self.assertFalse(hasattr(obj, 'identity'))

        # entry that is None
        obj.cfg['search']['identity'] = None
        get_config(obj, 'identity', 'search.identity')
        self.assertFalse(hasattr(obj, 'identity'))

        # variable type conversion
        obj.cfg['coverage'] = '95'
        get_config(obj, 'coverage', 'coverage', float)
        self.assertEqual(getattr(obj, 'coverage'), 95)

        # string manipulation
        obj.cfg['rank'] = 'genus'
        get_config(obj, 'code', 'rank', lambda x: x[0])
        self.assertEqual(getattr(obj, 'code'), 'g')

    def test_arg2bool(self):
        # true
        self.assertTrue(arg2bool(True))
        self.assertTrue(arg2bool('True'))
        self.assertTrue(arg2bool('Yes'))
        self.assertTrue(arg2bool('y'))
        self.assertTrue(arg2bool('1'))

        # false
        self.assertFalse(arg2bool(None))
        self.assertFalse(arg2bool(False))
        self.assertFalse(arg2bool('False'))
        self.assertFalse(arg2bool('No'))
        self.assertFalse(arg2bool('n'))
        self.assertFalse(arg2bool('0'))

        # invalid input
        with self.assertRaises(ValueError) as ctx:
            arg2bool('test')
        msg = 'Boolean value expected.'
        self.assertEqual(str(ctx.exception), msg)

    def test_read_file(self):
        text = 'Hello World!'

        # read regular file
        fp = join(self.tmpdir, 'test.txt')
        with open(fp, 'w') as f:
            f.write(text)
        with read_file(fp) as f:
            obs = f.read()
        self.assertEqual(obs, text)

        # read compressed file
        fp = join(self.tmpdir, 'test.txt.gz')
        with gzip.open(fp, 'wb') as f:
            f.write(text.encode())
        with read_file(fp) as f:
            obs = f.read()
        self.assertEqual(obs, text)
        remove(fp)

    def test_run_command(self):
        # simple command
        cmd = 'echo "This is a test!"'
        obs = run_command(cmd)
        self.assertEqual(obs[0], 0)
        self.assertListEqual(obs[1], ['This is a test!'])

        # multiple-line output
        cmd = 'for i in 1 2 3; do echo $i; done'
        obs = run_command(cmd)[1]
        self.assertListEqual(obs, ['1', '2', '3'])

        # capture stderr
        cmd = 'echo "This is an error!" >&2'
        obs = run_command(cmd)[1]
        self.assertListEqual(obs, ['This is an error!'])

    def test_list_from_param(self):
        # nothing
        self.assertEqual(list_from_param(None), [])
        self.assertEqual(list_from_param(''), [])

        # already list
        self.assertListEqual(list_from_param([1, 2, 3]), [1, 2, 3])

        # list string
        self.assertListEqual(list_from_param('test'), ['test'])
        self.assertListEqual(list_from_param('a,b,c'), ['a', 'b', 'c'])

        # list file
        exp = ['this', 'is', 'a', 'list']
        fp = join(self.tmpdir, 'test.txt')
        with open(fp, 'w') as f:
            for e in exp:
                print(e, file=f)
        obs = list_from_param(fp)
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_dict_from_param(self):
        # nothing
        self.assertEqual(dict_from_param(None), {})
        self.assertEqual(dict_from_param(''), {})

        # already dict
        self.assertEqual(dict_from_param({'a': 1}), {'a': 1})

        # dict string
        exp = {'a': '1', 'b': '2', 'c': '3'}
        self.assertDictEqual(dict_from_param('a:1,b:2,c:3'), exp)

        # invalid dict string
        with self.assertRaises(ValueError) as ctx:
            dict_from_param('test')
        msg = 'Invalid dictionary string: "test".'
        self.assertEqual(str(ctx.exception), msg)

        # dict file
        fp = join(self.tmpdir, 'test.txt')
        with open(fp, 'w') as f:
            for itm in exp.items():
                print('\t'.join(itm), file=f)
        obs = dict_from_param(fp)
        self.assertDictEqual(obs, exp)

        # invalid dict file
        with open(fp, 'w') as f:
            for itm in exp.items():
                print('test', file=f)
        with self.assertRaises(ValueError) as ctx:
            dict_from_param(fp)
        msg = f'Invalid dictionary file: "{fp}".'
        self.assertEqual(str(ctx.exception), msg)
        remove(fp)

    def test_file2id(self):
        self.assertEqual(file2id('sample.fna'), 'sample')
        self.assertEqual(file2id('sample.fna.gz'), 'sample')
        self.assertEqual(file2id('/path/to/sample.fna'), 'sample')
        self.assertEqual(file2id('sample'), 'sample')
        self.assertEqual(file2id(''), '')

    def test_id2file_map(self):
        ids = ['a', 'b', 'c']

        # regular Fasta files
        for id_ in ids:
            open(join(self.tmpdir, f'{id_}.faa'), 'a').close()
        obs = id2file_map(self.tmpdir)
        exp = {x: f'{x}.faa' for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, f'{id_}.faa'))

        # gzipped Fasta files
        for id_ in ids:
            open(join(self.tmpdir, f'{id_}.faa.gz'), 'a').close()
        obs = id2file_map(self.tmpdir)
        exp = {x: f'{x}.faa.gz' for x in ids}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, f'{id_}.faa.gz'))

        # user-defined extension filename
        for id_ in ids:
            open(join(self.tmpdir, f'{id_}.faa'), 'a').close()
        open(join(self.tmpdir, 'readme.txt'), 'a').close()
        open(join(self.tmpdir, 'taxdump'), 'a').close()
        obs = id2file_map(self.tmpdir, ext='faa')
        exp = {x: f'{x}.faa' for x in ids}
        self.assertDictEqual(obs, exp)
        remove(join(self.tmpdir, 'readme.txt'))
        remove(join(self.tmpdir, 'taxdump'))

        # user-defined ID list
        obs = id2file_map(self.tmpdir, ids=['a', 'b'])
        exp = {x: f'{x}.faa' for x in ['a', 'b']}
        self.assertDictEqual(obs, exp)
        for id_ in ids:
            remove(join(self.tmpdir, f'{id_}.faa'))

        # duplicated IDs
        open(join(self.tmpdir, 'x.faa'), 'a').close()
        open(join(self.tmpdir, 'x.tsv'), 'a').close()
        with self.assertRaises(ValueError) as ctx:
            id2file_map(self.tmpdir)
        msg = 'Ambiguous files for Id: x.'
        self.assertEqual(str(ctx.exception), msg)
        remove(join(self.tmpdir, 'x.faa'))
        remove(join(self.tmpdir, 'x.tsv'))

    def test_read_fasta(self):
        # a typical NCBI-style FASTA file
        lines = (
            '>NP_052663.1 transposase Tra5 (plasmid) [Escherichia coli O157:H7'
            ' str. Sakai]',
            'MDELRAQGYHFKVKTVTESLRRHGLRAKASWNFSPVCYRAHSQPVSENLLEQDFYASGPN',
            'QKWAGDITYLRTDEGWPYLAVVTCHYWLINVVTHDGTGAPGTLCYRMSFAKMTDRYTSI',
            '>NP_052627.1 PapX (plasmid) [Escherichia coli O157:H7'
            ' str. Sakai]',
            'MLLALLSSTDNFCLSSTELSERLDVSRTYITRACDSLEKFGFIKRMESKEDRRSKNIYLT',
            'SDGNLYLQRTTRIYGRYLKKYGATLQMMKSKHLK'
        )
        obs = read_fasta(lines)
        exp = [
            ['NP_052663.1', 'transposase Tra5 (plasmid)',
             'MDELRAQGYHFKVKTVTESLRRHGLRAKASWNFSPVCYRAHSQPVSENLLEQDFYASGPNQKWA'
             'GDITYLRTDEGWPYLAVVTCHYWLINVVTHDGTGAPGTLCYRMSFAKMTDRYTSI'],
            ['NP_052627.1', 'PapX (plasmid)',
             'MLLALLSSTDNFCLSSTELSERLDVSRTYITRACDSLEKFGFIKRMESKEDRRSKNIYLTSDGN'
             'LYLQRTTRIYGRYLKKYGATLQMMKSKHLK']]
        self.assertListEqual(obs, exp)

        # a simple example without title and with trailing "*"
        lines = ('>seq1', 'ABCDEFG*')
        obs = read_fasta(lines)
        exp = [['seq1', '', 'ABCDEFG']]
        self.assertListEqual(obs, exp)

    def test_read_input_prots(self):
        # FASTA format
        fd, fp = mkstemp()
        with open(fd, 'w') as f:
            f.write('>NP_052663.1 transposase Tra5 (plasmid) [Escherichia '
                    'coli O157:H7 str. Sakai]\n')
            f.write('MDELRAQGYHFKVKTVTESLRRHGLRAKASWNFSPVCYRAHSQPVSENLL\n')
            f.write('EQDFYASGPNQKWAGDITYLRTDEGWPYLAVVTCHYWLINVVTHDGTGAP\n')
            f.write('GTLCYRMSFAKMTDRYTSI\n')
            f.write('>NP_052627.1 PapX (plasmid) [Escherichia coli O157:H7 '
                    'str. Sakai]\n')
            f.write('MLLALLSSTDNFCLSSTELSERLDVSRTYITRACDSLEKFGFIKRMESKE\n')
            f.write('DRRSKNIYLTSDGNLYLQRTTRIYGRYLKKYGATLQMMKSKHLK\n')
        obs = read_input_prots(fp)
        exp = [{'id': 'NP_052663.1', 'product': 'transposase Tra5 (plasmid)',
                'seq': 'MDELRAQGYHFKVKTVTESLRRHGLRAKASWNFSPVCYRAHSQPVSENLLEQDF'
                       'YASGPNQKWAGDITYLRTDEGWPYLAVVTCHYWLINVVTHDGTGAPGTLCYRMS'
                       'FAKMTDRYTSI'},
               {'id': 'NP_052627.1', 'product': 'PapX (plasmid)',
                'seq': 'MLLALLSSTDNFCLSSTELSERLDVSRTYITRACDSLEKFGFIKRMESKEDRRS'
                       'KNIYLTSDGNLYLQRTTRIYGRYLKKYGATLQMMKSKHLK'}]
        self.assertListEqual(obs, exp)

        # with duplicated Ids
        with open(fp, 'a') as f:
            f.write('>NP_052663.1\n')
        obs = read_input_prots(fp)
        exp[0]['dups'] = 1
        self.assertListEqual(obs, exp)

        # plain list
        with open(fp, 'w') as f:
            f.write('NP_454622.1\nNP_230502.1\nNP_384288.1\n')
        obs = read_input_prots(fp)
        exp = [{'id': 'NP_454622.1', 'product': '', 'seq': ''},
               {'id': 'NP_230502.1', 'product': '', 'seq': ''},
               {'id': 'NP_384288.1', 'product': '', 'seq': ''}]
        self.assertListEqual(obs, exp)

        # with duplicated Ids
        with open(fp, 'w') as f:
            f.write('#ids\nseq1\nseq2\nseq3\nseq2\n')
        obs = read_input_prots(fp)
        exp = [{'id': 'seq1', 'product': '', 'seq': ''},
               {'id': 'seq2', 'product': '', 'seq': '', 'dups': 1},
               {'id': 'seq3', 'product': '', 'seq': ''}]
        self.assertListEqual(obs, exp)
        remove(fp)

    def test_is_capital(self):
        self.assertTrue(is_capital('Escherichia coli'))
        self.assertTrue(is_capital('[Clostridium] difficile'))
        self.assertFalse(is_capital('unknown bacterium'))
        self.assertFalse(is_capital('[]'))
        self.assertFalse(is_capital(''))

    def test_is_latin(self):
        self.assertTrue(is_latin('Escherichia coli'))
        self.assertTrue(is_latin('Rickettsia felis'))
        self.assertTrue(is_latin('[Clostridium] difficile'))
        self.assertFalse(is_latin('Enterobacteriaceae'))
        self.assertFalse(is_latin('Escherichia coli O157:H7'))
        self.assertFalse(is_latin('Citrobacter sp. A293'))
        self.assertFalse(is_latin('bacterium LF-3'))
        self.assertFalse(is_latin('Firmicutes bacterium CAG:129'))
        self.assertFalse(is_latin('unidentified enterobacterium'))
        self.assertFalse(is_latin('Staphylococcus aureus subsp. aureus'))
        self.assertFalse(is_latin('X123 Y456'))
        self.assertFalse(is_latin(''))

    def test_contain_words(self):
        text = 'The quick brown fox jumps over the lazy dog'
        self.assertTrue(contain_words(text, ['cat', 'fox', 'dog']))
        self.assertFalse(contain_words(text, ['rat', 'tiger', 'a']))
        self.assertTrue(contain_words(text, ['Dog', 'Woodpecker']))
        text = 'This is a complete "sentence" with symbols.'
        self.assertTrue(contain_words(text, ['this']))
        self.assertTrue(contain_words(text, ['sentence']))
        self.assertTrue(contain_words(text, ['symbols']))

    def test_get_product(self):
        input = 'hypothetical protein [Candidatus organisma]'
        obs = get_product(input)
        exp = 'hypothetical protein'
        self.assertEqual(obs, exp)

    def test_seqid2accver(self):
        self.assertEqual(seqid2accver('NP_123456.1'), 'NP_123456.1')
        self.assertEqual(seqid2accver('ref|NP_123456.1|'), 'NP_123456.1')

    def test_read_taxdump(self):
        tmpdir = mkdtemp()
        exp = {'1': {'name': 'root', 'parent': '1', 'rank': 'no rank'},
               '2': {'name': 'Bacteria', 'parent': '131567', 'rank':
                     'superkingdom'}}

        # custom short format
        with open(join(tmpdir, 'nodes.dmp'), 'w') as f:
            f.write('1	|	1	|	no rank	|\n'
                    '2	|	131567	|	superkingdom	|\n')
        with open(join(tmpdir, 'names.dmp'), 'w') as f:
            f.write('1	|	root	|\n'
                    '2	|	Bacteria	|\n'
                    '2157	|	Archaea	|\n')
        obs = read_taxdump(tmpdir)
        self.assertDictEqual(obs, exp)

        # original NCBI format
        with open(join(tmpdir, 'nodes.dmp'), 'w') as f:
            f.write(
                '1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0'
                '\t|\t0\t|\t0\t|\t\t|\n'
                '2\t|\t131567\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0'
                '\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\n')
        with open(join(tmpdir, 'names.dmp'), 'w') as f:
            f.write(
                '1\t|\tall\t|\t\t|\tsynonym\t|\n'
                '1\t|\troot\t|\t\t|\tscientific name\t|\n'
                '2\t|\tBacteria\t|\tBacteria <prokaryotes>\t|\tscientific name'
                '\t|\n'
                '2\t|\tMonera\t|\tMonera <Bacteria>\t|\tin-part\t|\n')
        obs = read_taxdump(tmpdir)
        self.assertDictEqual(obs, exp)
        rmtree(tmpdir)

    def test_read_prot2taxid(self):
        fd, fp = mkstemp()
        exp = {'NP_454622.1': '220341',
               'NP_230502.1': '243277',
               'NP_720561.1': '210007'}

        # plain format
        with open(fd, 'w') as f:
            f.write('NP_454622.1\t220341\n'
                    'NP_230502.1\t243277\n'
                    'NP_720561.1\t210007\n')
        obs = read_prot2taxid(fp)
        self.assertDictEqual(obs, exp)

        # NCBI format
        with open(fp, 'w') as f:
            f.write('accession\taccession.version\ttaxid\tgi\n'
                    'NP_454622\tNP_454622.1\t220341\t16759005\n'
                    'NP_230502\tNP_230502.1\t243277\t15640871\n'
                    'NP_720561\tNP_720561.1\t210007\t24378606\n'
                    '# this is a comment')
        obs = read_prot2taxid(fp)
        self.assertDictEqual(obs, exp)
        remove(fp)

    def test__get_taxon(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        with self.assertRaises(ValueError) as ctx:
            _get_taxon('12345', taxdump)
        msg = 'TaxID 12345 is not found in taxonomy database.'
        self.assertEqual(str(ctx.exception), msg)

    def test_get_lineage(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        self.assertListEqual(get_lineage('2157', taxdump),
                             ['2157', '131567', '1'])
        obs = get_lineage('1538547', taxdump)
        exp = ['1538547', '1655637', '1655434', '1935183', '2157', '131567',
               '1']
        self.assertListEqual(obs, exp)

    def test_is_ancestral(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        self.assertTrue(is_ancestral('1538547', {'2157'}, taxdump))
        self.assertFalse(is_ancestral('1538547', {'2'}, taxdump))

    def test_taxid_at_rank(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        self.assertEqual(
            taxid_at_rank('1538547', 'genus', taxdump), '1655637')
        self.assertEqual(
            taxid_at_rank('1538547', 'phylum', taxdump), '1655434')

    def test_taxids_at_ranks(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        ranks = ['phylum', 'class', 'genus', 'species']
        obs = taxids_at_ranks('1538547', ranks, taxdump)
        exp = {'phylum': '1655434', 'class': None, 'genus': '1655637',
               'species': '1538547'}
        self.assertDictEqual(obs, exp)

    def test_find_lca(self):
        taxdump = taxdump_from_text(taxdump_archaea)

        self.assertEqual(find_lca(['131567'], taxdump), '131567')
        self.assertEqual(find_lca(['1935183', '1783276'], taxdump), '2157')
        self.assertEqual(find_lca([
            '1935183', '1783276', '1655434'], taxdump), '2157')
        self.assertEqual(find_lca([
            '1935183', '1783276', '2157'], taxdump), '2157')
        self.assertEqual(find_lca(['1935183', '2'], taxdump), '131567')
        self.assertEqual(find_lca(['1', '2', '1'], taxdump), '1')

        taxdump['x'] = {'name': 'x', 'parent': 'x'}
        with self.assertRaises(ValueError) as ctx:
            find_lca(['2', 'x'], taxdump)
        msg = 'Cannot find LCA of taxIds in database.'
        self.assertEqual(str(ctx.exception), msg)

    def test_sort_by_hierarchy(self):
        taxdump = taxdump_from_text(taxdump_archaea)

        # sort by hierarchy from low to high:
        # Lokiarchaeum sp. GC14_75, Lokiarchaeum, Candidatus Lokiarchaeota,
        # Asgard group, Archaea
        tids = ['1935183', '1655637', '2157', '1538547', '1655434']
        obs = sort_by_hierarchy(tids, taxdump)
        exp = ['1538547', '1655637', '1655434', '1935183', '2157']
        self.assertListEqual(obs, exp)

        # add DPANN group, which cannot be sorted
        tids.append('1783276')
        with self.assertRaises(ValueError) as ctx:
            sort_by_hierarchy(tids, taxdump)
        msg = 'Cannot sort taxIds by hierarchy.'
        self.assertEqual(str(ctx.exception), msg)

        # take away Candidatus Lokiarchaeota to break the sequence
        tids.pop()
        tids.pop()
        with self.assertRaises(ValueError) as ctx:
            sort_by_hierarchy(tids, taxdump)
        msg = 'Cannot sort taxIds by hierarchy.'
        self.assertEqual(str(ctx.exception), msg)

    def test_refine_taxdump(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        tids = ['2', '2157']  # Bacteria and Archaea
        obs = set(refine_taxdump(tids, taxdump))
        exp = {'1', '131567', '2', '2157'}  # plus cellular organisms and root
        self.assertSetEqual(obs, exp)

    def test_add_children(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        add_children(taxdump)
        self.assertSetEqual(
            set(taxdump['1']['children']), {'131567'})
        self.assertSetEqual(
            set(taxdump['2157']['children']), {'1935183', '1783276'})
        self.assertSetEqual(
            set(taxdump['1655434']['children']), {'1655637'})
        self.assertListEqual(taxdump['2']['children'], [])

    def test_get_descendants(self):
        taxdump = taxdump_from_text(taxdump_archaea)
        add_children(taxdump)
        obs = get_descendants('1935183', taxdump)  # Asgard group
        exp = ['1655434', '1655637', '1538547']
        self.assertListEqual(obs, exp)

    def test_taxdump_from_text(self):
        obs = taxdump_from_text(taxdump_archaea)
        self.assertEqual(len(obs), 9)
        self.assertDictEqual(
            obs['1'], {'name': 'root', 'parent': '1', 'rank': 'no rank'})
        self.assertDictEqual(
            obs['2157'], {'name': 'Archaea', 'parent': '131567',
                          'rank': 'superkingdom'})

    def test_rank_plural(self):
        self.assertEqual(rank_plural('subspecies'), 'subspecies')
        self.assertEqual(rank_plural('superphylum'), 'superphyla')
        self.assertEqual(rank_plural('infraorder'), 'infraorders')
        self.assertEqual(rank_plural('class'), 'classes')
        self.assertEqual(rank_plural('family'), 'families')
        self.assertEqual(rank_plural('genus'), 'genera')


"""Constants"""

taxdump_archaea = (
    '1,root,1,no rank',
    '131567,cellular organisms,1,no rank',
    '2,Bacteria,131567,superkingdom',
    '2157,Archaea,131567,superkingdom',
    '1935183,Asgard group,2157,no rank',
    '1783276,DPANN group,2157,no rank',
    '1655434,Candidatus Lokiarchaeota,1935183,phylum',
    '1655637,Lokiarchaeum,1655434,genus',
    '1538547,Lokiarchaeum sp. GC14_75,1655637,species')


if __name__ == '__main__':
    main()

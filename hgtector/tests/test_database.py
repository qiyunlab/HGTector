#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove, makedirs, listdir
from os.path import join, dirname, realpath, isfile, isdir
from shutil import rmtree, copy
from tempfile import mkdtemp
import gzip

import pandas as pd

from hgtector.database import Database
from hgtector.util import taxdump_from_text


class DatabaseTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

        # whether to test remote functions, which highly depend on the network
        # connection and the current status of the NCBI server
        self.test_remote = False

    def tearDown(self):
        rmtree(self.tmpdir)

    def test___call__(self):
        # TODO
        pass

    def test_set_parameters(self):
        # TODO
        pass

    def test_connect_server(self):
        # TODO
        pass

    def test_retrieve_taxdump(self):
        # TODO
        pass

    def test_retrieve_summary(self):
        # TODO
        pass

    def test_retrieve_categories(self):
        # TODO
        pass

    def test_filter_genomes(self):
        me = Database()
        header = ('# assembly_accession', 'assembly_level')
        data = (('GCF_000000001.1', 'Chromosome'),
                ('GCF_000000002.1', 'Complete Genome'),
                ('GCF_000000003.2', 'Scaffold'),
                ('GCF_000000004.1', 'Contig'),
                ('GCA_000000004.1', 'Contig'))
        df = pd.DataFrame(data, columns=header)
        me.complete = False
        me.genoids = None
        me.exclude = False

        # drop duplicates
        me.df = df.copy()
        me.filter_genomes()
        self.assertEqual(me.df.shape[0], 4)
        self.assertListEqual(me.df['genome'].tolist(), [
            'G000000001', 'G000000002', 'G000000003', 'G000000004'])
        self.assertEqual(me.df.query(
            'accession == "GCF_000000004.1"').shape[0], 1)

        # complete genomes only
        me.complete = True
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accnov'].tolist(), [
            'GCF_000000001', 'GCF_000000002'])

        # include certain genomes
        me.complete = False
        me.genoids = 'G000000001,G000000003'
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accession'].tolist(), [
            'GCF_000000001.1', 'GCF_000000003.2'])

        # exclude certain genomes
        me.genoids = ['GCF_000000002.1', 'GCF_000000004']
        me.exclude = True
        me.df = df.copy()
        me.filter_genomes()
        self.assertListEqual(me.df['accession'].tolist(), [
            'GCF_000000001.1', 'GCF_000000003.2'])

    def test_identify_taxonomy(self):
        me = Database()
        header = ('organism_name', 'taxid', 'species', 'species_taxid')
        data = (('Escherichia coli UMN026', '585056', 'E. coli', '562'),
                ('Escherichia coli O104:H4', '1038927', 'E. coli', '562'),
                ('Klebsiella aerogenes', '548', 'Klebsiella aerogenes', '548'),
                ('unclassified Gammaproteobacteria', '118884', '', ''),
                ('Plasmid pPY113', '126792', '', ''))
        df = pd.DataFrame(data, columns=header)

        # organism names must be capital and latinate
        me.capital = True
        me.block = None
        me.latin = True
        me.taxids = None
        me.exclude = False
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.df = df.copy()
        me.identify_taxonomy()
        self.assertNotIn('species_taxid', me.df.columns)
        self.assertListEqual(me.df.index.tolist(), [0, 1, 2])
        self.assertListEqual(me.df['species'].tolist(), ['562', '562', '548'])

        # block word
        me.block = 'plasmid'
        me.latin = False
        me.df = df.copy()
        me.identify_taxonomy()
        self.assertListEqual(me.df.index.tolist(), [0, 1, 2])

        # no Escherichia
        me.taxids = '561'
        me.exclude = True
        me.df = df.copy()
        me.identify_taxonomy()
        self.assertListEqual(me.df.index.tolist(), [2])

    def test_sample_by_taxonomy(self):
        me = Database()

        # do nothing
        me.sample = None
        self.assertIsNone(me.sample_by_taxonomy())

        # xxx
        header = ('genome', 'taxid', 'refseq_category', 'assembly_level')
        data = (('G1', '585056', '', 'Chromosome'),  # E. coli UMN026
                ('G2', '1038927', 'representative genome', 'Chromosome'),
                # E. coli O104:H4 (rep. genome to be prioritized over G1)
                ('G3', '2580236', '', 'Contig'),  # sync E. coli
                ('G4', '622', '', 'Scaffold'),  # Shigella
                ('G5', '548', '', 'Scaffold'),  # Klebsiella
                ('G6', '126792', 'reference genome', 'Contig'))  # plasmid
        df = pd.DataFrame(data, columns=header)
        me.reference = False
        me.representative = False
        me.taxdump = taxdump_from_text(taxdump_proteo)

        # up to one genome per genus
        me.rank = 'genus'
        me.sample = 1
        me.df = df.copy()
        me.sample_by_taxonomy()
        self.assertListEqual(me.df.columns.tolist(), list(header) + ['genus'])
        self.assertListEqual(me.df['genome'].tolist(), ['G2', 'G4', 'G5'])

        # include reference genome (plasmid)
        me.reference = True
        me.df = df.copy()
        me.sample_by_taxonomy()
        self.assertEqual(me.df['genome'].tolist()[-1], 'G6')

        # up to two genomes for entire cellular life
        me.rank = 'superkingdom'
        me.sample = 2
        me.reference = False
        me.df = df.copy()
        me.sample_by_taxonomy()
        self.assertListEqual(me.df['genome'].tolist(), ['G1', 'G2'])

    def test_download_genomes(self):
        # TODO
        pass

    def test_extract_genomes(self):
        # TODO
        pass

    def test_genome_lineages(self):
        me = Database()
        me.output = self.tmpdir
        me.taxdump = taxdump_from_text(taxdump_proteo)
        data = (('G1', '1224', ''),     # Proteobacteria
                ('G2', '562',  '562'),  # Escherichia coli
                ('G3', '622',  '622'),  # Shigella dysenteriae
                ('G4', '548',  '548'))  # Klebsiella aerogenes
        me.df = pd.DataFrame(data, columns=[
            'genome', 'taxid', 'species']).set_index('genome')
        for rank in ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                     'family', 'genus']:
            me.df[rank] = ''
        me.genome_lineages()
        with open(join(self.tmpdir, 'lineages.txt'), 'r') as f:
            obs = dict(x.split('\t') for x in f.read().splitlines())
        proteo = 'k__Bacteria; p__Proteobacteria;'
        self.assertEqual(obs['G1'], proteo + ' c__; o__; f__; g__; s__')
        entero = proteo + ' c__Gammaproteobacteria; o__Enterobacterales;' +\
            ' f__Enterobacteriaceae;'
        self.assertEqual(
            obs['G2'], entero + ' g__Escherichia; s__Escherichia coli')
        self.assertEqual(
            obs['G3'], entero + ' g__Shigella; s__Shigella dysenteriae')
        self.assertEqual(
            obs['G4'], entero + ' g__Klebsiella; s__Klebsiella aerogenes')
        remove(join(self.tmpdir, 'lineages.txt'))

    def test_genome_metadata(self):
        me = Database()
        me.output = self.tmpdir
        me.df = pd.Series({
            'genome': 'G1',
            'accession': 'GCF_000123456.1',
            'asm_name': 'ASM123v1',
            'bioproject': 'PRJNA123456',
            'biosample': 'SAMN00123456',
            'assembly_level': 'Chromosome',
            'organism_name': 'hypothetical organism',
            'infraspecific_name': '',
            'isolate': '',
            'taxid': '12345',
            'ftp_path': ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/123/'
                         '456/GCF_000123456.1_ASM123v1'),
            'proteins': 100,
            'residues': 12500,
            'whatever': 'nonsense'}).to_frame().T
        me.genome_metadata()
        with open(join(self.tmpdir, 'genomes.tsv'), 'r') as f:
            obs = f.read().splitlines()
        exp = ('genome', 'proteins', 'residues', 'assembly_level', 'accession',
               'bioproject', 'biosample', 'asm_name', 'organism_name',
               'infraspecific_name', 'isolate', 'taxid', 'ftp_path')
        self.assertEqual(obs[0], '\t'.join(exp))
        exp = ('G1', '100', '12500', 'Chromosome', 'GCF_000123456.1',
               'PRJNA123456', 'SAMN00123456', 'ASM123v1',
               'hypothetical organism', '', '', '12345',
               ('ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/123/456/'
                'GCF_000123456.1_ASM123v1'))
        self.assertEqual(obs[1], '\t'.join(exp))
        remove(join(self.tmpdir, 'genomes.tsv'))

    def test_build_taxdump(self):
        me = Database()
        me.output = self.tmpdir
        me.tmpdir = join(self.datadir, 'DnaK', 'taxdump')
        me.taxdump = taxdump_from_text(taxdump_proteo)
        data = (('G1', '1224'),     # Proteobacteria
                ('G2', '562'),      # Escherichia coli
                ('G3', '585056'),   # E. coli UMN026
                ('G4', '1038927'))  # E. coli O104:H4
        me.df = pd.DataFrame(data, columns=[
            'genome', 'taxid']).set_index('genome')
        me.build_taxdump()
        with open(join(self.tmpdir, 'taxdump', 'nodes.dmp'), 'r') as f:
            obs = set(x.split('\t')[0] for x in f.read().splitlines())
        exp = {'1', '131567', '2', '1224', '1236', '91347', '543', '561',
               '562', '585056', '1038927'}
        self.assertSetEqual(obs, exp)
        rmtree(join(self.tmpdir, 'taxdump'))

    def test_build_taxonmap(self):
        me = Database()
        me.output = self.tmpdir
        me.taxdump = taxdump_from_text(taxdump_proteo)
        me.p2tids = {'P1': {'766'},  # Rickettsiales
                     'P2': {'570', '548'},  # Klebsiella
                     'P3': {'620', '622'},  # Shigella
                     'P4': {'561', '562'},  # Escherichia
                     'P5': {'126792', '28211'}}  # root
        me.build_taxonmap()
        exp = {'P1': '766', 'P2': '570', 'P3': '620', 'P4': '561', 'P5': '1'}
        self.assertDictEqual(me.taxonmap, exp)
        with gzip.open(join(self.tmpdir, 'taxon.map.gz'), 'rt') as f:
            obs = dict(x.split('\t') for x in f.read().splitlines())
        self.assertDictEqual(obs, exp)
        remove(join(self.tmpdir, 'taxon.map.gz'))

    def test_compile_database(self):
        me = Database()
        me.output = self.tmpdir

        # don't compile
        me.compile = 'none'
        me.compile_database()
        self.assertListEqual(listdir(self.tmpdir), [])

        # get database files
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        makedirs(join(self.tmpdir, 'taxdump'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
             join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
             join(self.tmpdir, 'taxdump', 'names.dmp'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())

        # set parameters
        me.threads = 1
        me.tmpdir = self.tmpdir
        me.makeblastdb = 'makeblastdb'
        me.diamond = 'diamond'

        # compile blast database
        me.compile = 'blast'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        rmtree(join(self.tmpdir, 'blast'))

        # compile diamond database
        me.compile = 'diamond'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'diamond'))

        # compile both databases
        me.compile = 'both'
        me.compile_database()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'blast'))
        rmtree(join(self.tmpdir, 'diamond'))

        # clean up
        remove(join(self.tmpdir, 'db.faa'))
        rmtree(join(self.tmpdir, 'taxdump'))

    def test_build_blast_db(self):
        me = Database()
        me.output = self.tmpdir
        me.makeblastdb = 'makeblastdb'
        me.tmpdir = self.tmpdir
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())
        me.build_blast_db()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', f'db.{ext}')))
        rmtree(join(self.tmpdir, 'blast'))
        remove(join(self.tmpdir, 'db.faa'))

    def test_build_diamond_db(self):
        me = Database()
        me.output = self.tmpdir
        me.diamond = 'diamond'
        me.threads = 1
        me.tmpdir = self.tmpdir
        copy(join(self.datadir, 'DnaK', 'linear.faa'),
             join(self.tmpdir, 'db.faa'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())
        makedirs(join(self.tmpdir, 'taxdump'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
             join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        copy(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
             join(self.tmpdir, 'taxdump', 'names.dmp'))
        me.build_diamond_db()
        self.assertTrue(isdir(join(self.tmpdir, 'diamond')))
        self.assertTrue(isfile(join(self.tmpdir, 'diamond', 'db.dmnd')))
        rmtree(join(self.tmpdir, 'diamond'))
        remove(join(self.tmpdir, 'db.faa'))
        remove(join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        remove(join(self.tmpdir, 'taxdump', 'names.dmp'))

    def test_check_local_file(self):
        me = Database()

        # file does not exist
        file = join(self.tmpdir, 'tmp.in')
        self.assertFalse(me.check_local_file(file))

        # empty file will be deleted
        open(file, 'w').close()
        self.assertFalse(me.check_local_file(file))
        self.assertFalse(isfile(file))

        # file exists and has content
        with open(file, 'w') as f:
            f.write('Hello world!')
        self.assertTrue(isfile(file))
        self.assertTrue(me.check_local_file(file))

        # overwrite existing file
        self.assertFalse(me.check_local_file(file, overwrite=True))
        self.assertFalse(isfile(file))


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

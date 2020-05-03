#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkstemp, mkdtemp

from hgtector.search import Search


class SearchTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datadir = join(dirname(realpath(__file__)), 'data')

        # whether to test remote functions, which highly depends on network
        # connection and current status of NCBI server
        self.test_remote = False

    def tearDown(self):
        rmtree(self.tmpdir)

    def test___call__(self):
        # TODO
        pass

    """workflows"""

    def test_args_wf(self):
        # TODO
        pass

    def test_input_wf(self):
        # TODO
        pass

    def test_search_wf(self):
        me = Search()
        me.evalue = 1e-10
        me.identity = 30
        me.coverage = 30
        me.maxseqs = 5
        me.threads = 1
        me.tmpdir = self.tmpdir
        me.extrargs = None
        seqs = [(
            'WP_000516135.1',  # Proteobacteria DnaK protein
            'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEGDRTTPSIIAYTQDGETLVGQPAKRQAVT'
            'NPQNTLFAIKRLIGRRFQDEEVQRDVSIMPFKIIAADNGDAWVEVKGQKMAPPQISAEVL'
            'KKMKKTAEDYLGEPVTEAVITVPAYFNDAQRQATKDAGRIAGLEVKRIINEPTAAALAYG'
            'LDKGTGNRTIAVYDLGGGTFDISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDSRLINYL'
            'VEEFKKDQGIDLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
            'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDIDDVILVGGQTRMPMVQKKVAEFFGKE'
            'PRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLGIETMGGVMTTLIAKNTTIPT'
            'KHSQVFSTAEDNQSAVTIHVLQGERKRAADNKSLGQFNLDGINPAPRGMPQIEVTFDIDA'
            'DGILHVSAKDKNSGKEQKITIKASSGLNEDEIQKMVRDAEANAEADRKFEELVQTRNQGD'
            'HLLHSTRKQVEEAGDKLPADDKTAIESALTALETALKGEDKAAIEAKMQELAQVSQKLME'
            'IAQQQHAQQQTAGADASANNAKDDDVVDAEFEEVKDKK')]

        # pre-computed search result
        me.method = 'precomp'
        obs = me.search_wf(seqs, join(self.datadir, 'DnaK', 'blast.1.m8'))
        # note: standard BLAST output (m8) does not contain taxId(s)
        exp = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100.00', 'taxid': ''},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100.00', 'taxid': ''},
            {'id': 'NP_454622.1', 'identity': '96.865', 'evalue': '0.0',
             'score': '1253', 'coverage': '100.00', 'taxid': ''},
            {'id': 'YP_003611339.1', 'identity': '96.238', 'evalue': '0.0',
             'score': '1245', 'coverage': '100.00', 'taxid': ''},
            {'id': 'YP_005225024.1', 'identity': '95.611', 'evalue': '0.0',
             'score': '1219', 'coverage': '100.00', 'taxid': ''}]}
        self.assertDictEqual(obs, exp)

        # blast search
        me.method = 'blast'
        me.blastp = 'blastp'
        me.db = join(self.datadir, 'DnaK', 'blast', 'db')
        obs = me.search_wf(seqs)
        exp = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '1133852'},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '585056'},
            {'id': 'NP_454622.1', 'identity': '96.865', 'evalue': '0.0',
             'score': '1253', 'coverage': '100', 'taxid': '220341'},
            {'id': 'YP_003611339.1', 'identity': '96.238', 'evalue': '0.0',
             'score': '1245', 'coverage': '100', 'taxid': '716541'},
            {'id': 'YP_005225024.1', 'identity': '95.611', 'evalue': '0.0',
             'score': '1219', 'coverage': '100', 'taxid': '1125630'}]}
        self.assertDictEqual(obs, exp)

        # diamond search
        me.method = 'diamond'
        me.diamond = 'diamond'
        me.db = join(self.datadir, 'DnaK', 'diamond', 'db')
        obs = me.search_wf(seqs)
        exp = {'WP_000516135.1': [
            {'id': 'YP_002410793.1', 'identity': '100.0', 'evalue': '0.0e+00',
             'score': '1092.8', 'coverage': '100.0', 'taxid': '585056'},
            {'id': 'YP_006780959.1', 'identity': '100.0', 'evalue': '0.0e+00',
             'score': '1092.8', 'coverage': '100.0', 'taxid': '1133852'},
            {'id': 'NP_454622.1', 'identity': '97.2', 'evalue': '2.2e-313',
             'score': '1060.8', 'coverage': '100.0', 'taxid': '220341'},
            {'id': 'YP_005225024.1', 'identity': '96.7', 'evalue': '5.5e-312',
             'score': '1056.2', 'coverage': '100.0', 'taxid': '1125630'},
            {'id': 'YP_003611339.1', 'identity': '96.9', 'evalue': '2.1e-311',
             'score': '1054.3', 'coverage': '100.0', 'taxid': '716541'}]}
        self.assertDictEqual(obs, exp)

        # remote search
        if self.test_remote is False:
            return
        me.method = 'remote'
        me.server = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
        me.db = 'nr'
        me.retries = 5
        me.delay = 60
        me.timeout = 1800
        me.algorithm = 'kmerBlastp'
        me.entrez = None
        obs = me.search_wf(seqs)
        self.assertEqual(len(obs['WP_000516135.1']), 5)

    def test_taxid_wf(self):
        me = Search()
        me.method = 'precomp'
        me.taxmap = None

        # this case has two hits that miss taxIds: one is available from taxon
        # map and one needs to be looked-up from database; and another hit with
        # invalid sequence Id
        prots = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '1133852'},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': ''},
            {'id': 'NP_454622.1', 'identity': '96.865', 'evalue': '0.0',
             'score': '1253', 'coverage': '100', 'taxid': '220341'},
            {'id': 'YP_003611339.1', 'identity': '96.238', 'evalue': '0.0',
             'score': '1245', 'coverage': '100', 'taxid': ''},
            {'id': 'im_not_id', 'taxid': ''},
            {'id': 'YP_005225024.1', 'identity': '95.611', 'evalue': '0.0',
             'score': '1219', 'coverage': '100', 'taxid': '1125630'}]}
        me.prot2tid = {'YP_006780959.1': '1133852',
                       'YP_002410793.1': '585056',
                       'NP_454622.1': '220341'}

        # look up taxIds from local BLAST database
        me.blastdbcmd = 'blastdbcmd'
        me.db = join(self.datadir, 'DnaK', 'blast', 'db')
        me.fetch_enable = 'no'
        me.taxid_wf(prots)

        # confirm hit deletion
        self.assertEqual('YP_005225024.1', prots['WP_000516135.1'][4]['id'])

        # check updated hits
        self.assertListEqual(['585056', '716541'], [prots['WP_000516135.1'][
            i]['taxid'] for i in (1, 3)])

        # check updated taxId map
        exp = {'YP_006780959.1': '1133852',
               'YP_002410793.1': '585056',
               'NP_454622.1': '220341',
               'YP_003611339.1': '716541',
               'YP_005225024.1': '1125630'}
        self.assertDictEqual(me.prot2tid, exp)

        # look up a taxId from remote server
        if self.test_remote is False:
            return
        me.db = ''
        me.fetch_enable = 'yes'
        me.fetch_queries = 100
        me.fetch_retries = 3
        me.fetch_delay = 5
        me.fetch_timeout = 60
        me.fetch_server = (
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')
        prots['WP_000516135.1'][2]['taxid'] = ''
        del me.prot2tid['NP_454622.1']
        me.taxid_wf(prots)

        # check retrieved taxId
        self.assertEqual(prots['WP_000516135.1'][2]['taxid'], '220341')

    def test_taxinfo_wf(self):
        # three valid hits, one of which has a taxId that is not in taxdump
        # yet; fourth hit is fake
        me = Search()
        prots = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '1133852'},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '585056'},
            {'id': 'NP_454622.1', 'identity': '96.865', 'evalue': '0.0',
             'score': '1253', 'coverage': '100', 'taxid': '220341'},
            {'id': 'im_not_seqid', 'taxid': 'im_not_taxid'}]}
        taxdump_text = (
            '1,root,1,no rank',
            '131567,cellular organisms,1,no rank',
            '2,Bacteria,131567,superkingdom',
            '1224,Proteobacteria,2,phylum',
            '1236,Gammaproteobacteria,1224,class',
            '91347,Enterobacterales,1236,order',
            '543,Enterobacteriaceae,91347,family',
            '561,Escherichia,543,genus',
            '562,Escherichia coli,561,species',
            '585056,Escherichia coli UMN026,562,no rank',
            '1038927,Escherichia coli O104:H4,562,no rank',
            '1133852,Escherichia coli O104:H4 str. 2011C-3493,1038927,no rank')
        me.taxdump = {}
        me.badtaxids = set()
        for line in taxdump_text:
            x = line.split(',')
            me.taxdump[x[0]] = {'name': x[1], 'parent': x[2], 'rank': x[3]}
        me.fetch_enable = 'yes'
        me.output = self.tmpdir

        # update taxonomic information from remote server
        if self.test_remote is False:
            return
        me.fetch_queries = 100
        me.fetch_retries = 3
        me.fetch_delay = 5
        me.fetch_timeout = 60
        me.fetch_server = (
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')
        me.taxinfo_wf(prots)

        # confirm hit deletion
        self.assertEqual(len(prots['WP_000516135.1']), 3)
        self.assertEqual(len(prots['WP_000516135.1'][2]['id']), 'NP_454622.1')
        exp = {'name': 'Salmonella enterica',
               'parent': '590', 'rank': 'species'}
        self.assertEqual(me.taxdump['28901'], exp)

    def test_taxfilt_wf(self):
        me = Search()
        me.taxdump = {}
        taxdump_text = (
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
        for line in taxdump_text:
            x = line.split(',')
            me.taxdump[x[0]] = {'name': x[1], 'parent': x[2], 'rank': x[3]}

        me.tax_include = set(['1224'])  # Gammaproteobacteria only
        me.tax_exclude = set(['620'])  # exclude Shigella
        me.tax_capital = True  # taxon name must be capital
        me.tax_block = set(['unclassified'])  # drop unclassified organisms
        me.tax_latin = False  # species name must be Latin
        me.tax_unique = True  # taxId must be unique
        me.tax_unirank = 'species'  # species must be unique
        me.badtaxids = set(['126792'])  # drop a plasmid

        prots_text = (
            'q1	s1	100	0.0	240	100	585056\n',  # E. coli UMN026
            'q1	s2	99	1e-50	180	95	1038927\n',  # E. coli O104:H4
            'q1	s3	95	2e-40	120	93	2580236\n',  # synthetic
            'q1	s4	90	3e-30	100	90	126792\n',  # plasmid
            'q1	s5	80	4e-10	80	78	622\n',  # S. dysenteriae
            'q2	s6	99	1.5e-90	99	99	548\n',  # K. aerogenes
            'q2	s7	92	5.0e-40	90	94	548\n',  # 2nd K. aerogenes
            'q2	s8	83	2.7e-20	90	88	118884\n')  # unclassified
        prots = me.parse_def_table(prots_text)
        me.taxfilt_wf(prots)
        exp = {'q1': [{'id': 's1', 'identity': '100', 'evalue': '0.0',
                       'score': '240', 'coverage': '100', 'taxid': '585056'}],
               'q2': [{'id': 's6', 'identity': '99', 'evalue': '1.5e-90',
                       'score': '99', 'coverage': '99', 'taxid': '548'}]}
        self.assertDictEqual(prots, exp)

    def test_selfaln_wf(self):
        me = Search()
        seqs = [('NP_269215.1',
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR')]

        # look up from hit table
        me.aln_method = 'lookup'
        hits = {'NP_269215.1': [
            {'id': 'NP_269215.1', 'evalue': '9.36e-77', 'score': '207'}]}
        obs = me.selfaln_wf(seqs, hits)
        exp = {'NP_269215.1': '207'}
        self.assertDictEqual(obs, exp)

        # built-in algorithm
        me.aln_method = 'fast'
        obs = me.selfaln_wf(seqs)
        exp = {'NP_269215.1': '203.8'}
        self.assertDictEqual(obs, exp)

        # blast
        me.aln_method = 'native'
        me.tmpdir = self.tmpdir
        me.threads = 1
        me.extrargs = None
        me.method = 'blast'
        me.blastp = 'blastp'
        obs = me.selfaln_wf(seqs)
        exp = {'NP_269215.1': '207'}
        self.assertDictEqual(obs, exp)

        # diamond
        me.method = 'diamond'
        me.diamond = 'diamond'
        obs = me.selfaln_wf(seqs)
        exp = {'NP_269215.1': '203.8'}
        self.assertDictEqual(obs, exp)

        # diamond can't align this sequence, do built-in instead
        seqs_ = [('YP_208127.1',
                  'MKKLLIAAMMAAALAACSQEAKQEVKEAAQAVESDVKDTAASAAESAASA'
                  'VEEAKGQVKDAAADAKASAEEAVTEAKDAAAETKEAVSEAAKDTLNKAAD'
                  'AAQEAADKMKDAAK')]
        obs = me.selfaln_wf(seqs_)
        exp = {'YP_208127.1': '206.8'}
        self.assertDictEqual(obs, exp)

        # remote
        if self.test_remote is False:
            return
        me.method = 'remote'
        me.aln_server = 'https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi'
        me.retries = 5
        me.delay = 60
        me.timeout = 1800
        me.maxchars = 2000
        me.db = 'nr'
        obs = me.selfaln_wf(seqs)
        exp = {'NP_269215.1': '207'}
        self.assertDictEqual(obs, exp)

    """input/output functions"""

    def test_subset_seqs(self):

        def helper(res):
            return [','.join([x[0] for x in y]) for y in res]

        seqs = [('seq1', 'MGKIIGIDLGTTNSCVAIMDGTQARVLENAEGDRTTPSIIAYTQDGETLV'),
                ('seq2', 'GQPAKRQAVTNPQNTLFAIKRLIGRRFQDEEVQRDVSIMPYKIIGADNGD'),
                ('seq3', 'AWLDVKGQKMAPPQISAEVLKKMKKTAEDYLGEPVTEAVITVPAYFNDAQ'),
                ('seq4', 'RQATKDAGRIAGLEVKRIINEPTAAALAYGLDKEVGNRTI'),
                ('seq5', 'AVYDLGGGTFDISIIEIDEVDGEKTFEVLATNGDTHLGGE'),
                ('seq6', 'DFDTRLINYLVDEFKKDQGIDLRNDPLAMQ'),
                ('seq7', 'RLKEAAEKAKIELSSAQQTDVNLPYITADA'),
                ('seq8', 'TGPKHMNIKVTRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDINDVILV'),
                ('seq9', 'GGQTRMPMVQKKVAEFFGKEPRKDVNPDEAVAIGAAVQGGVLTGDVKDVL'),
                ('seq10', 'LLDVTPLSLGIETMGGVMTPLITKNTTIPTKHSQV')]

        # no subsetting
        obs = helper(Search.subset_seqs(seqs))
        exp = ['seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8,seq9,seq10']
        self.assertListEqual(obs, exp)

        # by queries
        obs = helper(Search.subset_seqs(seqs, 3, None))
        exp = ['seq1,seq2,seq3', 'seq4,seq5,seq6', 'seq7,seq8,seq9', 'seq10']
        self.assertListEqual(obs, exp)

        # by maxchars
        obs = helper(Search.subset_seqs(seqs, None, 100))
        exp = [
            'seq1,seq2', 'seq3,seq4', 'seq5,seq6,seq7', 'seq8,seq9', 'seq10']
        self.assertListEqual(obs, exp)

        # by queries and maxchars
        obs = helper(Search.subset_seqs(seqs, 3, 135))
        exp = ['seq1,seq2', 'seq3,seq4,seq5', 'seq6,seq7,seq8', 'seq9,seq10']
        self.assertListEqual(obs, exp)

        # sequence length exceeds maxchars
        msg = 'Sequence seq1 exceeds maximum allowed length 40 for search.'
        with self.assertRaisesRegex(ValueError, msg):
            Search.subset_seqs(seqs, 3, 40)

    def test_update_search_results(self):
        me = Search()
        me.maxhits = None
        prots = [{'id': 'prot1'},
                 {'id': 'prot2', 'hits': []},
                 {'id': 'prot3'},
                 {'id': 'prot4'},
                 {'id': 'prot5'}]
        res = {'prot1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                          'score': '240', 'coverage': '100', 'taxid': '1234'},
                         {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                          'score': '180', 'coverage': '90', 'taxid': '567'},
                         {'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                          'score': '185', 'coverage': '90', 'taxid': ''}],
               'prot2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                          'score': '99', 'coverage': '100', 'taxid': '534'}],
               'prot3': [{'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                          'score': '185', 'coverage': '90', 'taxid': ''},
                         {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                          'score': '90', 'coverage': '80', 'taxid': '882'}]}

        # default behavior
        me.update_search_results(prots, res)
        exp = [{'id': 'prot1', 'hits': [
                   {'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                    'score': '240', 'coverage': '100', 'taxid': '1234'},
                   {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                    'score': '180', 'coverage': '90', 'taxid': '567'},
                   {'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                    'score': '185', 'coverage': '90', 'taxid': ''}]},
               {'id': 'prot2', 'hits': []},
               {'id': 'prot3', 'hits': [
                   {'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                    'score': '185', 'coverage': '90', 'taxid': ''},
                   {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                    'score': '90', 'coverage': '80', 'taxid': '882'}]},
               {'id': 'prot4'},
               {'id': 'prot5'}]
        self.assertListEqual(prots, exp)

        # with indices, one of which has no hit, plus maxhits
        for i in (0, 2):
            del prots[i]['hits']
        indices = set([0, 2, 3])
        me.maxhits = 2
        me.update_search_results(prots, res, indices)
        exp = [{'id': 'prot1', 'hits': [
                   {'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                    'score': '240', 'coverage': '100', 'taxid': '1234'},
                   {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                    'score': '180', 'coverage': '90', 'taxid': '567'}]},
               {'id': 'prot2', 'hits': []},
               {'id': 'prot3', 'hits': [
                   {'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                    'score': '185', 'coverage': '90', 'taxid': ''},
                   {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                    'score': '90', 'coverage': '80', 'taxid': '882'}]},
               {'id': 'prot4', 'hits': []},
               {'id': 'prot5'}]
        self.assertListEqual(prots, exp)

    def test_write_search_results(self):
        prots = [{'id': 'prot1', 'product': 'my protein 1', 'score': 125,
                  'seq': 'ABCDEFGHIJKLMN', 'hits': [
                      {'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                       'score': '240', 'coverage': '100', 'taxid': '12345'},
                      {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                       'score': '180', 'coverage': '100', 'taxid': '678'}]},
                 {'id': 'prot3', 'product': 'my protein 2', 'score': 96,
                  'seq': 'OPQRSTUVWXYZ', 'hits': [
                      {'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                       'score': '99', 'coverage': '93', 'taxid': '910'},
                      {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                       'score': '90', 'coverage': '63', 'taxid': '1112'}]}]
        file = join(self.tmpdir, 'tmp.out')
        with open(file, 'w') as f:
            Search.write_search_results(f, prots)
        with open(file, 'r') as f:
            obs = f.read()
        remove(file)
        exp = ('# ID: prot1\n'
               '# Length: 14\n'
               '# Product: my protein 1\n'
               '# Score: 125\n'
               '# Hits: 2\n'
               'sub1	100	0.0	240	100	12345\n'
               'sub2	95	1.2e-30	180	100	678\n'
               '# ID: prot3\n'
               '# Length: 12\n'
               '# Product: my protein 2\n'
               '# Score: 96\n'
               '# Hits: 2\n'
               'sub4	99	1.5e-90	99	93	910\n'
               'sub5	85	5.0e-20	90	63	1112\n')
        self.assertEqual(obs, exp)

    def test_parse_prev_results(self):
        fd, fp = mkstemp()
        with open(fd, 'w') as f:
            f.write(
                '# ID: prot1\n'
                '# Length: 100\n'
                '# Product: im a protein\n'
                '# Score: 500\n'
                '# Hits: 1\n'
                'seq1  100.000 0.0     533     100.00  12345\n'
                '# ID: prot2\n'
                '# Length: 85\n'
                '# Product: im another protein\n'
                '# Score: 275\n'
                '# Hits: 3\n'
                'seq2  100.000 3.01e-36        125     100.00  12345\n'
                'seq3  71.212  2.44e-18        80.1    100.00  678\n'
                'seq4  66.154  4.88e-17        76.6    100.00  910\n')
        prots = [{'id': 'prot1'}, {'id': 'prot2'}, {'id': 'prot3'},
                 {'id': 'prot4'}, {'id': 'prot5'}]
        obs = Search.parse_prev_results(fp, prots)
        exp = ['prot1', 'prot2']
        self.assertListEqual(obs, exp)
        exp = [{'id': 'prot1', 'score': 0, 'hits': []},
               {'id': 'prot2', 'score': 0, 'hits': []},
               {'id': 'prot3'}, {'id': 'prot4'}, {'id': 'prot5'}]
        self.assertListEqual(prots, exp)
        remove(fp)

    def test_check_missing_seqs(self):
        data = {
            'set1': {'prots': [
                {'id': 'id1', 'product': 'prod1', 'seq': 'ABCDEFG'},
                {'id': 'id2', 'product': 'prod2', 'seq': 'HIJKLMN'},
                {'id': 'id3', 'product': '', 'seq': ''}]},
            'set2': {'prots': [
                {'id': 'id4', 'product': 'prod4', 'seq': ''},
                {'id': 'id5', 'product': '', 'seq': '', 'hits': []}]},
            'set3': {'prots': [
                {'id': 'id6', 'product': '', 'seq': ''}], 'done': 1}}
        obs = Search.check_missing_seqs(data)
        exp = ['id3', 'id4']
        self.assertListEqual(obs, exp)

    def test_update_dmp_files(self):
        me = Search()
        me.taxdump = {}
        me.output = self.tmpdir
        with open(join(self.datadir, 'efetch_taxonomy.xml'), 'r') as f:
            xml = f.read()
        me.parse_taxonomy_xml(xml)
        tids = ['2', '561', '1301', '1496', '1239']
        me.update_dmp_files(tids)
        file = join(self.tmpdir, 'nodes.dmp')
        with open(file, 'r') as f:
            obs = f.read()
        exp = ('2	|	131567	|	superkingdom	|\n'
               '561	|	543	|	genus	|\n'
               '1239	|	1783272	|	phylum	|\n'
               '1301	|	1300	|	genus	|\n'
               '1496	|	1870884	|	species	|\n')
        self.assertEqual(obs, exp)
        remove(file)
        fp = join(self.tmpdir, 'names.dmp')
        with open(fp, 'r') as f:
            obs = f.read()
        exp = ('2	|	Bacteria	|\n'
               '561	|	Escherichia	|\n'
               '1239	|	Firmicutes	|\n'
               '1301	|	Streptococcus	|\n'
               '1496	|	Clostridioides difficile	|\n')
        self.assertEqual(obs, exp)
        remove(fp)

    """sequence query functions"""

    def test_blast_seqinfo(self):
        me = Search()
        me.blastdbcmd = 'blastdbcmd'
        me.db = join(self.datadir, 'DnaK', 'blast', 'db')

        # a normal case
        ids = ['NP_454622.1', 'NP_230502.1']
        obs = me.blast_seqinfo(ids)
        exp = [('NP_454622.1', '220341', 'DnaK protein',
                'MGKIIGIDLGTTNSCVAIMDGTQARVLENAEGDRTTPSIIAYTQDGETLV'
                'GQPAKRQAVTNPQNTLFAIKRLIGRRFQDEEVQRDVSIMPYKIIGADNGD'
                'AWLDVKGQKMAPPQISAEVLKKMKKTAEDYLGEPVTEAVITVPAYFNDAQ'
                'RQATKDAGRIAGLEVKRIINEPTAAALAYGLDKEVGNRTIAVYDLGGGTF'
                'DISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDTRLINYLVDEFKKDQGI'
                'DLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
                'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDINDVILVGGQTRMPMVQ'
                'KKVAEFFGKEPRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLG'
                'IETMGGVMTPLITKNTTIPTKHSQVFSTAEDNQSAVTIHVLQGERKRASD'
                'NKSLGQFNLDGINPAPRGMPQIEVTFDIDADGILHVSAKDKNSGKEQKIT'
                'IKASSGLNEEEIQKMVRDAEANAESDRKFEELVQTRNQGDHLLHSTRKQV'
                'EEAGDKLPADDKTAIESALSALETALKGEDKAAIEAKMQELAQVSQKLME'
                'IAQQQHAQQQAGSADASANNAKDDDVVDAEFEEVKDKK'),
               ('NP_230502.1', '243277', 'molecular chaperone DnaK',
                'MGKIIGIDLGTTNSCVAVLDGDKPRVIENAEGERTTPSVIAYTDGETLVG'
                'QPAKRQAVTNPQNTLFAIKRLIGRRFEDEEVQRDIKIMPYKIVKADNGDA'
                'WVEAKGQKMAAPQVSAEVLKKMKKTAEDFLGEPVTAAVITVPAYFNDAQR'
                'QATKDAGRIAGLEVKRIINEPTAAALAYGLDKQGGDRTIAVYDLGGGTFD'
                'ISIIEIDEVEGEKTFEVLSTNGDTHLGGEDFDNRMINYLVDEFKKDQGID'
                'LRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKVT'
                'RAKLEALVEDLVQRSLEPLKVALADADLSVNDITDVILVGGQTRMPMVQK'
                'KVAEFFGKEPRKDVNPDEAVAVGAAVQGGVLAGEVKDVLLLDVTPLSLGI'
                'ETMGGVMTKLIEKNTTIPTKANQVFSTAEDNQSAVTIHVLQGERKQAMYN'
                'KSLGQFNLEGINPAPRGMPQIEVIFDLDADGILHVSAKDKQTGKEQKITI'
                'QASGGLSDAEIEKMVQEAEANKEADKKFEELATARNQADQMIHATRKQIT'
                'EAGEALPADEKAKIETAINELETAKKGEDKAEIDAKVQALMAAAQKLMEI'
                'AQQQAQAQGANAGQSSAKEDDVVDAEFEEVNDDKK')]
        self.assertListEqual(obs, exp)

        # invalid sequence Id
        obs = me.blast_seqinfo(['im_not_seqid', 'nor_am_i'])
        self.assertEqual(len(obs), 0)

        # invalid database
        me.db = 'im_not_db'
        msg = 'Invalid BLAST database'
        with self.assertRaisesRegex(ValueError, msg):
            me.blast_seqinfo(ids)

    def test_remote_seqinfo(self):
        if self.test_remote is False:
            return
        me = Search()
        me.fetch_queries = 100
        me.fetch_retries = 3
        me.fetch_delay = 5
        me.fetch_timeout = 60
        me.fetch_server = (
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')
        ids = ['NP_454622.1', 'NP_230502.1']
        obs = me.remote_seqinfo(ids)
        exp = [['NP_454622.1', '220341', 'DnaK protein',
                'MGKIIGIDLGTTNSCVAIMDGTQARVLENAEGDRTTPSIIAYTQDGETLV'
                'GQPAKRQAVTNPQNTLFAIKRLIGRRFQDEEVQRDVSIMPYKIIGADNGD'
                'AWLDVKGQKMAPPQISAEVLKKMKKTAEDYLGEPVTEAVITVPAYFNDAQ'
                'RQATKDAGRIAGLEVKRIINEPTAAALAYGLDKEVGNRTIAVYDLGGGTF'
                'DISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDTRLINYLVDEFKKDQGI'
                'DLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
                'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDINDVILVGGQTRMPMVQ'
                'KKVAEFFGKEPRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLG'
                'IETMGGVMTPLITKNTTIPTKHSQVFSTAEDNQSAVTIHVLQGERKRASD'
                'NKSLGQFNLDGINPAPRGMPQIEVTFDIDADGILHVSAKDKNSGKEQKIT'
                'IKASSGLNEEEIQKMVRDAEANAESDRKFEELVQTRNQGDHLLHSTRKQV'
                'EEAGDKLPADDKTAIESALSALETALKGEDKAAIEAKMQELAQVSQKLME'
                'IAQQQHAQQQAGSADASANNAKDDDVVDAEFEEVKDKK'],
               ['NP_230502.1', '243277', 'molecular chaperone DnaK',
                'MGKIIGIDLGTTNSCVAVLDGDKPRVIENAEGERTTPSVIAYTDGETLVG'
                'QPAKRQAVTNPQNTLFAIKRLIGRRFEDEEVQRDIKIMPYKIVKADNGDA'
                'WVEAKGQKMAAPQVSAEVLKKMKKTAEDFLGEPVTAAVITVPAYFNDAQR'
                'QATKDAGRIAGLEVKRIINEPTAAALAYGLDKQGGDRTIAVYDLGGGTFD'
                'ISIIEIDEVEGEKTFEVLSTNGDTHLGGEDFDNRMINYLVDEFKKDQGID'
                'LRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKVT'
                'RAKLEALVEDLVQRSLEPLKVALADADLSVNDITDVILVGGQTRMPMVQK'
                'KVAEFFGKEPRKDVNPDEAVAVGAAVQGGVLAGEVKDVLLLDVTPLSLGI'
                'ETMGGVMTKLIEKNTTIPTKANQVFSTAEDNQSAVTIHVLQGERKQAMYN'
                'KSLGQFNLEGINPAPRGMPQIEVIFDLDADGILHVSAKDKQTGKEQKITI'
                'QASGGLSDAEIEKMVQEAEANKEADKKFEELATARNQADQMIHATRKQIT'
                'EAGEALPADEKAKIETAINELETAKKGEDKAEIDAKVQALMAAAQKLMEI'
                'AQQQAQAQGANAGQSSAKEDDVVDAEFEEVNDDKK']]
        self.assertListEqual(obs, exp)

    def test_update_prot_seqs(self):
        me = Search()
        me.data = {
            'set1': {'prots': [
                {'id': 'id1', 'product': 'prod1', 'seq': 'ABCDEFG'},
                {'id': 'id2', 'product': 'prod2', 'seq': 'HIJKLMN'},
                {'id': 'id3', 'product': '', 'seq': ''}]},
            'set2': {'prots': [
                {'id': 'id4', 'product': 'prod4', 'seq': ''},
                {'id': 'id5', 'product': '', 'seq': ''}]},
            'set3': {'prots': [
                {'id': 'id6', 'product': '', 'seq': ''}]}}
        seqs = [('id1', '0', 'prod1', 'ABCDEFG'),
                ('id3', '0', 'prod3', 'OPQRSTU'),
                ('id6.1', '0', 'prod6', 'VWXYZAB')]
        obs = me.update_prot_seqs(seqs)
        self.assertEqual(obs, 2)
        exp = {
            'set1': {'prots': [
                {'id': 'id1', 'product': 'prod1', 'seq': 'ABCDEFG'},
                {'id': 'id2', 'product': 'prod2', 'seq': 'HIJKLMN'},
                {'id': 'id3', 'product': 'prod3', 'seq': 'OPQRSTU'}]},
            'set2': {'prots': [
                {'id': 'id4', 'product': 'prod4', 'seq': ''},
                {'id': 'id5', 'product': '', 'seq': ''}]},
            'set3': {'prots': [
                {'id': 'id6', 'product': 'prod6', 'seq': 'VWXYZAB'}]}}
        self.assertDictEqual(me.data, exp)

    def test_parse_fasta_xml(self):
        with open(join(self.datadir, 'efetch_fasta.xml'), 'r') as f:
            xml = f.read()
        obs = Search.parse_fasta_xml(xml)
        exp = [['NP_454622.1', '220341', 'DnaK protein',
                'MGKIIGIDLGTTNSCVAIMDGTQARVLENAEGDRTTPSIIAYTQDGETLV'
                'GQPAKRQAVTNPQNTLFAIKRLIGRRFQDEEVQRDVSIMPYKIIGADNGD'
                'AWLDVKGQKMAPPQISAEVLKKMKKTAEDYLGEPVTEAVITVPAYFNDAQ'
                'RQATKDAGRIAGLEVKRIINEPTAAALAYGLDKEVGNRTIAVYDLGGGTF'
                'DISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDTRLINYLVDEFKKDQGI'
                'DLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
                'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDINDVILVGGQTRMPMVQ'
                'KKVAEFFGKEPRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLG'
                'IETMGGVMTPLITKNTTIPTKHSQVFSTAEDNQSAVTIHVLQGERKRASD'
                'NKSLGQFNLDGINPAPRGMPQIEVTFDIDADGILHVSAKDKNSGKEQKIT'
                'IKASSGLNEEEIQKMVRDAEANAESDRKFEELVQTRNQGDHLLHSTRKQV'
                'EEAGDKLPADDKTAIESALSALETALKGEDKAAIEAKMQELAQVSQKLME'
                'IAQQQHAQQQAGSADASANNAKDDDVVDAEFEEVKDKK'],
               ['NP_230502.1', '243277', 'molecular chaperone DnaK',
                'MGKIIGIDLGTTNSCVAVLDGDKPRVIENAEGERTTPSVIAYTDGETLVG'
                'QPAKRQAVTNPQNTLFAIKRLIGRRFEDEEVQRDIKIMPYKIVKADNGDA'
                'WVEAKGQKMAAPQVSAEVLKKMKKTAEDFLGEPVTAAVITVPAYFNDAQR'
                'QATKDAGRIAGLEVKRIINEPTAAALAYGLDKQGGDRTIAVYDLGGGTFD'
                'ISIIEIDEVEGEKTFEVLSTNGDTHLGGEDFDNRMINYLVDEFKKDQGID'
                'LRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKVT'
                'RAKLEALVEDLVQRSLEPLKVALADADLSVNDITDVILVGGQTRMPMVQK'
                'KVAEFFGKEPRKDVNPDEAVAVGAAVQGGVLAGEVKDVLLLDVTPLSLGI'
                'ETMGGVMTKLIEKNTTIPTKANQVFSTAEDNQSAVTIHVLQGERKQAMYN'
                'KSLGQFNLEGINPAPRGMPQIEVIFDLDADGILHVSAKDKQTGKEQKITI'
                'QASGGLSDAEIEKMVQEAEANKEADKKFEELATARNQADQMIHATRKQIT'
                'EAGEALPADEKAKIETAINELETAKKGEDKAEIDAKVQALMAAAQKLMEI'
                'AQQQAQAQGANAGQSSAKEDDVVDAEFEEVNDDKK']]
        self.assertListEqual(obs, exp)

    """homology search functions"""

    def test_blast_search(self):
        me = Search()
        me.db = join(self.datadir, 'DnaK', 'blast', 'db')

        # single query sequence
        seqs = [(
            'WP_000516135.1',  # Proteobacteria DnaK protein
            'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEGDRTTPSIIAYTQDGETLVGQPAKRQAVT'
            'NPQNTLFAIKRLIGRRFQDEEVQRDVSIMPFKIIAADNGDAWVEVKGQKMAPPQISAEVL'
            'KKMKKTAEDYLGEPVTEAVITVPAYFNDAQRQATKDAGRIAGLEVKRIINEPTAAALAYG'
            'LDKGTGNRTIAVYDLGGGTFDISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDSRLINYL'
            'VEEFKKDQGIDLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
            'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDIDDVILVGGQTRMPMVQKKVAEFFGKE'
            'PRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLGIETMGGVMTTLIAKNTTIPT'
            'KHSQVFSTAEDNQSAVTIHVLQGERKRAADNKSLGQFNLDGINPAPRGMPQIEVTFDIDA'
            'DGILHVSAKDKNSGKEQKITIKASSGLNEDEIQKMVRDAEANAEADRKFEELVQTRNQGD'
            'HLLHSTRKQVEEAGDKLPADDKTAIESALTALETALKGEDKAAIEAKMQELAQVSQKLME'
            'IAQQQHAQQQTAGADASANNAKDDDVVDAEFEEVKDKK')]
        me.evalue = 1e-10
        me.identity = 30
        me.coverage = 30
        me.maxseqs = 5
        me.blastp = 'blastp'
        me.threads = 1
        me.tmpdir = self.tmpdir
        obs = me.blast_search(seqs)
        exp = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '1133852'},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '1290', 'coverage': '100', 'taxid': '585056'},
            {'id': 'NP_454622.1', 'identity': '96.865', 'evalue': '0.0',
             'score': '1253', 'coverage': '100', 'taxid': '220341'},
            {'id': 'YP_003611339.1', 'identity': '96.238', 'evalue': '0.0',
             'score': '1245', 'coverage': '100', 'taxid': '716541'},
            {'id': 'YP_005225024.1', 'identity': '95.611', 'evalue': '0.0',
             'score': '1219', 'coverage': '100', 'taxid': '1125630'}]}
        self.assertDictEqual(obs, exp)

        # with extra arguments
        me.extrargs = '-word_size 2 -gapopen 9 -gapextend 1'
        obs = me.blast_search(seqs)
        # note: actual numbers depend on software environment
        exp = {'WP_000516135.1': [
            {'id': 'YP_006780959.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '974', 'coverage': '100', 'taxid': '1133852'},
            {'id': 'YP_002410793.1', 'identity': '100.000', 'evalue': '0.0',
             'score': '974', 'coverage': '100', 'taxid': '585056'},
            {'id': 'NP_454622.1',    'identity': '96.865', 'evalue': '0.0',
             'score': '954', 'coverage': '100', 'taxid': '220341'},
            {'id': 'YP_005225024.1', 'identity': '95.611', 'evalue': '0.0',
             'score': '943', 'coverage': '100', 'taxid': '1125630'},
            {'id': 'YP_003611339.1', 'identity': '95.925', 'evalue': '0.0',
             'score': '935', 'coverage': '100', 'taxid': '716541'}]}
        self.assertDictEqual(obs, exp)

        # invalid blastp executable
        me.blastp = 'im_not_blastp'
        msg = 'blastp failed with error code'
        with self.assertRaisesRegex(ValueError, msg):
            me.blast_search(seqs)

    def test_diamond_search(self):
        me = Search()
        me.db = join(self.datadir, 'DnaK', 'diamond', 'db')

        # single query sequence
        seqs = [(
            'WP_000516135.1',  # Proteobacteria DnaK protein
            'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEGDRTTPSIIAYTQDGETLVGQPAKRQAVT'
            'NPQNTLFAIKRLIGRRFQDEEVQRDVSIMPFKIIAADNGDAWVEVKGQKMAPPQISAEVL'
            'KKMKKTAEDYLGEPVTEAVITVPAYFNDAQRQATKDAGRIAGLEVKRIINEPTAAALAYG'
            'LDKGTGNRTIAVYDLGGGTFDISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDSRLINYL'
            'VEEFKKDQGIDLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
            'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDIDDVILVGGQTRMPMVQKKVAEFFGKE'
            'PRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLGIETMGGVMTTLIAKNTTIPT'
            'KHSQVFSTAEDNQSAVTIHVLQGERKRAADNKSLGQFNLDGINPAPRGMPQIEVTFDIDA'
            'DGILHVSAKDKNSGKEQKITIKASSGLNEDEIQKMVRDAEANAEADRKFEELVQTRNQGD'
            'HLLHSTRKQVEEAGDKLPADDKTAIESALTALETALKGEDKAAIEAKMQELAQVSQKLME'
            'IAQQQHAQQQTAGADASANNAKDDDVVDAEFEEVKDKK')]
        me.evalue = 1e-10
        me.identity = 30
        me.coverage = 30
        me.maxseqs = 5
        me.diamond = 'diamond'
        me.threads = 1
        me.tmpdir = self.tmpdir
        obs = me.diamond_search(seqs)
        # note: actual numbers depend on software environment
        exp = {'WP_000516135.1': [
            {'id': 'YP_002410793.1', 'identity': '100.0', 'evalue': '0.0e+00',
             'score': '1092.8', 'coverage': '100.0', 'taxid': '585056'},
            {'id': 'YP_006780959.1', 'identity': '100.0', 'evalue': '0.0e+00',
             'score': '1092.8', 'coverage': '100.0', 'taxid': '1133852'},
            {'id': 'NP_454622.1', 'identity': '97.2', 'evalue': '2.2e-313',
             'score': '1060.8', 'coverage': '100.0', 'taxid': '220341'},
            {'id': 'YP_005225024.1', 'identity': '96.7', 'evalue': '5.5e-312',
             'score': '1056.2', 'coverage': '100.0', 'taxid': '1125630'},
            {'id': 'YP_003611339.1', 'identity': '96.9', 'evalue': '2.1e-311',
             'score': '1054.3', 'coverage': '100.0', 'taxid': '716541'}]}
        self.assertDictEqual(obs, exp)

        # with extra arguments
        me.extrargs = '--block-size 3.0 --gapopen 9 --gapextend 1'
        obs = me.diamond_search(seqs)
        # note: actual numbers depend on software environment
        exp = {'WP_000516135.1': [
            {'id': 'YP_002410793.1', 'identity': '100.0', 'evalue': '1.0e-241',
             'score': '822.7', 'coverage': '100.0', 'taxid': '585056'},
            {'id': 'YP_006780959.1', 'identity': '100.0', 'evalue': '1.0e-241',
             'score': '822.7', 'coverage': '100.0', 'taxid': '1133852'},
            {'id': 'NP_454622.1', 'identity': '97.0', 'evalue': '7.7e-234',
             'score': '796.6', 'coverage': '100.0', 'taxid': '220341'},
            {'id': 'YP_005225024.1', 'identity': '96.6', 'evalue': '2.7e-233',
             'score': '794.8', 'coverage': '100.0', 'taxid': '1125630'},
            {'id': 'YP_003611339.1', 'identity': '96.9', 'evalue': '1.4e-232',
             'score': '792.4', 'coverage': '100.0', 'taxid': '716541'}]}
        self.assertDictEqual(obs, exp)

        # invalid diamond executable
        me.diamond = 'im_not_diamond'
        msg = 'diamond failed with error code'
        with self.assertRaisesRegex(ValueError, msg):
            me.diamond_search(seqs)

    def test_remote_search(self):
        if self.test_remote is False:
            return
        me = Search()
        me.server = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
        me.db = 'nr'
        me.retries = 5
        me.delay = 60
        me.timeout = 1800
        me.algorithm = 'kmerBlastp'
        me.evalue = 1e-10
        me.identity = 30
        me.coverage = 30
        me.maxseqs = 5
        me.extrargs = None
        me.entrez = None

        # single query sequence
        seqs = [('WP_000516135.1',
                 'MGKIIGIDLGTTNSCVAIMDGTTPRVLENAEGDRTTPSIIAYTQDGETLVGQPAKRQAVT'
                 'NPQNTLFAIKRLIGRRFQDEEVQRDVSIMPFKIIAADNGDAWVEVKGQKMAPPQISAEVL'
                 'KKMKKTAEDYLGEPVTEAVITVPAYFNDAQRQATKDAGRIAGLEVKRIINEPTAAALAYG'
                 'LDKGTGNRTIAVYDLGGGTFDISIIEIDEVDGEKTFEVLATNGDTHLGGEDFDSRLINYL'
                 'VEEFKKDQGIDLRNDPLAMQRLKEAAEKAKIELSSAQQTDVNLPYITADATGPKHMNIKV'
                 'TRAKLESLVEDLVNRSIEPLKVALQDAGLSVSDIDDVILVGGQTRMPMVQKKVAEFFGKE'
                 'PRKDVNPDEAVAIGAAVQGGVLTGDVKDVLLLDVTPLSLGIETMGGVMTTLIAKNTTIPT'
                 'KHSQVFSTAEDNQSAVTIHVLQGERKRAADNKSLGQFNLDGINPAPRGMPQIEVTFDIDA'
                 'DGILHVSAKDKNSGKEQKITIKASSGLNEDEIQKMVRDAEANAEADRKFEELVQTRNQGD'
                 'HLLHSTRKQVEEAGDKLPADDKTAIESALTALETALKGEDKAAIEAKMQELAQVSQKLME'
                 'IAQQQHAQQQTAGADASANNAKDDDVVDAEFEEVKDKK')]
        obs = me.remote_search(seqs)
        self.assertEqual(len(obs['WP_000516135.1']), 5)

        # multiple query sequences and more parameters
        seqs = [('NP_052663.1',
                 'MDELRAQGYHFKVKTVTESLRRHGLRAKASWNFSPVCYRAHSQPVSENLLEQDFYA'
                 'SGPNQKWAGDITYLRTDEGWPYLAVVTCHYWLINVVTHDGTGAPGTLCYRMSFAKM'
                 'TDRYTSI'),
                ('NP_052627.1',
                 'MLLALLSSTDNFCLSSTELSERLDVSRTYITRACDSLEKFGFIKRMESKEDRRSKN'
                 'IYLTSDGNLYLQRTTRIYGRYLKKYGATLQMMKSKHLK')]
        me.maxseqs = 3
        me.entrez = 'all [filter] NOT(metagenomes[orgn])'
        me.extrargs = 'BLAST_SPEC=MicrobialGenomes'
        obs = me.remote_search(seqs)
        self.assertEqual(len(obs['NP_052663.1']), 3)
        self.assertEqual(len(obs['NP_052627.1']), 3)

    def test_parse_def_table(self):
        me = Search()

        # a simple case
        me.maxhits = 2
        lines = ('qry1	sub1	100	0.0	240	100	1234\n',
                 'qry1	sub2	95	1.2e-30	180	90	567\n',
                 'qry1	sub3	80	4.5e-10	96	75	89\n',
                 'qry2	sub4	99	1.5e-90	99	100	534\n',
                 'qry2	sub5	85	5.0e-20	90	80	882\n')
        obs = me.parse_def_table(lines)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100', 'taxid': '1234'},
                        {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                         'score': '180', 'coverage': '90', 'taxid': '567'}],
               'qry2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                         'score': '99', 'coverage': '100', 'taxid': '534'},
                        {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                         'score': '90', 'coverage': '80', 'taxid': '882'}]}
        self.assertDictEqual(obs, exp)

        # with more filters
        me.maxhits = 0
        me.evalue = 1e-10
        me.identity = 60
        me.coverage = 90
        lines = ('qry1	sub1	100	0.0	240	100	N/A\n',
                 'qry1	sub2	99	1.1e-85	200	100	N/A\n',
                 'qry1	sub3	97	2.2e-60	185	90	N/A\n',
                 'qry1	sub4	95	3.3e-40	150	80	N/A\n',
                 'qry1	sub5	90	4.4e-10	125	70	N/A\n',
                 'qry1	sub5	80	5.5e-5	100	60	N/A\n',
                 'qry1	sub5	50	1.23	75	50	N/A\n')
        obs = me.parse_def_table(lines)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100', 'taxid': ''},
                        {'id': 'sub2', 'identity': '99', 'evalue': '1.1e-85',
                         'score': '200', 'coverage': '100', 'taxid': ''},
                        {'id': 'sub3', 'identity': '97', 'evalue': '2.2e-60',
                         'score': '185', 'coverage': '90', 'taxid': ''}]}
        self.assertDictEqual(obs, exp)

    def test_parse_m8_table(self):
        me = Search()

        # default (no filter)
        lenmap = {'qry1': 110, 'qry2': 105}
        lines = ('qry1	ref|sub1|	100	108	0	0	1	110	1	108	0.0	240\n',
                 'qry1	sub2	95	105	3	0	1	110	5	112	1.2e-30	180\n',
                 'qry1	sub3	80	96	8	1	5	110	5	107	4.5e-10	96\n',
                 'qry2	sub4	99	90	11	1	5	102	3	130	1.5e-90	99\n',
                 'qry2	sub5	85	52	13	3	45	110	148	192	5.0e-20	90\n')
        obs = me.parse_m8_table(lines, lenmap)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100.00', 'taxid': ''},
                        {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                         'score': '180', 'coverage': '100.00', 'taxid': ''},
                        {'id': 'sub3', 'identity': '80', 'evalue': '4.5e-10',
                         'score': '96', 'coverage': '96.36', 'taxid': ''}],
               'qry2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                         'score': '99', 'coverage': '93.33', 'taxid': ''},
                        {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                         'score': '90', 'coverage': '62.86', 'taxid': ''}]}
        self.assertDictEqual(obs, exp)

        # with filters
        me.maxhits = 2
        me.coverage = 90
        obs = me.parse_m8_table(lines, lenmap)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100.00', 'taxid': ''},
                        {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                         'score': '180', 'coverage': '100.00', 'taxid': ''}],
               'qry2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                         'score': '99', 'coverage': '93.33', 'taxid': ''}]}
        self.assertDictEqual(obs, exp)

    def test_parse_hit_table(self):
        me = Search()
        fd, fp = mkstemp()

        # default format
        with open(fd, 'w') as f:
            f.write('qry1	sub1	100	0.0	240	100	1234\n'
                    'qry1	sub2	95	1.2e-30	180	90	567\n'
                    'qry1	sub3	80	4.5e-10	96	75	89\n'
                    'qry2	sub4	99	1.5e-90	99	100	534\n'
                    'qry2	sub5	85	5.0e-20	90	80	882\n')
        obs = me.parse_hit_table(fp)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100', 'taxid': '1234'},
                        {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                         'score': '180', 'coverage': '90', 'taxid': '567'},
                        {'id': 'sub3', 'identity': '80', 'evalue': '4.5e-10',
                         'score': '96', 'coverage': '75', 'taxid': '89'}],
               'qry2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                         'score': '99', 'coverage': '100', 'taxid': '534'},
                        {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                         'score': '90', 'coverage': '80', 'taxid': '882'}]}
        self.assertDictEqual(obs, exp)

        # m8 format
        with open(fp, 'w') as f:
            f.write('# I am a comment line\n'
                    'qry1	sub1	100	108	0	0	1	110	1	108	0.0	'
                    '240	some\n'
                    'qry1	sub2	95	105	3	0	1	110	5	112	1.2e-30	'
                    '180	thing\n'
                    'qry1	sub3	80	96	8	1	5	110	5	107	4.5e-10	'
                    '96	extra\n'
                    'qry2	sub4	99	90	11	1	5	102	3	130	1.5e-90	'
                    '99	is\n'
                    'qry2	sub5	85	52	13	3	45	110	148	192	5.0e-20	'
                    '90	here\n')
        lenmap = {'qry1': 110, 'qry2': 105}
        me.maxhits = 2
        obs = me.parse_hit_table(fp, lenmap)
        exp = {'qry1': [{'id': 'sub1', 'identity': '100', 'evalue': '0.0',
                         'score': '240', 'coverage': '100.00', 'taxid': ''},
                        {'id': 'sub2', 'identity': '95', 'evalue': '1.2e-30',
                         'score': '180', 'coverage': '100.00', 'taxid': ''}],
               'qry2': [{'id': 'sub4', 'identity': '99', 'evalue': '1.5e-90',
                         'score': '99', 'coverage': '93.33', 'taxid': ''},
                        {'id': 'sub5', 'identity': '85', 'evalue': '5.0e-20',
                         'score': '90', 'coverage': '62.86', 'taxid': ''}]}
        self.assertDictEqual(obs, exp)
        remove(fp)

    """taxonomy query functions"""

    def test_update_hit_taxids(self):
        me = Search()
        me.prot2tid = {'seq1': '561', 'seq3': '590', 'seq6': '258'}
        prots = {'prot1': [{'id': 'seq1', 'taxid': '561'},
                           {'id': 'seq2', 'taxid': '562'},
                           {'id': 'seq3', 'taxid': ''},
                           {'id': 'seq4', 'taxid': ''}],
                 'prot2': [{'id': 'seq5', 'taxid': '2'},
                           {'id': 'seq6', 'taxid': ''},
                           {'id': 'seq7', 'taxid': ''}]}

        # basic
        obs = me.update_hit_taxids(prots)
        exp = (['seq4', 'seq7'], ['2', '562'])
        self.assertTupleEqual(obs, exp)
        exp = {'seq1': '561', 'seq2': '562', 'seq3': '590', 'seq5': '2',
               'seq6': '258'}
        self.assertDictEqual(me.prot2tid, exp)
        exp = {'prot1': [{'id': 'seq1', 'taxid': '561'},
                         {'id': 'seq2', 'taxid': '562'},
                         {'id': 'seq3', 'taxid': '590'},
                         {'id': 'seq4', 'taxid': ''}],
               'prot2': [{'id': 'seq5', 'taxid': '2'},
                         {'id': 'seq6', 'taxid': '258'},
                         {'id': 'seq7', 'taxid': ''}]}
        self.assertDictEqual(prots, exp)

        # with reference taxon map
        taxmap = {'seq4': '816', 'seq7': '838'}
        obs = me.update_hit_taxids(prots, taxmap)
        exp = ([], ['816', '838'])
        self.assertTupleEqual(obs, exp)

    def test_remote_taxinfo(self):
        if self.test_remote is False:
            return
        me = Search()
        me.taxdump = {}
        me.fetch_queries = 100
        me.fetch_retries = 3
        me.fetch_delay = 5
        me.fetch_timeout = 60
        me.fetch_server = (
            'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi')

        # Escherichia coli, Streptococcus pneumoniae and Clostridioides
        # difficile
        tids = ['561', '1313', '1496']
        xml = me.remote_taxinfo(tids)
        me.parse_taxonomy_xml(xml)
        obs = []
        for tid, d in sorted(me.taxdump.items(), key=lambda x: int(x[0])):
            obs.append(f'{tid}:{d["name"]},{d["parent"]},{d["rank"]}')
        exp = ['2:Bacteria,131567,superkingdom',
               '543:Enterobacteriaceae,91347,family',
               '561:Escherichia,543,genus',
               '1224:Proteobacteria,2,phylum',
               '1236:Gammaproteobacteria,1224,class',
               '1239:Firmicutes,1783272,phylum',
               '1300:Streptococcaceae,186826,family',
               '1301:Streptococcus,1300,genus',
               '1313:Streptococcus pneumoniae,1301,species',
               '1496:Clostridioides difficile,1870884,species',
               '91061:Bacilli,1239,class',
               '91347:Enterobacterales,1236,order',
               '131567:cellular organisms,1,no rank',
               '186801:Clostridia,1239,class',
               '186802:Clostridiales,186801,order',
               '186804:Peptostreptococcaceae,186802,family',
               '186826:Lactobacillales,91061,order',
               '1783272:Terrabacteria group,2,no rank',
               '1870884:Clostridioides,186804,genus']
        self.assertListEqual(obs, exp)

    def test_parse_taxonomy_xml(self):
        me = Search()
        me.taxdump = {}
        with open(join(self.datadir, 'efetch_taxonomy.xml'), 'r') as f:
            xml = f.read()
        obs = me.parse_taxonomy_xml(xml)
        exp = ['561', '543', '91347', '1236', '1224', '2', '131567',
               '1313', '1301', '1300', '186826', '91061', '1239', '1783272',
               '1496', '1870884', '186804', '186802', '186801']
        self.assertListEqual(obs, exp)
        obs = []
        for tid, d in sorted(me.taxdump.items(), key=lambda x: int(x[0])):
            obs.append(f'{tid}:{d["name"]},{d["parent"]},{d["rank"]}')
        exp = ['2:Bacteria,131567,superkingdom',
               '543:Enterobacteriaceae,91347,family',
               '561:Escherichia,543,genus',
               '1224:Proteobacteria,2,phylum',
               '1236:Gammaproteobacteria,1224,class',
               '1239:Firmicutes,1783272,phylum',
               '1300:Streptococcaceae,186826,family',
               '1301:Streptococcus,1300,genus',
               '1313:Streptococcus pneumoniae,1301,species',
               '1496:Clostridioides difficile,1870884,species',
               '91061:Bacilli,1239,class',
               '91347:Enterobacterales,1236,order',
               '131567:cellular organisms,1,no rank',
               '186801:Clostridia,1239,class',
               '186802:Clostridiales,186801,order',
               '186804:Peptostreptococcaceae,186802,family',
               '186826:Lactobacillales,91061,order',
               '1783272:Terrabacteria group,2,no rank',
               '1870884:Clostridioides,186804,genus']
        self.assertListEqual(obs, exp)

    """self-alignment functions"""

    def test_lookup_selfaln(self):
        seqs = [('seq1', 'ABCDE'), ('seq2', 'FGHIJ'), ('seq3', 'KLMNO')]
        hits = {'seq1': [
            {'id': 'seq1', 'evalue': '1e-50', 'score': '100'},
            {'id': 'seq2', 'evalue': '1e-30', 'score': '90'},
            {'id': 'seq4', 'evalue': '1e-20', 'score': '75'}],
                'seq2': [
            {'id': 'seq2', 'evalue': '1e-45', 'score': '95'},
            {'id': 'seq1', 'evalue': '1e-35', 'score': '90'},
            {'id': 'seq5', 'evalue': '1e-10', 'score': '25'}],
                'seq3': [
            {'id': 'seq3', 'evalue': '1e-60', 'score': '135'}]}
        obs = Search.lookup_selfaln(seqs, hits)
        exp = [('seq1', '100', '1e-50'),
               ('seq2', '95', '1e-45'),
               ('seq3', '135', '1e-60')]
        self.assertListEqual(obs, exp)

    def test_fast_selfaln(self):
        # query sequence: NP_269215.1 (Cas9 protein), 1-100 aa
        seq = ('MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
               'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR')
        obs = Search.fast_selfaln(seq)
        exp = ('203.8', '4.6e-58')
        self.assertTupleEqual(obs, exp)

    def test_blast_selfaln(self):
        # note: with one sequence or multiple sequences,
        # the e-values are different
        me = Search()
        me.blastp = 'blastp'
        me.threads = 1
        me.tmpdir = self.tmpdir

        # single sequence
        seqs = [('NP_269215.1',
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR')]
        obs = me.blast_selfaln(seqs)
        exp = [('NP_269215.1', '207', '9.36e-77')]
        self.assertListEqual(obs, exp)

        # with extra arguments
        me.extrargs = '-word_size 2 -gapopen 9 -gapextend 1'
        obs = me.blast_selfaln(seqs)
        exp = [('NP_269215.1', '140', '2.24e-53')]
        self.assertListEqual(obs, exp)

        # multiple sequences
        me.extrargs = None
        seqs = [('NP_269215.1',    # Cas9, 1-100aa
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR'),
                ('NP_009225.1',    # BRCA1, 1-150aa
                 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLK'
                 'LLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLE'
                 'YANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQET'),
                ('NP_000537.3',    # p53, 1-200aa
                 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'
                 'EQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ'
                 'KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDST'
                 'PPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGN')]
        obs = me.blast_selfaln(seqs)
        exp = [('NP_269215.1', '207', '4.21e-76'),
               ('NP_009225.1', '305', '2.86e-113'),
               ('NP_000537.3', '407', '4.69e-152')]
        self.assertListEqual(obs, exp)

    def test_diamond_selfaln(self):
        me = Search()
        me.diamond = 'diamond'
        me.threads = 1
        me.tmpdir = self.tmpdir

        # single sequence
        seqs = [('NP_269215.1',
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR')]
        obs = me.diamond_selfaln(seqs)
        exp = [('NP_269215.1', '203.8', '4.6e-58')]
        self.assertListEqual(obs, exp)

        # with extra arguments
        me.extrargs = '--block-size 3.0 --gapopen 9 --gapextend 1'
        obs = me.diamond_selfaln(seqs)
        exp = [('NP_269215.1', '160.3', '5.6e-45')]
        self.assertListEqual(obs, exp)

        # multiple sequences
        me.extrargs = None
        seqs = [('NP_269215.1',    # Cas9, 1-100aa
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR'),
                ('NP_009225.1',    # BRCA1, 1-150aa
                 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLK'
                 'LLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLE'
                 'YANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQET'),
                ('NP_000537.3',    # p53, 1-200aa
                 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'
                 'EQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ'
                 'KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDST'
                 'PPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGN')]
        obs = me.diamond_selfaln(seqs)
        exp = [('NP_269215.1', '203.8', '2.1e-57'),
               ('NP_009225.1', '300.4', '2.4e-86'),
               ('NP_000537.3', '370.5', '2.6e-107')]
        self.assertListEqual(obs, exp)

    def test_remote_selfaln(self):
        if self.test_remote is False:
            return
        me = Search()
        me.aln_server = 'https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi'
        me.retries = 5
        me.delay = 60
        me.timeout = 1800
        me.db = 'nr'
        me.extrargs = None
        me.maxchars = None

        # single sequence
        seqs = [('NP_269215.1',
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR')]
        obs = me.remote_selfaln(seqs)
        exp = [('NP_269215.1', '207', '9.36e-77')]
        self.assertListEqual(obs, exp)

        # multiple sequences
        seqs = [('NP_269215.1',
                 'MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGA'
                 'LLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHR'),
                ('NP_009225.1',
                 'MDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLK'
                 'LLNQKKGPSQCPLCKNDITKRSLQESTRFSQLVEELLKIICAFQLDTGLE'
                 'YANSYNFAKKENNSPEHLKDEVSIIQSMGYRNRAKRLLQSEPENPSLQET'),
                ('NP_000537.3',
                 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDI'
                 'EQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ'
                 'KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDST'
                 'PPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGN')]
        obs = me.remote_selfaln(seqs)
        exp = [('NP_269215.1', '207', '4.21e-76'),
               ('NP_009225.1', '305', '2.86e-113'),
               ('NP_000537.3', '407', '4.69e-152')]
        self.assertListEqual(obs, exp)

    def test_parse_self_m8(self):
        # normal case
        lines = ['NP_269215.1	NP_269215.1	100.000	100	0	0	1	100	1	'
                 '100	4.21e-76	207',
                 'NP_009225.1	NP_009225.1	100.000	150	0	0	1	150	1	'
                 '150	2.86e-113	305',
                 'NP_000537.3	NP_000537.3	100.000	200	0	0	1	200	1	'
                 '200	4.69e-152	407']
        obs = Search.parse_self_m8(lines)
        exp = [('NP_269215.1', '207', '4.21e-76'),
               ('NP_009225.1', '305', '2.86e-113'),
               ('NP_000537.3', '407', '4.69e-152')]
        self.assertListEqual(obs, exp)

        # unusual case
        lines = ['# header',
                 'seq1 	seq1	100	*	*',
                 'seq1	seq1	100	200	0	0	1	200	1	200	1e-100	500',
                 'seq1	seq1	100	50	0	0	1	50	1	50	1e-10	100',
                 'seq1	seq2	95	150	20	0	1	175	25	200	1e-50	300']
        obs = Search.parse_self_m8(lines)
        exp = [('seq1', '500', '1e-100')]
        self.assertListEqual(obs, exp)


if __name__ == '__main__':
    main()

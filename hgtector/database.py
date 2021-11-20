#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os import remove, makedirs, cpu_count
from os.path import join, isfile, isdir, getsize, basename
from shutil import which, rmtree
from time import sleep
import ftplib
import gzip
import tarfile
from tempfile import mkdtemp

import pandas as pd

from hgtector.util import (
    timestamp, load_configs, get_config, arg2bool, run_command,
    list_from_param, read_taxdump, contain_words, is_latin, is_capital,
    taxid_at_rank, taxids_at_ranks, is_ancestral, find_lca)


description = """build reference protein sequence and taxonomy database"""

arguments = [
    'default protocol',
    ['--default',    'apply the default protocol', {'action': 'store_true'}],

    'basic',
    ['-o|--output',  'output database directory', {'required': True}],
    ['-c|--cats',    'enter one or more of the following RefSeq genome '
                     'categories, separated by comma: archaea, bacteria, '
                     'fungi, invertebrate, plant, protozoa, vertebrate_'
                     'mammalian, vertebrate_other, viral; "all" for all '
                     'categories, or "microbe" for archaea, bacteria, fungi, '
                     'and protozoa', {'default': 'microbe'}],

    'custom list',
    ['-t|--taxids',  'provide a custom list of TaxIDs; only genomes under '
                     'these taxonomic groups will be included'],
    ['-g|--genoids', 'provide a custom list of genome IDs, which are NCBI '
                     'assembly accessions (e.g., "GCF_000123456.1"), with or '
                     'without version number, or simply "G000123456"'],
    ['-e|--exclude', 'exclude instead of include TaxIDs and genomes defined '
                     'by custom lists', {'action': 'store_true'}],

    'taxon sampling',
    ['-s|--sample',  'sample up to this number of genomes per taxonomic '
                     'group at given rank (0 for all)', {'type': int}],
    ['-r|--rank',    'taxonomic rank at which subsampling will be performed',
                     {'default': 'species'}],

    'genome sampling',
    ['--genbank',    'also search GenBank if a genome is not found in RefSeq; '
                     'otherwise only include RefSeq genomes',
                     {'action': 'store_true'}],
    ['--complete',   'only include complete genomes or chromosomes',
                     {'action': 'store_true'}],
    ['--reference',  'include NCBI-defined reference genomes',
                     {'action': 'store_true'}],
    ['--represent',  'include NCBI-defined representative genomes',
                     {'action': 'store_true'}],
    ['--typemater',  'include NCBI-defined type material genomes',
                     {'action': 'store_true'}],

    'taxonomic filter',
    ['--capital',    'organism name must be capitalized',
                     {'choices': ['yes', 'no']}],
    ['--block',      'ignore organism names that contain any of these words'],
    ['--latin',      'genomes must have Latinate species names, e.g. "Vibrio '
                     'cholerae", but not "Vibrio sp. 123"',
                     {'choices': ['yes', 'no']}],

    'download',
    ['--overwrite',  'overwrite existing files with newly downloaded ones; '
                     'othewise use existing files whenever available, i.e., '
                     'resume an interrupted run', {'action': 'store_true'}],
    ['--retries',    'number of trials for downloading each data file',
                     {'type': int}],
    ['--delay',      'seconds between two retries', {'type': int}],
    ['--timeout',    'seconds seconds before program gives up waiting',
                     {'type': int}],

    'compile',
    ['--compile',    'compile database using this program if available',
                     {'choices': ['diamond', 'blast', 'both', 'none']}],
    ['--diamond',    'diamond executable'],
    ['--makeblastdb', 'makeblastdb executable'],
    ['--threads',    'number of threads for diamond makedb', {'type': int}],
    ['--tmpdir',     'temporary directory for diamond makedb']
]


class Database(object):

    def __init__(self):
        self.arguments = arguments
        self.description = description

    def __call__(self, args):
        print(f'Database building started at {timestamp()}.')

        # load configurations
        self.cfg = load_configs()

        # read and validate arguments
        self.set_parameters(args)

        # connect to NCBI FTP server
        self.connect_server()

        # create target directory
        makedirs(self.output, exist_ok=True)

        # retrieve taxonomy database
        self.retrieve_taxdump()

        # retrieve genome assembly summary
        self.retrieve_summary()

        # retrieve genome categories
        self.retrieve_categories()

        # filter genomes
        self.filter_genomes()

        # identify taxonomy of genomes
        self.identify_taxonomy()

        # sample genomes by taxonomy
        self.sample_by_taxonomy()

        # download genomes
        self.download_genomes()

        # extract genomes
        self.extract_genomes()

        # identify genome lineages
        self.genome_lineages()

        # write genome metadata
        self.genome_metadata()

        # build taxonomy database
        self.build_taxdump()

        # build protein-to-taxonomy map
        self.build_taxonmap()

        # build protein sequence database
        self.compile_database()

        # clean up
        if hasattr(self, 'mkdtemp'):
            rmtree(self.tmpdir)

        print(f'Database building finished at {timestamp()}.')

    def set_parameters(self, args):
        """Workflow for validating and setting arguments.

        Parameters
        ----------
        args : dict
            command-line arguments
        """
        # load arguments
        for key, val in vars(args).items():
            setattr(self, key, val)

        # load configurations
        for key in 'capital', 'block', 'latin':
            get_config(self, key, f'taxonomy.{key}')
        for key in 'retries', 'delay', 'timeout':
            get_config(self, key, f'download.{key}')
        for key in 'diamond', 'makeblastdb':
            get_config(self, key, f'program.{key}')
        for key in 'threads', 'tmpdir':
            get_config(self, key, f'local.{key}')

        # convert boolean values
        for key in 'capital', 'latin':
            setattr(self, key, arg2bool(getattr(self, key, None)))

        # make temporary directory
        if not self.tmpdir:
            self.tmpdir = mkdtemp()
            setattr(self, 'mkdtemp', True)  # mark for cleanup
        if not isdir(self.tmpdir):
            raise ValueError(f'Invalid temporary directory: {self.tmpdir}.')

        # check local executables
        for key, exe in {'blast': 'makeblastdb', 'diamond': 'diamond'}.items():
            if self.compile in (key, 'both'):
                if getattr(self, exe) is None:
                    setattr(self, exe, exe)
                if which(getattr(self, exe)) is None:
                    raise ValueError(
                        f'Invalid {exe} executable: {getattr(self, exe)}.')

        # default protocol
        if self.default:
            print('The default protocol is selected for database building.')
            print('The program will download all protein sequences of NCBI '
                  'RefSeq genomes of bacteria, archaea, fungi and protozoa, '
                  'keep one genome per species, plus all NCBI-defined '
                  'reference and representative genomes.')
            self.cats = 'microbe'
            self.sample = 1
            self.rank = 'species'
            self.reference = True
            self.represent = True

            if self.diamond is None:
                self.diamond = 'diamond'
            if which(self.diamond) is None:
                print('WARNING: Cannot find DIAMOND in this computer. '
                      'You will need to manually compile the database '
                      'after download is complete.')
                self.compile = 'none'
            else:
                self.compile = 'diamond'

        # determine number of CPUs to use
        if self.compile in ('diamond', 'both') and not self.threads:
            self.threads = cpu_count()
            if self.threads is None:
                self.threads = 1

        makedirs(self.output, exist_ok=True)

    def connect_server(self):
        """Connect to the NCBI FTP server.
        """
        print('Connecting to the NCBI FTP server...', end='', flush=True)
        makedirs(join(self.output, 'download'), exist_ok=True)
        self.ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov', timeout=self.timeout)
        self.ftp.login()
        print(' done.')

    def retrieve_taxdump(self):
        """Retrieve NCBI taxdump."""
        fname = 'taxdump.tar.gz'
        remote_file = f'/pub/taxonomy/{fname}'
        local_file = join(self.output, 'download', fname)

        # download taxdump
        if not self.check_local_file(local_file, self.overwrite):
            print('Downloading NCBI taxonomy database...', end='', flush=True)
            with open(local_file, 'wb') as f:
                self.ftp.retrbinary('RETR ' + remote_file, f.write)
            print(' done.')

        # read taxdump
        print('Reading NCBI taxonomy database...', end='', flush=True)
        with tarfile.open(local_file, 'r:gz') as f:
            f.extract('names.dmp', self.tmpdir)
            f.extract('nodes.dmp', self.tmpdir)
        self.taxdump = read_taxdump(self.tmpdir)
        print(' done.')
        print(f'Total number of TaxIDs: {len(self.taxdump)}.')

    def retrieve_summary(self, genbank=False):
        """Retrieve genome assembly summary.
        """

        def get_summary(target):
            key = target.lower()
            fname = f'assembly_summary_{key}.txt'
            remote_file = f'/genomes/{key}/{fname}'
            local_file = join(self.output, 'download', fname)

            # download summary
            if not self.check_local_file(local_file, self.overwrite):
                print(f'Downloading {target} assembly summary...', end='',
                      flush=True)
                with open(local_file, 'wb') as f:
                    self.ftp.retrbinary('RETR ' + remote_file, f.write)
                print(' done.')

            # read summary
            print(f'Reading {target} assembly summary...', end='', flush=True)
            df = pd.read_table(local_file, skiprows=1, low_memory=False)
            print(' done.')
            return df

        self.df = get_summary('RefSeq')
        if self.genbank:
            self.df = pd.concat([self.df, get_summary('GenBank')])
        print(f'Total number of genomes: {self.df.shape[0]}.')

    def retrieve_categories(self):
        """Retrieve genome categories.
        """
        # parse categories
        cats = self.cats
        if cats == 'all':
            return
        elif cats == 'microbe':
            self.cats = ['archaea', 'bacteria', 'fungi', 'protozoa']
        else:
            self.cats = cats.split(',')
        print(f'Genome categories:  {", ".join(self.cats)}')

        def get_categories(target):
            key = target.lower()

            # validate categories
            self.ftp.cwd(f'/genomes/{key}')
            dirs = [x[0] for x in self.ftp.mlsd() if x[1]['type'] == 'dir']
            for cat in self.cats:
                if cat not in dirs:
                    raise ValueError(
                        f'"{cat}" is not a valid {target} genome category.')

            # get genome list per category
            print(f'Downloading genome list per {target} category...')
            dir_ = join(self.output, 'download', 'cats')
            makedirs(dir_, exist_ok=True)
            fname = 'assembly_summary.txt'
            asms = []
            file_ = join(self.tmpdir, 'tmp.txt')

            for cat in self.cats:
                file = join(dir_, f'{key}_{cat}.txt')
                islocal = self.check_local_file(file)

                # use local file
                if islocal:
                    with open(file, 'r') as f:
                        asms_ = f.read().splitlines()

                # download remote file
                else:
                    with open(file_, 'wb') as f:
                        self.ftp.retrbinary(f'RETR {cat}/{fname}', f.write)
                    with open(file_, 'r') as f:
                        asms_ = [x.split('\t', 1)[0] for x in f.read().
                                 splitlines() if not x.startswith('#')]
                    with open(file, 'w') as f:
                        f.write(''.join([x + '\n' for x in asms_]))

                print(f'  {cat}: {len(asms_)}')
                asms += asms_

            if isfile(file_):
                remove(file_)
            print('Done.')
            return asms

        asmset = set(get_categories('RefSeq'))
        if self.genbank:
            asmset.update(get_categories('GenBank'))

        # filter genomes by category
        self.df = self.df[self.df['# assembly_accession'].isin(asmset)]
        print(f'Total number of genomes in categories: {self.df.shape[0]}.')

        # close and reconnect later to avoid some problems
        self.ftp.close()

    def filter_genomes(self):
        """Filter genomes based on genome metadata.
        """
        print('Filtering genomes...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        # complete genomes only
        if self.complete:
            self.df = self.df[self.df['assembly_level'].isin(
                {'Complete Genome', 'Chromosome'})]
            report_diff('Dropped {} non-complete genomes.')

        # non-redundant genome IDs
        # typically not necessary, just in case
        self.df.rename(columns={'# assembly_accession': 'accession'},
                       inplace=True)
        self.df['accnov'] = self.df['accession'].str.split('.', 1).str[0]
        self.df['genome'] = 'G' + self.df['accnov'].str.split('_', 1).str[-1]
        self.df.drop_duplicates(subset=['genome'], inplace=True)

        # include/exclude genome Ids
        if self.genoids:
            self.genoids = set(list_from_param(self.genoids))
            print(f'{"Ex" if self.exclude else "In"}cluding '
                  f'{len(self.genoids)} custom genome IDs...')
            self.df = self.df[(self.df['accession'].isin(self.genoids) |
                               self.df['accnov'].isin(self.genoids) |
                               self.df['genome'].isin(self.genoids)) !=
                              self.exclude]
            report_diff('Dropped {} genomes.')
        print('Done.')

    def identify_taxonomy(self):
        """Identify taxonomy of genomes.
        """
        print('Identifying taxonomy of genomes...')
        n = self.df.shape[0]

        def report_diff(msg):
            nonlocal n
            n_ = self.df.shape[0]
            if n_ < n:
                print('  ' + msg.format(n - n_))
            n = n_

        # remove non-capitalized organism names
        if self.capital:
            self.df = self.df[self.df['organism_name'].apply(is_capital)]
            report_diff('Dropped {} genomes without captalized organism name.')

        # block certain words in organism names
        if self.block:
            self.block = list_from_param(self.block)
            self.df = self.df[~self.df['organism_name'].apply(
                contain_words, args=(self.block,))]
            report_diff('Dropped {} genomes with one or more blocked words in '
                        'organism name.')

        # remove original species information
        self.df.drop(columns=['species_taxid'], inplace=True)

        # drop genomes whose taxIds are not in taxdump
        self.df.dropna(subset=['taxid'], inplace=True)
        self.df['taxid'] = self.df['taxid'].astype(str)
        self.df = self.df[self.df['taxid'].isin(self.taxdump)]
        report_diff('Dropped {} genomes without valid taxId.')

        # assign genomes to species (represented by taxId not name)
        self.df['species'] = self.df['taxid'].apply(
            taxid_at_rank, rank='species', taxdump=self.taxdump)

        # drop genomes without species taxId
        self.df.dropna(subset=['species'], inplace=True)
        report_diff('Dropped {} genomes without valid species taxId.')

        # drop genomes without Latinate species name
        if self.latin:
            self.df = self.df[self.df['species'].apply(
                lambda x: is_latin(self.taxdump[x]['name']))]
            report_diff('Dropped {} genomes without Latinate species name.')
        print('Done.')

        # include/exclude taxIds
        if self.taxids:
            self.taxids = set(list_from_param(self.taxids))
            print(f'{"Ex" if self.exclude else "In"}cluding '
                  f'{len(self.taxids)} custom TaxIDs...')

            self.df = self.df[self.df['taxid'].apply(
                lambda x: is_ancestral(x, self.taxids, self.taxdump)
                != self.exclude)]
            report_diff('Dropped {} genomes.')

    def sample_by_taxonomy(self):
        """Taxonomy-based sampling.
        """
        if not self.sample:
            return
        print('Sampling genomes based on taxonomy...')
        print(f'Up to {self.sample} genome(s) will be sampled per {self.rank}'
              '.')

        # assign genomes to given rank
        if self.rank not in self.df.columns:
            self.df[self.rank] = self.df['taxid'].apply(
                taxid_at_rank, rank=self.rank, taxdump=self.taxdump)

        # sort by reference > representative > type material > other
        self.df['rc_seq'] = self.df.apply(
            lambda x: 0 if x['refseq_category'] == 'reference genome'
            else (1 if x['refseq_category'] == 'representative genome'
                  else (2 if pd.notnull(x['relation_to_type_material'])
                  else 3)), axis=1)

        # sort by complete > scaffold > contig
        self.df['al_seq'] = self.df['assembly_level'].map(
            {'Chromosome': 0, 'Complete Genome': 0, 'Scaffold': 1,
             'Contig': 2})

        # sort genomes by three criteria
        self.df.sort_values(by=['rc_seq', 'al_seq', 'genome'], inplace=True)

        # select up to given number of genomes of each taxonomic group
        selected = set(self.df.groupby(self.rank).head(self.sample)['genome'])
        if not selected:
            raise ValueError(f'No genome is classified at rank "{self.rank}".')

        # add reference genomes
        if self.reference:
            print('Add reference genomes to selection.')
            selected.update(self.df.query(
                'refseq_category == "reference genome"')['genome'].tolist())

        # add representative genomes
        if self.represent:
            print('Add representative genomes to selection.')
            selected.update(self.df.query(
                'refseq_category == "representative genome"')[
                    'genome'].tolist())

        # add type material genomes
        if self.typemater:
            print('Add type material genomes to selection.')
            selected.update(self.df[self.df[
                'relation_to_type_material'].notna()]['genome'].tolist())

        # filter genomes to selected
        self.df.query('genome in @selected', inplace=True)
        n = self.df.shape[0]
        if n == 0:
            raise ValueError('No genome is retained after sampling.')
        print(f'Total number of sampled genomes: {n}.')

        # sort by genome ID
        self.df.sort_values('genome', inplace=True)

        # clean up temporary columns
        self.df.drop(columns=['al_seq', 'rc_seq'], inplace=True)

    def download_genomes(self):
        """Download genomes from NCBI.
        """
        # reconnect to avoid server timeout problem
        self.ftp = ftplib.FTP('ftp.ncbi.nlm.nih.gov', timeout=self.timeout)
        self.ftp.login()
        self.ftp.cwd('/genomes/all')

        print('Downloading non-redundant genomic data from NCBI...',
              flush=True)
        dir_ = join(self.output, 'download', 'faa')
        makedirs(dir_, exist_ok=True)
        stems = {}
        failures = []
        for row in self.df.itertuples():
            g = row.genome
            remote_dir = row.ftp_path.split('/', 5)[-1]
            stem = remote_dir.rsplit('/', 1)[-1]
            stems[g] = stem
            fname = f'{stem}_protein.faa.gz'
            file = join(dir_, fname)
            if self.check_local_file(file):
                success = True
            else:
                success = False
                for i in range(self.retries):
                    try:
                        with open(file, 'wb') as f:
                            self.ftp.retrbinary(
                                f'RETR {remote_dir}/{fname}', f.write)
                        print('  ' + g, flush=True)
                        success = True
                    except ftplib.error_perm as resp:
                        if str(resp).split()[0] == '550':
                            break
                    except ftplib.error_temp:
                        sleep(self.delay)
                        continue
                    except EOFError:
                        sleep(self.delay)
                        self.ftp = ftplib.FTP(
                            'ftp.ncbi.nlm.nih.gov', timeout=self.timeout)
                        self.ftp.login()
                        self.ftp.cwd('/genomes/all')
                        continue
                    else:
                        break
            if not success:
                print(f'WARNING: Cannot retrieve {fname}.')
                failures.append(g)
        print('Done.')

        # drop genomes that cannot be retrieved
        if len(failures):
            print('Failed to retrieve the following genomes:')
            print('  ' + ', '.join(failures))
            failures = set(failures)
            self.df = self.df[~self.df['genome'].isin(failures)]

    def extract_genomes(self):
        """Extract proteins from genomes.

        Notes
        -----
        Write protein sequences of all genomes into db.faa.
        Write protein to genomes(s) map to genome.map.gz.
        """
        print('Extracting downloaded genomic data...', end='', flush=True)
        dir_ = join(self.output, 'download', 'faa')
        prots = {}
        cp = None
        fout = open(join(self.output, 'db.faa'), 'w')

        def write_prot():
            if cp:
                fout.write(f'>{cp} {prots[cp]["name"]}\n{prots[cp]["seq"]}\n')

        g2n, g2aa = {}, {}
        for row in self.df.itertuples():
            g, tid = row.genome, row.taxid
            g2n[g], g2aa[g] = 0, 0
            stem = row.ftp_path.rsplit('/', 1)[-1]
            fp = join(dir_, f'{stem}_protein.faa.gz')
            with gzip.open(fp, 'rb') as f:
                try:
                    content = f.read().decode().splitlines()
                except (OSError, EOFError, TypeError):
                    print(f' skipping corrupted file {stem}.', end='',
                          flush=True)
                    continue
            cp = None
            for line in content:
                line = line.rstrip('\r\n')
                if line.startswith('>'):
                    write_prot()
                    p, name = line[1:].split(None, 1)
                    g2n[g] += 1
                    if p in prots:
                        cp = None
                        prots[p]['gs'].add(g)
                        prots[p]['tids'].add(tid)
                        g2aa[g] += len(prots[p]['seq'])
                    else:
                        cp = p
                        prots[p] = {
                            'name': name, 'gs': {g}, 'tids': {tid}, 'seq': ''}
                elif cp:
                    line = line.rstrip('*')
                    prots[cp]['seq'] += line
                    g2aa[g] += len(line)
            write_prot()
        fout.close()
        print(' done.')
        print('Combined protein sequences written to db.faa.')

        self.df['proteins'] = self.df['genome'].map(g2n)
        self.df['residues'] = self.df['genome'].map(g2aa)
        self.p2tids = {k: v['tids'] for k, v in prots.items()}

        # report summary
        print(f'Total number of genomes extracted: {len(g2n)}.')
        print(f'Total number of unique proteins extracted: {len(prots)}.')
        print('Number of proteins shared by multiple genomes: {}.'.format(
            sum([1 for k, v in prots.items() if len(v['gs']) > 1])))
        print('Number of proteins shared by multiple TaxIDs: {}.'.format(
            sum([1 for k, v in prots.items() if len(v['tids']) > 1])))

        # write protein to genomes(s) map
        fname = 'genome.map.gz'
        with gzip.open(join(self.output, fname), 'wb') as f:
            for k, v in sorted(prots.items()):
                f.write('{}\t{}\n'.format(k, ','.join(
                    sorted(v['gs']))).encode())
        print(f'Protein-to-genome(s) map written to {fname}.')

    def genome_lineages(self):
        """Generate lineage information for genomes.

        Notes
        -----
        Write genome lineages to lineages.txt.
        """
        # identify taxa at standard ranks
        ranks = ['superkingdom', 'kingdom', 'phylum', 'class', 'order',
                 'family', 'genus', 'species']
        self.df[ranks[:-1]] = self.df['taxid'].apply(lambda x: pd.Series(
            taxids_at_ranks(x, ranks[:-1], self.taxdump)))

        # report number of taxa represented at each rank
        print('Number of taxonomic groups represented:')
        for rank in ranks:
            print(f'  {rank}: {self.df[rank].nunique()}.')

        # merge superkingdom and kingdom
        self.df['kingdom'] = self.df[['superkingdom', 'kingdom']].apply(
            lambda x: x[1] if x[1] else x[0], axis=1)

        # generate lineage string
        self.df['lineage'] = self.df[ranks[1:]].fillna('').apply(
            lambda col: col.apply(lambda val, name: '{}__{}'.format(
                name[0], self.taxdump[val]['name'] if val else ''),
                args=(col.name,))).apply('; '.join, axis=1)

        # write table
        fname = 'lineages.txt'
        self.df['lineage'].to_csv(
            join(self.output, fname), sep='\t', header=False)
        print(f'Genome lineages written to {fname}.')

    def genome_metadata(self):
        """Write and report genome metadata.

        Notes
        -----
        Write genome metadata to genomes.tsv.
        """
        self.df.set_index('genome', inplace=True)
        self.df = self.df[[
            'proteins', 'residues', 'assembly_level',
            'accession', 'bioproject', 'biosample', 'asm_name',
            'organism_name', 'infraspecific_name', 'isolate', 'taxid',
            'ftp_path']]
        fname = 'genomes.tsv'
        self.df.to_csv(join(self.output, fname), sep='\t')
        print(f'Genome metadata written to {fname}.')

    def build_taxdump(self):
        """Build taxonomy database.

        Notes
        -----
        Write refined taxonomy database to taxdump/nodes.dmp and names.dmp.
        """
        tids = self.df['taxid'].unique().tolist()

        # identify all ancestral TaxIDs
        ancs = set(tids)
        for tid in tids:
            cid = tid
            while True:
                pid = self.taxdump[cid]['parent']
                if cid == pid or pid in ancs:
                    break
                ancs.add(pid)
                cid = pid

        # shrink taxdump files
        dir_ = join(self.output, 'taxdump')
        makedirs(dir_, exist_ok=True)
        for fname in ('nodes.dmp', 'names.dmp'):
            fo = open(join(dir_, fname), 'w')
            fi = open(join(self.tmpdir, fname), 'r')
            for line in fi:
                row = line.rstrip('\r\n').replace('\t|', '').split('\t')
                if row[0] in ancs:
                    if fname == 'nodes.dmp' or 'scientific name' in row:
                        fo.write(line)
            fi.close()
            fo.close()
        print('Taxonomy database written to taxdump/.')

    def build_taxonmap(self):
        """Build protein-to-taxonomy map.
        """
        # assign shared protein to lowest common ancestor (LCA)
        self.taxonmap = {p: max(tids) if len(tids) == 1 else find_lca(
            tids, self.taxdump) for p, tids in self.p2tids.items()}

        # write taxonomy map
        fname = 'taxon.map.gz'
        with gzip.open(join(self.output, fname), 'wb') as f:
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p}\t{tid}\n'.encode())
        print(f'Protein-to-taxonomy map written to {fname}.')

    def compile_database(self):
        """Compile database using external programs.
        """
        if self.compile == 'none':
            return
        if self.compile in ('blast', 'both'):
            self.build_blast_db()
        if self.compile in ('diamond', 'both'):
            self.build_diamond_db()

    def build_blast_db(self):
        """Build BLAST database.
        """
        # create temporary taxon map
        taxonmap = join(self.tmpdir, 'taxon.map')
        with open(taxonmap, 'w') as f:
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p}\t{tid}\n')

        # build BLAST database
        makedirs(join(self.output, 'blast'), exist_ok=True)
        cmd = ' '.join((
            self.makeblastdb,
            '-dbtype', 'prot',
            '-in', join(self.output, 'db.faa'),
            '-out', join(self.output, 'blast', 'db'),
            '-title', 'db',
            '-parse_seqids',
            '-taxid_map', taxonmap))
        print('Build BLAST database...', flush=True)
        ec, out = run_command(cmd)
        if ec:
            raise ValueError(f'makeblastdb failed with error code {ec}.')

        # clean up
        remove(taxonmap)
        print('Done.')

    def build_diamond_db(self):
        """Build DIAMOND database.
        """
        # create temporary taxon map
        taxonmap = join(self.tmpdir, 'prot.accession2taxid')
        with open(taxonmap, 'w') as f:
            f.write('accession\taccession.version\ttaxid\n')
            for p, tid in sorted(self.taxonmap.items()):
                f.write(f'{p.rsplit(".", 1)[0]}\t{p}\t{tid}\n')

        # build DIAMOND database
        makedirs(join(self.output, 'diamond'), exist_ok=True)
        cmd = ' '.join((
            self.diamond, 'makedb',
            '--threads', str(self.threads),
            '--in', join(self.output, 'db.faa'),
            '--taxonmap', taxonmap,
            '--taxonnodes', join(self.output, 'taxdump', 'nodes.dmp'),
            '--taxonnames', join(self.output, 'taxdump', 'names.dmp'),
            '--db', join(self.output, 'diamond', 'db'),
            '--tmpdir', self.tmpdir))
        print('Build DIAMOND database...', flush=True)
        ec, out = run_command(cmd)
        if ec:
            raise ValueError(f'diamond failed with error code {ec}.')

        # clean up
        remove(taxonmap)
        print('Done.')

    @staticmethod
    def check_local_file(file, overwrite=False):
        """Check existing local file.

        Parameters
        ----------
        file : str
            expected local file path and name
        overwrite : bool, optional
            whether overwrite existing file (default: False)

        Returns
        -------
        bool
            whether file exists
        """
        if isfile(file):
            if getsize(file) == 0:
                remove(file)
            elif overwrite:
                print(f'  WARNING: existing local file {basename(file)} will '
                      'be overwritten.')
                remove(file)
            else:
                print(f'  Using local file {basename(file)}.')
                return True
        return False

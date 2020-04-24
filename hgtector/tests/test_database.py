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
from os.path import join, dirname, realpath, isfile, isdir
from shutil import rmtree, copyfile
from tempfile import mkdtemp

from hgtector.database import Database


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

    def test_build_blast_db(self):
        me = Database()
        me.output = self.tmpdir
        me.makeblastdb = 'makeblastdb'
        copyfile(join(self.datadir, 'DnaK', 'linear.faa'),
                 join(self.tmpdir, 'db.faa'))
        copyfile(join(self.datadir, 'DnaK', 'prot2tid.txt'),
                 join(self.tmpdir, 'taxon.map'))
        me.build_blast_db()
        self.assertTrue(isdir(join(self.tmpdir, 'blast')))
        for ext in ('phr', 'pin', 'pog', 'psd', 'psi', 'psq'):
            self.assertTrue(isfile(join(self.tmpdir, 'blast', 'db.{}'.format(
                ext))))
        rmtree(join(self.tmpdir, 'blast'))
        remove(join(self.tmpdir, 'db.faa'))
        remove(join(self.tmpdir, 'taxon.map'))

    def test_build_diamond_db(self):
        me = Database()
        me.output = self.tmpdir
        me.diamond = 'diamond'
        me.threads = 1
        me.tmpdir = self.tmpdir
        copyfile(join(self.datadir, 'DnaK', 'linear.faa'),
                 join(self.tmpdir, 'db.faa'))
        with open(join(self.datadir, 'DnaK', 'prot2tid.txt'), 'r') as f:
            me.taxonmap = dict(x.split('\t') for x in f.read().splitlines())
        makedirs(join(self.tmpdir, 'taxdump'))
        copyfile(join(self.datadir, 'DnaK', 'taxdump', 'nodes.dmp'),
                 join(self.tmpdir, 'taxdump', 'nodes.dmp'))
        copyfile(join(self.datadir, 'DnaK', 'taxdump', 'names.dmp'),
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
        file = join(self.tmpdir, 'tmp.in')
        self.assertFalse(me.check_local_file(file))
        with open(file, 'w') as f:
            f.write('Hello world!')
        self.assertTrue(isfile(file))
        self.assertTrue(me.check_local_file(file))
        self.assertFalse(me.check_local_file(file, overwrite=True))
        self.assertFalse(isfile(file))


if __name__ == '__main__':
    main()

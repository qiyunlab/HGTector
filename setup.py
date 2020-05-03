#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, Qiyun Zhu and Katharina Dittmar.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from glob import glob

import hgtector.__init__ as init


classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS
    Operating System :: Windows
"""

params = {
    'name':             init.__name__,
    'version':          init.__version__,
    'license':          init.__license__,
    'description':      init.__description__,
    'long_description': open('README.md').read(),
    'long_description_content_type': 'text/markdown',
    'author':           init.__author__,
    'author_email':     init.__email__,
    'url':              init.__url__,
    'install_requires': ['pyyaml',
                         'numpy',
                         'scipy',
                         'pandas',
                         'scikit-learn',
                         'matplotlib'],
    'classifiers':      classes.strip().split('\n    '),
    'python_requires':  '>=3.6',
    'entry_points':     '[console_scripts]',
    'scripts':          glob(f'scripts/{init.__name__}'),
    'package_data':     {init.__name__: ['config.yml']},
    'include_package_data': True
}

setup(**params, packages=find_packages())

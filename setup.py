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

classifiers = [s.strip() for s in classes.split('\n') if s]


description = (
    'Genome-wide detection of putative horizontal gene transfer (HGT) events '
    'based on sequence homology search hit distribution statistics')

long_description = open('README.md').read()


setup(
    name='hgtector',
    version='2.0b1',
    license='BSD-3-Clause',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Qiyun Zhu',
    author_email='qiyunzhu@gmail.com',
    url='https://github.com/DittmarLab/HGTector',
    packages=find_packages(),
    scripts=glob('scripts/hgtector'),
    package_data={'hgtector': ['../config.yml']},
    include_package_data=True,
    install_requires=[
        'pyyaml',
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn',
        'matplotlib'
    ],
    classifiers=classifiers,
    python_requires='>=3.6',
    entry_points='[console_scripts]'
)

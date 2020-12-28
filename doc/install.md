Installation
============

## Operating system

There is no restriction as far as I am aware of. Tested and working on Linux, macOS and Windows systems.

(But Windows users will need extra configurations...)


## Software environment

HGTector is written in Python 3. One needs at least Python 3.6 to run the program. I recommend [Conda](https://docs.conda.io/en/latest/) for managing Python version and packages.


## Installation

### Option 1: Through Conda (recommended)

```bash
conda create -n hgtector python=3 pyyaml pandas matplotlib scikit-learn
conda activate hgtector
pip install git+https://github.com/qiyunlab/HGTector.git
```

### Option 2: Native installation

Download this [repository](https://github.com/qiyunlab/HGTector/archive/master.zip). Unzip. Then execute:

```bash
python setup.py install
```

Type `hgtector` to check if installation is successful, in which case command-line help information will be displayed on the screen.

You may now read [first run](1strun.md) and [second run](2ndrun.md) before proceeding with aligner and database installation.


## Aligner

One may use choice of [DIAMOND](https://github.com/bbuchfink/diamond) or [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) for sequence homology search. If you have already installed them, make sure they are callable from the environment, or use command-line [arguments](search.md#Local-search-behaviors) to point the executables to HGTector. Alternatively, you may install them via Conda:

```bash
conda install -c bioconda diamond blast
```


## Database

HGTector has a command `database` for automated database construction. It defaults to the **NCBI** RefSeq microbial genomes and taxonomy. Meanwhile, we also provide instructions for using **GTDB** and custom databases. See [details](database.md).

A standard database built using the default protocol on 2019-10-21 is available for [download](https://www.dropbox.com/s/qdnfgzdcjadlm4i/hgtdb_20191021.tar.xz?dl=0), together with [instruction](database.md#Manual-compiling) for compiling.

A small, pre-compiled test database is also available for [download](https://www.dropbox.com/s/46v3uc708rvc5rc/ref107.tar.xz?dl=0).


## Upgrade

Just add `--upgrade` or `-U` to the pip command:

```bash
pip install -U git+https://github.com/qiyunlab/HGTector.git
```

Note: You can only upgrade from HGTector 2.0b1 or above. You cannot upgrade from older versions, which were written in Perl.


## Uninstallation

```bash
pip uninstall hgtector
```

If you no longer need the conda environment:

```bash
conda env remove -n hgtector --all
```


## Compatibility

If in the future some dependencies have changes that are not compatible with the current release of HGTector, the following "safe" command can be used to install the current versions of dependencies (note: DIAMOND version is too tricky to specify).

```bash
conda create -n hgtector -c bioconda python=3.7.7 pyyaml=5.3.1 pandas=1.0.3 matplotlib=3.1.3 scikit-learn=0.22.1 diamond blast=2.9.0
```

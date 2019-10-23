Installation
============

## Operating system

There is no restriction as far as I am aware of. Tested and working on Linux, macOS and Windows systems.

## Software environment

HGTector is written in Python 3. One needs at least Python 3.6 to run the program. I recommend [Conda](https://docs.conda.io/en/latest/) for managing Python version and packages.

## Installation

### Option 1: Through Conda (recommended)

```bash
conda create -n hgtector python=3 pyyaml pandas matplotlib scikit-learn
conda activate hgtector
pip install git+https://github.com/DittmarLab/HGTector.git
```

### Option 2: Native installtion

Download this [repository](https://github.com/DittmarLab/HGTector/archive/master.zip). Unzip. Then execute:

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

HGTector has a command `database` for automated database construction. See [details](database.md).

A standard database built using the default protocol on 2019-11-22 is available for [download](https://www.dropbox.com/s/qdnfgzdcjadlm4i/hgtdb_20191021.tar.xz?dl=0), together with [instruction](database.md#Manual-compiling) for compiling.

A small, pre-compiled test database is also available for [download](https://www.dropbox.com/s/46v3uc708rvc5rc/ref107.tar.xz?dl=0).

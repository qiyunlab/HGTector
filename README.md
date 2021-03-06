HGTector2
=========

The development of HGTector is now at [qiyunlab](https://qiyunlab.github.io/). Versions starting from 2.0b3 will be released from this repo. Please access HGTector using the new URL: https://github.com/qiyunlab/HGTector.

**HGTector2** is a completely re-engineered software tool, featuring a fully automated analytical pipeline with smart determination of parameters which requires minimum human involvement, a re-designed command-line interface which facilitates standardized scientific computing, and a high-quality Python 3 codebase.

**HGTector** is a computational pipeline for genome-wide detection of putative horizontal gene transfer (HGT) events based on sequence homology search hit distribution statistics.

## Documentation

[What's New](CHANGELOG.md)

[Installation](doc/install.md)

Tutorials
- [First Run](doc/1strun.md)
- [Second Run](doc/2ndrun.md)
- [Real Runs](doc/realrun.md)

References
- [Search](doc/search.md)
- [Analyze](doc/analyze.md)
- [Database](doc/database.md)
- [Configure](doc/config.md)


## Quick start

Set up a Conda environment and install dependencies:

```bash
conda create -n hgtector python=3 pyyaml pandas matplotlib scikit-learn bioconda::diamond
conda activate hgtector
```

Install HGTector2:

```bash
pip install git+https://github.com/qiyunlab/HGTector.git
```

Build a reference database using the default protocol:

```bash
hgtector database -o db_dir --default
```

This will retrieve the latest genomic data from NCBI. If this does not work (e.g., due to network issues), or you need some customization, please read the [database](doc/database.md) page.

Prepare input file(s). They should be multi-Fasta files of amino acid sequences (faa). Each file represents the whole protein set of a complete or partial genome.

Perform homology [search](doc/search.md):

```bash
hgtector search -i input.faa -o search_dir -m diamond -p 16 -d db_dir/diamond/db -t db_dir/taxdump
```

Perform HGT [prediction](doc/analyze.md):

```bash
hgtector analyze -i search_dir -o analyze_dir -t hgtdb/taxdump
```

Examine the prediction results under the `analyze_dir` directory.

It is recommended that you read the [first run](doc/1strun.md), [second run](doc/2ndrun.md) and [real runs](doc/realrun.md) pages to get familiar with the pipeline, the underlying methodology, and the customization options.


## License

Copyright (c) 2013-2020, [Qiyun Zhu](mailto:qiyunzhu@gmail.com) and [Katharina Dittmar](mailto:katharinad@gmail.com). Licensed under [BSD 3-clause](http://opensource.org/licenses/BSD-3-Clause). See full license [statement](LICENSE).


## Citation

> Zhu Q, Kosoy M, Dittmar K. HGTector: [an automated method facilitating genome-wide discovery of putative horizontal gene transfers](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-15-717). *BMC Genomics*. 2014. 15:717.

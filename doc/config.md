Configuration
=============

## Overview

The default parameter settings are defined in a configuration file, `config.yml`. The original file is located in the program directory. The file is in [YAML](https://en.wikipedia.org/wiki/YAML) format, which is human-readable and editable. For example:

```yaml
## Database locations
database:
  # reference protein sequence database for DIAMOND
  diamond: /home/me/dbs/hgtector/2019-11-22/diamond/db.dmnd

## Sequence homology search
search:
  # search cutoffs
  evalue: 1.0e-20
```

One may directly modify it using a text editor. But don't forget to make a backup before saving!

- A typical use case is to enter the paths to the databases in the configuration file, so that one does not need to specify them in future uses.

HGTector will sequentially look for this file in the following locations:

1. Current directory (`.`)
2. Home directory (`~`), under subdirectory `.hgtector`
3. Program directory

Therefore one may prepare analysis-specific configuration files, if the same settings are to be executed multiple times and records are needed.

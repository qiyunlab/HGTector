Analyze
=======

## Overview

The `analyze` command predicts HGT-derived genes based on the distribution of the patterns of sequence homology search results of all proteins (genes) in one or more samples (genomes). Several statistical approaches are involved in this workflow.

```bash
hgtector analyze -i sample.tsv -o <output_dir> <parameters...>
```

The underlying procedures and the guideline for tuning parameters are explained in detail in [second run](2ndrun.md).

### Input

The input may be a file or a directory of files which are the output of the `search` command.

### Output

File | Description
--- | ---
`scores.tsv` | A single multi-Fasta file containing all protein sequences.
`hgts/<sample>.txt` | Predicted HGT-derived genes and their silhouette scores.
`<group>.hist.png` | Histogram of scores of a group.
`<group>.kde.png` | Density function of scores of a group, and clustering threshold (grey line).
`scatter.png` | Scatter plot of distal vs. close scores. Each point represents a gene (protein).

Format of `scores.tsv`:

Field | Description
--- | ---
`sample` | Sample ID
`protein` | Protein ID
`length` | Protein length (aa)
`hits` | Total number of hits after filtering
`self`, `close`, `distal` | Score (sum of normalized bit scores) of each group
`match` | Best match in "distal" group which implicates potential donor for predicted HGTs


## Command-line reference

### Basic

Option | Default | Description
--- | --- | ---
`-i`, `--input` | - | Input search result file, or directory where one or more input files are located.
`-o`, `--output` | - | Directory where analysis results will be saved.
`-t`, `--taxdump` | - | Directory of taxonomy database files (`nodes.dmp` and `names.dmp`). Required if they are not found in input directory.

### Hit filtering

Option | Default | Description
--- | --- | ---
`-k`, `--maxhits` | - | Maximum number of sequence similarity search hits per gene (protein) to preserve.
`--evalue` | - | Maximum E-value cutoff.
`--identity` | - | Minimum percent identity cutoff.
`--coverage` | - | Minimum percent query coverage cutoff.

### Taxonomic assignment

Option | Default | Description
--- | --- | ---
`--input-tax` | - | TaxIDs of input samples. Will auto-infer if omitted. Can be comma-delimited, colon-inserted string (`sample1:taxId1,sample2:taxId2...`), or a file in the format of `sample <tab> taxId`.
`--input-cov` | 75 | For auto-inference: an input sample will be assigned to a taxon if it represents this percentage or more best hits.

### Taxonomic grouping

Option | Default | Description
--- | --- | ---
`--self-tax` | - | TaxIDs of "self" group (a comma-delimited string, or a file of one TaxID per line). Will auto-infer if omitted.
`--close-tax` | - | TaxIDs of "close" group (a comma-delimited string, or a file of one TaxID per line). Will auto-infer if omitted.
`--self-rank` | - | For auto-inference: "self" group must be at or above this rank (e.g., species, genus, family...).
`--close-size` | 10 | For auto-inference: "close" group must have at least this number of taxa.
`--distal-top` | 10 | Find a match in "distal" group which is LCA of hits with bit score at most this percentage lower than the best hit. The behavior is consistent with DIAMOND's `--top` parameter. This match implicates the potential donor of an HGT-derived gene.


### Scoring

Option | Default | Description
--- | --- | ---
`--weighted` | yes | Score is sum of weighted bit scores; otherwise simple counts.

### Filtering

Option | Default | Description
--- | --- | ---
`--outliers` | zscore | Detect and remove outliers using selected method (none, zscore, boxplot).
`--orphans` | no | Keep orphans (proteins without non-self hits) in statistical analysis.

### Prediction

Option | Default | Description
--- | --- | ---
`--bandwidth` | auto | Bandwidth for Gaussian KDE (auto, grid, silverman, or a number between 0.1 and 1.0).
`--bw-steps` | 20 | Number of steps for auto and grid kernel bandwidth optimization.
`--low-part` | 75 | Maximum percentage below threshold for automatic bandwidth optimization.
`--noise` | 50 | Percent valley-to-peak distance to exclude from cluster.
`--fixed` | 25 | Use this percentage as threshold if KDE clustering fails.
`--silhouette` | 0.5 | Silhouette score threshold for cluster refinement.
`--self-low` | no | HGT has low "self" score (an optional criterion).

### Program behavior

Option | Default | Description
--- | --- | ---
`--from-scores` | - | If score table already exists, use it and skip search result parsing and taxonomy inference. Otherwise overwrite it.

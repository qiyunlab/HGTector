Search
======

## Overview

The `search` command performs batch homology searching on input protein sequences against a reference sequence database, and filters hits based on alignment metrics and taxonomic information.

```bash
hgtector search -i sample.fa -o <output_dir> <parameters...>
```

### Remote search

The simplest (but not recommended) usage of the `search` command is to have it communicate with the NCBI server to complete [searching](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins), [self-aligning](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&BLAST_SPEC=blast2seq&LINK_LOC=align2seq), [sequence](https://www.ncbi.nlm.nih.gov/protein) and [taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) information retrieval. This does not require any local database or compute resource. The command is as simple as:

```bash
hgtector search -i sample.fa -o .
```

### Local search

In this recommended mode, the program requires an external aligner to execute searching. Two common aligners are supported: [**DIAMOND**](https://github.com/bbuchfink/diamond) ([Buchfink et al., 2015](https://www.nature.com/articles/nmeth.3176)) and [**BLASTp**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) ([Altschul et al., 1990](https://www.sciencedirect.com/science/article/pii/S0022283605803602)). We recommend DIAMOND for normal use cases, since it is much faster whereas retaining comparable accuracy versus the classical BLAST.

The program also needs a reference protein sequence database, and a taxonomy database. Please refer to the [database](database.md) section.

```bash
hgtector search -i sample.faa.gz -o . -m diamond -p 32 -d <diamond_db> -t <taxdump_dir>
```

### Output file

One output file, `sample.tsv` is created for each sample to host basic properties and hit table for each protein. The header section contains the following keys:

```
ID, Length, Product, Score, Hits
```

The hit table has the following columns, which are a subset of the [standard BLAST output metrics](https://www.ncbi.nlm.nih.gov/books/NBK279684/):

Code | Description
--- | ---
`qaccver` | Query accesion.version
`saccver` | Subject accession.version
`pident` | Percentage of identical matches
`evalue` | Expect value
`bitscore` | Bit score
`qcovhsp` | Query coverage per HSP
`staxids` | Unique subject taxonomy ID(s)

### Search thresholds

Multiple search parameters can be specified (see below for details). Here is an example command:

```bash
hgtector search ... --maxhits 500 --minsize 50 --evalue 1e-20 --identity 50 --coverage 50 --maxseqs 1000 --extrargs "--sensitive --range-culling"
```

Note that `maxhits` and `maxseqs` are different. The later controls how many hits to return in a DIAMOND/BLAST search; the former controls how many hits to preserve after filtering.

### Taxonomic filtering

Multiple parameters are available to control the inclusion and exclusion of taxonomic groups. For example:

```bash
hgtector search ... --tax-include 2,2157 --tax-exclude 1117,766 --tax-unirank genus --tax-block environmental,uncultured,...
```

This will limit the search result within domains Bacteria (TaxID: [2](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=2)) and Archaea (TaxID: [2157](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=2157)), but exclude Cyanobacteria (TaxID: [1117](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1117)) and Rickettsiales (TaxID: [766](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=780)). It will only retain the top hit within each genus. It will also remove hits with taxon names reading "environmental XXX", "uncultured XXX", etc.



### Pre-computed search result

In some cases, sequence homology search has already been executed for other tasks in the research project, and it is favorable to re-use the search result to ensure consistency and to save compute.

```bash
hgtector search -i sample.faa -m precomp -s sample.blast6out -o . -t <taxdump_dir> [--taxmap taxon.map.gz]
```

One still needs the taxonomy database. If the pre-computed search result does not contain subject TaxIDs, an additional taxon map is also needed. See [database](database.md).

The native HGTector-style hit table can be manually generated using the following commands:

BLAST:

```bash
blastp -query sample.fa -db <blast_db> -outfmt "6 qaccver saccver pident evalue bitscore qcovhsp staxids"
```

DIAMOND:

```bash
diamond blastp --query sample.fa --db <diamond_db> --outfmt 6 qseqid sseqid pident evalue bitscore qcovhsp staxids
```

But this is not the typical command one would execute before meeting HGTector. In many cases, the standard tabular output (`-outfmt 6`) is used:

```
diamond blastp --query sample.fa --db <diamond_db> --outfmt 6
blastp -query sample.fa -db <blast_db> -outfmt 6
```

HGTector automatically recognizes and parses this format. Note that subject TaxID is not part of the standard tabular output, so a taxon map file is necessary in this case.

### Self-alignment

HGTector computes the bit score of each protein sequence aligned against itself, and uses it as a baseline for normalizing bits scores of subjects. Several methods are available for doing this self-alignment, and the priority is automatically determined by default. The typical scenario is as follow (in case one wants to override):

```bash
hgtector search ... -m diamond --aln-method native
```

This will use the same method (DIAMOND) for both homology search and self-alignment.

Other options are:

`lookup`: The search result already contains a self-hit. This only works if the query protein is also in the database.

`fast`: Use a built-in implementation of a simplified BLAST algorithm to calculate this score. The result should be identical to the result of DIAMOND in its default cost scoring setting, but it is slightly different from the result of modern releases of BLAST. If the native method fails (it happens occassionally), the program will automatically use this algorithm instead.

`precomp`: Provide precomputed self-alignment results, in a simple format of `protein <tab> score`. For example:

```bash
hgtector search ... -m precomp -s sample.b6 --aln-method precomp --aln-precomp sample.txt
```

### Multiple samples

One can perform search on multiple samples under the same directory:

```bash
hgtector search -i <input_dir> -o <output_dir> [-s <precompute_dir>] <parameters...>
```

This will generate one output file per input file. The filename without extension (e.g., "Bin1" of `Bin1.faa`) will be considered as sample ID and carried to the output file (`Bin1.tsv`).


## Command-line reference

### Basic

Option | Default | Description
--- | --- | ---
`-i`, `--input` | - | A file or a directory of files which store query protein sequences in multi-Fasta format. Compressed files (`.gz`, `.bz2`, `.xz` and `.lz`) are supported. Plain lists of protein IDs without sequences are supported if the IDs can be found in the database.
`-o`, `--output` | - | Directory where search results are to be saved.
`-m`, `--method` | auto | Search method (auto, diamond, blast, remote, precomp).
`-s`, `--precomp` | - | File or directory of precomputed search results (when method = precomp).

### Database

Option | Default | Description
--- | --- | ---
`-d`, `--db` | - | Reference protein sequence database.
`-t`, `--taxdump` | - | Directory of taxonomy database files (`nodes.dmp` and `names.dmp`).
`--taxmap` | - | Sequence Id to taxId mapping file (not necessary if protein database already contains taxonomy).

### Search behaviors

Option | Default | Description
--- | --- | ---
`-k`, `--maxhits` | 0 | Maximum number of hits to preserve per query (0 for unlimited)
`--minsize` | 30 | Minimum length of query sequence (aa).
`--queries` | (depends) | Number of queries per run (0 for whole sample).
`--maxchars` | (depends) | Maximum number of characters per run (0 for unlimited).

### Search cutoffs

Option | Default | Description
--- | --- | ---
`--maxseqs` | 500 | Maximum number of sequences to return.
`--evalue` | 1e-5 | Maximum E-value cutoff.
`--identity` | 0 | Minimum percent identity cutoff.
`--coverage` | 0 | Minimum percent query coverage cutoff.
`--extrargs` | (depends) | Extra arguments for choice of search method. This can be command-line arguments for DIAMOND or BLAST, or URL API for remote search.

### Taxonomic filters

Option | Default | Description
--- | --- | ---
`--tax-include` | - | Include taxa under those TaxIDs (a comma-delimited string, or a file of one TaxID per line).
`--tax-exclude` | - | Exclude taxa under those TaxIDs (a comma-delimited string, or a file of one TaxID per line).
`--tax-unique` | yes | Ignore more than one hit with same TaxId.
`--tax-unirank` | - | Ignore more than one hit under same taxon at this rank (species, genus, family, etc.). Recommended if the database is highly taxonomically imbalanced (e.g., contains thousands of _E. coli_ strains).
`--tax-capital` | yes | Ignore taxon names that are not capitalized.
`--tax-latin` | no | Ignore species names that are not Latinate.
`--tax-block` | unknown,<br>uncultured,<br>unidentified,<br>unclassified,<br>unresolved,<br>environmental,<br>plasmid,<br>vector,<br>synthetic,<br>phage | Ignore taxon names containing any of these words.

### Local search behaviors

Option | Default | Description
--- | --- | ---
`-p`, `--threads` | 0 | Number of threads (0 for all CPU cores).
`--tmpdir` | - | Temporary directory.
`--diamond` | diamond | diamond executable.
`--blastp` | blastp | blastp executable.
`--blastdbcmd` | blastdbcmd | blastdbcmd executable.

### Remote search behaviors

Option | Default | Description
--- | --- | ---
`--algorithm` | kmerBlastp | Remote search algorithm (blastp, psiBlast, deltaBlast, kmerBlastp (a.k.a. Quick BLASTP), phiBlast, etc., see [here](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins)).
`--retries` | 5 | Maximum number of retries per search.
`--delay` | 60 | Seconds between two search requests.
`--timeout` | 900 | Seconds before program gives up waiting.
`--entrez` | all [filter] NOT(environmental samples[filter] OR metagenomes[orgn]) txid131567[orgn] | Entrez query text.
`--server` | [link](https://blast.ncbi.nlm.nih.gov/Blast.cgi) | Remote search server URL.

### Self-alignment options

Option | Default | Description
--- | --- | ---
`--aln-method` | auto | Self-alignment method (auto, native, fast, lookup, precomp).
`--aln-precomp` | - | File or directory of precomputed sequence Id to score maps (when self-alignment method = precomp).
`--aln-server` | [link](https://blast.ncbi.nlm.nih.gov/BlastAlign.cgi) | Remote align server URL.

### Remote fetch options

Option | Default | Description
--- | --- | ---
`--fetch-enable` | auto | Whether to enable remote fetch (yes, no, auto).
`--fetch-queries` | 100 | Maximum number of query entries per search.
`--fetch-retries` | 3 | Maximum number of retries per search.
`--fetch-delay` | 5 | Seconds between two fetch requests.
`--fetch-timeout` | 60 | Seconds before program gives up waiting.
`--fetch-server` | [link](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi) | Remote fetch server URL.

Change Log
==========

## Version 2.0b3 (2021, ongoing)

### Added
- Predicted HGT list now includes potential donors. Users can optionally specify a taxonomic rank at which they will be reported.
- A quick-start guide added to the homepage.
- Let the database construction task resume from an interrupted task and use already downloaded files, while detecting and skipping corrupt files.
- Added an option `--above` to sample genomes at all ranks from the designated one to phylum during database construction.
- Added an option `--typemater` to include NCBI-defined type materials in genome sampling during database construction.
- Added an option `--manual` to export URLs of sampled genomes during database construction, and let the user download them manually.

### Changed
- Updated pre-built default database to 2021-11-21 (after NCBI RefSeq 209)
- Repository transferred from [DittmarLab](https://github.com/DittmarLab) to [qiyunlab](https://github.com/qiyunlab).
- Updated recommended dependency versions, however the program should continue to be compatible with previous versions.
- Minor tweaks with no visible impact on program behavior.

### Fixed
- Fixed an issue with the NCBI FTP server connection during database construction. NCBI now recommends rsync over ftp. Therefore the protocol has been updated accordingly.
- Fixed compatibility with latest DIAMOND (2.0.13).


## Version 2.0b2 (5/3/2020)

Note: The command-line interface and the format and content of input, output and database files are identical to last version (2.0b1). You may safely combine analysis results of this and last versions without the need for re-analysis.

### Changed
- Increased Python unit testing coverage.
- Re-styled Python code to adopt recent features (3.6+).
- Updated installation tutorial, script and dependency versions.
- Minor tweaks of code.

### Added
- Instructions for using the GTDB database.

### Fixed
- Removed deprecated parameter "iid" for grid search since scikit-learn 0.22.
- Fixed a bug in custom genome exclusion during database building.
- Fixed a bug which breaks the program when the temporary directory (tmpdir) is not specified.
- Allowed matplotlib to generate figures when it cannot open X window (e.g., in a remote server).


## Version 2.0b1 (complete rework) (10/23/2019)

### Changed
- Completely re-written the program using Python 3, following modern software engineering practices.
- Changed license to BSD 3-clause.
- Combined search results into a single file to improve disk I/O performance.
- New format of taxonomy database which is identical to NCBI taxdump.
- Relaxed from dependency on a sequence-to-taxon mapping file to improve disk and memory performance.
- Re-designed module organization. Now only two modules: "search" and "analyze" are required for a typical analysis.

### Added
- Standardized installation method using pip / Conda.
- Standardized command-line user interface.
- Complete documentation in GitHub markdown format.
- Automated decision of search strategy.
  - Support for extra search arguments.
- Support for pre-computed search results.
- Automated self-alignment workflow.
  - Built-in algorithm for calculating self-alignment score.
- Automated inference of input genome taxonomy.
- Automated design of grouping scenario.
- Automated optimization of kernel bandwidth.
  - Built-in algorithm for progressive kernel bandwidth calculation.
  - Grid search with cross validation for bandwidth optimization.
- Automated refinement of prediction result using silhouette scores in 2D space.
- Database download workflow can resume from an interrupted run.
- New plotting features.

### Fixed
- NCBI FTP server timeout issue.
- Compatible with latest DIAMOND (0.9.26).
- Compatible with latest NCBI BLASTp and eFetch servers.

### Removed
- Dependency on R
- Modules: reporter, treer, orthologer.
- HTML user interface.
- Step-by-step wizard.
- RAPSearch2 support.
- Dip test for unimodality.
- Histogram-based clustering.


## Version 0.2.2 (8/23/2017)

### Added
- Provided a pre-built standard database.
- Added an installation script installer.sh.
- Adopted new NCBI standard, i.e., using accession instead of GI as sequence identifier.
- Modified the mechanism of databaser.py. Now it downloads both DNA and protein sequences.
- Replaced 'http' with 'https' in NCBI server URLs.
- Adopted the command-line interface of DIAMOND 0.8+.

### Fixed
- Reformatted all Perl and Python scripts.
- Reworded multiple screen prompts.
- Multiple trivial bug fixes.


## Incremental update (1/28/2017)

### Fixed
- Fixed a bug in databaser.py which incorrectly deals the case with subsampling off.


## Incremental update (11/6/2016)

### Added
- Added Python 3 support for databaser.py

### Fixed
- NCBI has changed the FTP site file system structure (see [notification](https://www.ncbi.nlm.nih.gov/news/08-30-2016-genomes-ftp-reorganization/)). So databaser.py was modified accordingly.


## Version 0.2.1 (2/2/2016)

### Added
- Optimized the process of self-alignment by DIAMOND and RAPSearch2. Significantly faster. Results don't change.
- Added a smart decision-making process. If the user does not specify a proper grouping scenario, the program will try to generate one.
- Changed the way of creating graphics.

### Fixed
- Disabled FTP passive mode in databaser.py.


## Version 0.2.0 (major upgrade) (10/25/2015)

### Changed
- Major modification of the entire pipeline.

Structural change:
- Renamed blaster and summarizer as searcher and reporter, respectively.
- Removed filterer, taxonomer and validator.
- Renamed input / output files and directories.
- Re-wrote searcher.
- Re-wrote the main program.
- Re-wrote the GUI.

### Added
- New script databaser.py to create a customized, taxonomically balanced BLAST database, which significantly improves the efficiency and accuracy of HGT discovery, comparing to the original nr.
- Added native support for non-BLAST search tools, including RAPSearch2 and DIAMOND.
- Allowed one protein to correspond to multiple TaxIDs. This further balances the taxonomic distribution of reference proteins.

### Fixed
- Multiple bug fixes.


## Incremental update (9/1/2015)

### Added
- Re-wrote the GUI.

### Fixed
- Several bug fixes.


## Version 0.1.9 (6/1/2015)

### Added
- Accepting pre-computed BLAST results (by BLAST, Blat, Usearch, etc).
- Sending multiple query sequences per search, thus accelerating the batch BLAST task.
- Looking up taxonomy from a GI-to-TaxID dictionary, thus accelerating the process.
- Plotting reference gene set along with prediction results.
- Percentage coverage cutoff (for BLAST search).

### Fixed
- Updated the BLAST server URL.
- Improved efficiency of batch BLAST and result processing.
- Allowed taxonomy of input genomes defined in config.txt.
- Multiple bug fixes in Perl code and in GUI.

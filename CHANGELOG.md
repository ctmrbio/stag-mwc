# Changelog
All notable changes to this project will be documented in this file.

The format is inspired by [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project loosely adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).  
The version numbering scheme consists of three numbers separated by dots:
`major.minor.patch`. Major versions are incremented when a substantial change
in overall functionality of StaG-mwc is introduced. Minor versions are
incremented for any modifications to the interface or output files (i.e.
changes that would likely lead to different output for end-users when
re-running an analysis, either by giving error messages or changing output
files), and the patch version is typically incremented for any set of changes
committed to the master branch that does not trigger any of the aforementioned
situations.

## [Unreleased]
### Added
- Added CHANGELOG.md
- New functionality to run mappers several times against different databases,
  based on a list of reference databases to map against in the config file.
- Functional profiling using HUMAnN2. Will automatically run all
  MetaPhlAn2-associated rules to produce taxonomic profiles for use in HUMAnN2.
- Added Overview page to documentation that includes a draft of a simplified
  graph overview of the workflow (including some unfinished parts).
- Added rules to run Kraken2.

### Changed
- Substantial improvements to Rackham Slurm profile, focusing on better Slurm
  log handling.
- A few low-impact rules that can be run locally are now declared as localrules.
- Replaced MEGARes antibiotic resistance gene mapping with GROOT resistance
  gene profiling using gene variation graphs, using a default database based on
  arg-annot.
- Added clustered sketch comparison output heatmap.
- Updated MetaPhlAn2 to version 2.7.8, with corresponding changes to config file.

### Fixed
- Fixed error handling if hg19 database is missing for the remove human step.

### Removed
- Removed duplicated MetaPhlAn2 configuration parameters from HUMAnN2 config
  section.


## [0.1.1-dev] - 2018-04-30
### Added

### Changed
- Started using Python's pathlib module for Snakefile rule input, output, and
  log file declarations. Some unsightly explicit string conversions still remain,
  due to Snakemake not being fully compatible with pathlib (yet).
- Add details about branching structure/strategy to CONTRIBUTING.md

### Removed


## [0.1.0-dev] - 2018-04-30 
First public release

### Added
- First public release of StaG-mwc! 
- Snakemake workflow capable of read preprocessing, rudimentary sequencing
  depth assessment using kmer uniqueness counting, naive sample comparison using
  MinHash sketches, mapping to user-defined databases using BBMap or Bowtie2
  (with customizable read count summarization per annotated feature), taxonomic
  profiling using Centrifuge, Kaiju, or MetaPhlAn2, and basics required for
  antibiotic resistance gene detection using MEGARes.
- First public draft of docs, available at https://stag-mwc.readthedocs.org.

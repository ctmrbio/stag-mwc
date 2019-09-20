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

## [0.3.3] Unreleased
### Added
- Added resource limiter for HUMAnN2 due to its intense use of huge temporary
  files in the output folder. Activated with --resources humann2=X, where X is
  the max number of parallel instances of humann2 to run.
- Added extra argument to BBDuk for additional quality control flexibility.
- Added groot report parameters `covcutoff` and `lowcov` to config file.
- Added second FastQC run after quality trimming.
- Added automatic plot of proportion human reads.
- Added MEGAHIT assembly step.

### Fixed
- Fixed bug in Slurm profile handling of cancelled/failed jobs.

### Changed
- Added read length window filter before groot alignment step.
- Change logdir of remove_human rule to LOGDIR/remove_human instead of
  OUTDIR/logs/remove_human.
- Improved make_count_table.py so it can use TSV annotation files with multiple
  columns. Added config setting for which columns to include.
- Updated GROOT to v0.8.4.
- Cleaned up sketch comparison cluster heatmap plotting script, making it more 
  robust to variations in output from different BBTools versions.
- Updated MetaPhlAn2 to 2.9.12 from conda. This required adding a local version
  of metaphlan_hclust_heatmap.py as that has disappeared in recent conda version. 
  Also required changing the call of merge_metaphlan_tables.py due to undocumented
  CLI change in conda version.


## [0.3.2-dev]
### Added
- Added Zenodo DOI reference. https://zenodo.org/badge/latestdoi/125840716
- Add printout of citations of used tools after successful workflow completion.

### Changed
- Updated docs regarding HTML execution report generation.
- Updated GROOT to 0.8.3


## [0.3.1]
### Changed
- Fixed MinHash sketch sample similarity plots.


## [0.3.0]
### Added
- Added CHANGELOG.md
- New functionality to run mappers several times against different databases,
  based on a list of reference databases to map against in the config file.
- Functional profiling using HUMAnN2. Will automatically run all
  MetaPhlAn2-associated rules to produce taxonomic profiles for use in HUMAnN2.
- Added Overview page to documentation that includes a draft of a simplified
  graph overview of the workflow (including some unfinished parts).
- Added rules to run Kraken2.
- Added onstart, onerror, and onsuccess messages.
- Added `email` functionality. The workflow can now automatically send an email
  after a successful or failed run.
- Added automatic report generation upon successful workflow completion.

### Changed
- Substantial improvements to Rackham Slurm profile, focusing on better Slurm
  log handling.
- A few low-impact rules that can be run locally are now declared as localrules.
- Replaced MEGARes antibiotic resistance gene mapping with Groot resistance gene
  profiling using gene variation graphs.
- Increased resource requirements for remove_human step in Rackham cluster profile.
- Added clustered sketch comparison output heatmap.
- Updated MetaPhlAn2 to version 2.7.8, with corresponding changes to config file.

### Fixed
- Fixed error handling if hg19 database is missing for the remove human step.

### Removed


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

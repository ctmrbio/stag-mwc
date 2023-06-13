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

## [0.7.1] Unreleased
### Added

### Fixed

### Changed
  
### Deprecated

### Removed


## [0.7.0] 2023-06-13
### Added
- Host removal: Bowtie2 now available as an option for host removal.

### Fixed
- HUMAnN3: Fixed critical bug leading to entire system-wide temporary directory
  being emptied after successfull HUMAnN run.
- Singularity: All Singularity definition files should from now on get version
  bumps in the description labels when conda environments built inside them are
  updated to reduce the risk of Singularity reusing old cached copies of images
  instead of download the latest version.

### Changed
- Preprocessing summary: Preprocessing summary script can now output a table of
  read counts regardless of which combination of read QC and host removal is
  used.
  
### Deprecated

### Removed
- AMR++, Groot: All tools for antibiotic resistance gene profiling have been
  removed entirely because they were out of date and had few active users.
  Users wanting to perform antibiotic resistance gene profiling are suggested
  to use the mapper modules with a suitable reference database or run the
  latest version of AMR++ separately.
- Assembly: The MetaWRAP-based assembly parts of the workflow have been removed
  entirely.


## [0.6.1] 2023-06-01
### Added
- BBMap: now outputs sorted BAM file, added options `keep_sam` and `keep_bam`.
- Bowtie2: added option `keep_bam`.

### Fixed
- KrakenUniq: environment variable `LC_ALL` has been added to Singularity image
  to prevent unnecessary warning messages related to it being undefined.
- KrakenUniq: now able to run when host removal is skipped, solved by adding
  `krakenuniq_merge_reads` rule to create a temporary merged fasta file with
  input data for KrakenUniq to avoid giving KrakenUniq symlinks as input.

### Changed
- MetaPhlAn: Updated to v4.0.6
- HUMAnN3: Updated to v3.7
- HUMAnN3: Changed the way the temporary directory is resolved, now using
  Snakemake's built-in `resources.tmpdir`. This should prevent HUMAnN from
  creating large temporary directories outside of Slurm job folders so that
  they cannot be automatically cleaned up if the Slurm job times out or fails
  before HUMAnN can clean up after itself.
- KrakenUniq: Concatenate reads with BBMap's `fuse.sh` with a padding of one
  `N` instead of interleaving the paired inputs into a single FASTA to avoid
  KrakenUniq treating paired reads independently.

### Deprecated
- Kaiju, Kraken2, MetaPhlAn: area plot removed due to repeatedly leading to
  failed runs in cached Singularity containers. The script still works as
  intended in newer matplotlib versions and will remain in the scripts folder
  for potential manual use if desired.

### Removed
- Groot: Removed settings related to read length window as that feature was
  removed in a previous StaG release.


## [0.6.0] 2023-04-17
### Added
- Added a new Slurm profile for use on CTMR Gandalf, also intended to be useful
  as a starting point for creating custom Slurm profiles.
- Added a README with basic instructions for how to configure the workflow.
- Added function to disable MetaPhlAn heatmap plots, which may be useful when
  processing very large numbers of samples.
- Added MetaPhlAn-style output tables for KrakenUniq.
- Added Krona plot output for KrakenUniq.

### Fixed
- Fixed missing interactive Kaiju Krona plots for all samples in final report.
- Now reusing metaphlan conda environment and biobakery container for running
  bowtie2 mapping rule to ensure a consistent execution environment for bowtie2.
- KrakenUniq now works in Singularity, thanks to new custom Singularity image.
- KrakenUniq rule is now correctly not rerun if `keep_kraken` or `keep_kreport`
  settings are set to false when executing the workflow a second time.
- Kraken2 rule is now correctly not rerun if `keep_kraken` or `keep_kreport`
  settings are set to false when executing the workflow a second time.

### Changed
- Restructured repo to conform to modern Snakemake best practices. This also
  includes updates to documentation where needed.
- Hardcoded default thread values for all rules used during local execution
  without profile. Intended to be overridden by profile.
- Updated KrakenTools to its latest version (1.2), with a minor custom modification
  of `kreport2mpa.py`, changing output column names to `taxon_name` and `reads`.
- Updated all Python packages in the main stag-mwc conda environment to their
  latest version.
- Some minor scripts affected by Pandas and Matplotlib updates were modified to
  work with the latest versions of those libraries.
- Updated groot to 1.1.2 that brings many performance improvements, but removed
  the built-in plotting functionality so the groot module no longer produces
  any plots. Removed size window filtering with BBMap from groot alignment rule,
  and renamed the groot config variable `index` to `index_dir` to better map to
  `--indexDir` used in groot CLI.

### Deprecated
- Older Slurm profiles for CTMR Gandalf and UPPMAX Rackham are now considered
  deprecated and will be removed in a future release.
- MetaWRAP support is considered deprecated and will eventually be replaced by
  another solution for metagenome assembly in a future release of StaG.

### Removed


## [0.5.1] 2022-12-06
### Added
- Produce Snakemake report in zip format instead of HTML due to the HTML report being
  broken in the later versions of Snakemake.
- Add KrakenUniq as taxonomic profiler as an alternative with lower false
  positive rate than Kraken2.
- Added samplesheet as alternative input file selection method, this also
  enables providing custom sample names that are not based on pattern in input
  filenames.
- Samplesheet can be used to specify remote input files from S3 or HTTP/HTTPS sources.
- Added `run_krona` setting for taxonomic profilers to make it possible to disable Krona
  table and plot creation.

### Fixed
- Corrected typo in `host_removal` rule concerning `keep_kreport` config flag.
- Corrected typo in bowtie2 annotation counts output files leading to workflow
  complaining about missing output files.
- Removed unintended stdout printouts from various helper scripts and some
  MetaPhlAn related rules.
- Removed outdated mentions of MetaPhlAn2 in report.

### Changed
- Replaced CircleCI automatic testing workflow with one implemented with Github actions.
- Updated MetaPhlAn to version 4.0.3.
- Updated HUMAnN to version 3.6.
- Modified area and MetaPhlAn heatmap plotting scripts to better deal
  with MetaPhlAn 4 output formats.
- Updated the documentation to reflect recent changes in StaG.
- Updated KrakenTools to v1.2
- Updated `scripts/join_tables.py` to v1.1, which includes support for skipping
  lines before the header.
- Improved automatic report generation code in main Snakefile to be more
  robust. Now works well also when --use-singularity or --jobs are used
  simultaneously with --report.

### Removed
- Removed old unmaintained DB download rules for groot, kaiju, kraken2. 


## [0.5.0] 2021-11-18
### Added
- Biobakery update: updated MetaPhlAn and HUMAnN to version 3 as well as introducing 
  StrainPhlAn3 for strain-level genomics.
- Added `$TMPDIR` variable which can be specified in config.yaml. It is normally not required
  to specify `$TMPDIR` but might be necessary to run HUMAnN due to large intermediary files.
- New internal StaG feature to better handle user messages and defer them for printing after
  the workflow finished execution so they don't get lost in the verbose log printout from
  Snakemake.
- Added new feature to automatically remove intermediaries with the use of
  `keep_` flags in the config file. Currently available for Quality control
  (fastp), host removal (kraken2), taxonomic profiling with Kraken2 and
  MetaPhlAn3.
- Automatic builds of Singularity images using Github actions. Images get uploaded to
  the Github Package repository at `ghcr.io/ctmrbio/stag-mwc:<tag>-<branch>`.

### Fixed
- Updated pandas to 1.2.1 to fix issue with `preprocessing_summary.py` failing.
- `preprocessing_summary.py` now parses the correct number of unclassified reads.
- Rule `bracken_kreport` now correctly produces kreport output file for
  downstream processing.
- Fixed typo in count summary output filenames for BBMap.
- Increased time allocations for host removal in `ctmr_gandalf` cluster config.
- Limited job allocations to one node in `ctmr_gandalf` cluster config.
- Fixed bug where unspecified kraken2 database did not raise expected WorkflowError.

### Changed
- Updated Kraken2 to 2.1.2 and added `--minimum-hit-groups` setting in config file.
- Updated fastp to 0.23.0.
- Updated Kaiju to 1.8.2.
- Updated BBMap to 38.93.
- Updated MultiQC to 1.11.

### Removed


## [0.4.1] 2021-02-02
### Added
- Created Singularity images for all conda environments. Run with
  `--use-singularity` (do not combine with `--use-conda`).
- New cluster profile "pseudo-rules" for anonymous rules for mappers: `bbmap`
  and `bowtie2` can now accept threads from `n` in the cluster profile. They
  still use the time allocation for the `__default__` rule, however.
- Added possibility to use `extra:` to define additional arguments passed on to
  Slurm submissions. Useful to request e.g. fat nodes with `extra: "-C fat"` 
- Added custom reimplementation of AMRPlusPlus v2.0 which can be executed with 
  either `--use-singularity` or `--use-conda`.

### Fixed
- The host removal module now correctly identifies setting `host_removal: False`
  in the config file. Thank you chrsb!

### Changed
- Do not combine `--use-singularity` with `--use-conda` anymore. The new 
  Singularity images already contain all dependencies.
- All rules now define the number of threads from cluster_config if defined.
  Old defaults are still used for local execution.
- The shebang of `area_plot.py` has been changed to work in more environments.
- Implemented workaround for error caused by automatic report generation when
  using Singularity.
- Disabled taxonomic area plot for Kaiju outputs due to issues processing the
  output files.

### Removed


## [0.4.0] 2020-02-18
### Added
- Added resource limiter for HUMAnN2 due to its intense use of huge temporary
  files in the output folder. Activated with --resources humann2=X, where X is
  the max number of parallel instances of humann2 to run.
- Added groot report parameters `covcutoff` and `lowcov` to config file.
- Added automatic plot of proportion human reads. Included in report.
- Added assembly and binning using MetaWRAP.
- Added the possibility to run in Singularity with conda using 
  `--use-singularity --use-conda`.
- Added more MetaPhlAn2 data in output report.
- Added HUMAnN2 summary tables to output report.
- Added combined table and Krona plot for Kraken2 to output report.
- Added metagenomic assembly, binning and "blobology", using MEGAHIT or SPAdes,
  with binning using CONCOCT and MetaBat (MaxBin2 is not working), implemented
  via MetaWrap.
- Added KrakenTools from Jennifer Lu under the MIT license.
- Added basic syntax and DAG validation test in CircleCI
- Added possibility to skip read QC and/or host removal: will symlink relevant files
  files into the relevant output directories so Snakemake can continue
  without performing read QC and/or host removal.
- Added MultiQC, mainly to summarize fastp logs.
- Added Bracken abundance estimation on Kraken2 report files; added Bracken to
  the StaG conda environment. Also added Bracken abundance filtering rules so
  users can include/exclude certain taxa.
- Added "all_samples" summary output files in a more common format for all
  taxonomic profilers, called `mpa_style`. They are not identical, but very
  similar, with full lineage listings for all detected taxa.
- Added pigz v2.4 to the main StaG conda environment.
- Added summary with read counts passing preprocessing steps as a table and
  basic line plot. Only runs if both read QC and host removal are performed.

### Fixed
- Fixed bug in Slurm profile handling of cancelled/failed jobs.
- MetaPhlAn2 rule now correctly detect if no database path has been entered in
  the config file.
- HUMAnN2 rule now correctly detects if no database path has been entered in
  the config file.
- Kraken2 rule now correctly detects if no database is available at the path
  given in the config file.

### Changed
- Updated Python to 3.7 in main stag-mwc conda environment.
- Updated Kaiju to 1.7.2.
- Updated BBMap to 38.68.
- Updated sambamba to 0.7.0.
- Updated Kraken2 to 2.0.8_beta.
- Updated seaborn to 0.8.1, added fastcluster to main stag-mwc conda env, installed via pip.
- Updated MetaPhlAn2 to 2.96.1.
- Updated HUMAnN2 to 2.8.1.
- Updated GROOT to v0.8.5.
- Updated plot_metaphlan2_heatmap.py to 0.3.
- Changed output filenames ending with `.tsv` to `.txt` to avoid pretty HTML
  representations in report.
- Replaced BBMap-based host removal with Kraken2, substantially reducing time
  and resources requirements. 
- Added read length window filter before groot alignment step.
- Change logdir of remove_human rule to LOGDIR/remove_human instead of
  OUTDIR/logs/remove_human.
- Improved make_count_table.py so it can use TSV annotation files with multiple
  columns. Added config setting for which columns to include.
- Cleaned up sketch comparison cluster heatmap plotting script, making it more 
  robust to variations in output from different BBTools versions.
- Changed the call to merge_metaphlan_tables.py due to undocumented
  CLI change in latest conda version.
- Changed Kaiju summary report output filenames.
- Split biobakery environment into metaphlan2 and humann2 so users only
  interested in MetaPhlAn2 do not have to download the huge HUMAnN2.
- Replaced the outdated metaphlan_hclust_heatmap.py with a custom
  plot_metaphlan2_heatmap.py script.
- Defined some low-impact summary and plotting rules as localrules.
- Reworded all rules relating to human removal to "host removal", and changed output folder
  structure accordingly.
- Renamed output folder and file names for quality and adapter trimming.
- Set Kraken2 --confidence to 0.1 by default.
- Adjusted HUMAnN2 cores to 20 (up from 8).
- Adjusted MetaPhlAn2 cores to 5 (up from 4).


### Removed
- Removed outdated database download rules for Centrifuge, MetaPhlAn2, Kaiju,
  Kraken2, HUMAnN2.
- Replaced FastQC + BBDuk with fastp adapter trimming and quality filtering.
- Removed support for Centrifuge


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

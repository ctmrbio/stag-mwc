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


## [0.5.1] Unreleased
### Added

### Fixed
- Corrected typo in `host_removal` rule concerning `keep_kreport` config flag.
- Corrected typo in bowtie2 annotation counts output files leading to workflow
  complaining about missing output files.
- Fixed CheckM installation in assembly singularity file so assembled bin consolidation
  from MetaWRAP now works as intended.

### Changed
- Replaced CircleCI automatic testing workflow with one implemented with Github actions.
- Updated MetaWRAP to 1.3.2.

### Removed


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

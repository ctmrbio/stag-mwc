# Instructions for configuration of the StaG-mwc Snakemake workflow
This file contains a condensed version of the official online documentation
available at https://stag-mwc.readthedocs.io.

## Specify input files
There are two ways to specify input files for StaG:

  1. Point StaG to a folder containing paired-end FASTQ files with a structured filename pattern.
  2. Prepare a sample sheet with sample names and paths/URLs to paired-end FASTQ input files

### Option 1: Input folder
If you have all input files (or symlinks to input files) located in a single
folder, and your input files have a structured filename containing a unique
sample identifier, this method of picking input files is the most convenient.
Open `config.yaml` in your favorite editor and change input file settings
under the `Run configuration` heading: 

 - the input directory
 - the input filename pattern

They can be declared using absolute or relative filenames. 

### Option 2: Sample sheet
If your input FASTQ files are spread across several filesystem locations or
potentially exist in remote locations (e.g. S3), or your input FASTQ filenames
do not follow a common filename pattern, the samplesheet option is the most
convenient. The samplesheet input option also allows you to specify custom
sample names that are not derived from a substring of the input filenames.

The format of the samplesheet is tab-separated text and it must contain a
header line with at least the following three columns: `sample_id`,
`fastq_1`, and `fastq_2`. An example file could look like this (columns are
separated by TAB characters):

    sample_id  fastq_1                             fastq_2
    ABC123     /path/to/sample1_1.fq.gz            /path/to/sample1_2.fq.gz
    DEF456     s3://bucketname/sample_R1.fq.gz     s3://bucketname/sample_R2.fq.gz
    GHI789     http://domain.com/sample_R1.fq.gz   http://domain.com/sample_R2.fq.gz

Open `config.yaml` in your favorite editor and enter the path to a
samplesheet TSV file that you have prepared in advance in the `samplesheet`
field under the `Run configuration` heading. Input files can be located
anywhere, i.e. their locations are not restricted to the repository folder and
they can even be located in remote storage systems like S3 or a public HTTP
URL.

If the samplesheet setting is configured in the `config/config.yaml` file, or
provided on the command line when running by utilizing Snakemake’s built-in
functionality for modifying configuration settings via the command line
directive `--config samplesheet=path/to/samplesheet.tsv`, it will override any
input folder settings configured in the `config/config.yaml` file. 



## Select which tools to run
Next, configure the settings under the `Pipeline steps included` heading. This
is where you define what steps should be included in your workflow. Simply
assign `True` or `False` to the steps you want to include. The default
configuration file already sets both `qc_reads` and `host_removal` to `True`:
these two steps are the primary read processing steps and most other steps
depends on host filtered reads (i.e. the output of the `host_removal` step).
Note that these two steps will pretty much always run, regardless of their
setting in the config file, because they produce output files that almost all
other workflow steps depend on. 

## Fill in required settings for the selected tools
Further down in `config.yaml` are sections for each individual tool, and most
tools have some required settings that need to be configured, typically paths
to reference databases. These are marked with `[Required]` and there are
comments explaining what is expected for each setting.

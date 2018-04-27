# Contributing guidelines for StaG-mwc

## Issue tracker
We use the issue tracker in Github. Submit issues for things such as
bug reports, feature requests, or general improvement discussion topics.

## Submitting changes
The typical procedure to develop new features or fix bugs in StaG-mwc looks
something like this:

1. Fork or clone the repository.
2. Create a branch with a descriptive name based on your intended changes using
   dashes to separate words, e.g. `branch-to-add-megahit-assembly-step`.
3. Insert your code into the respective folders, i.e. scripts, rules and envs.
   Define the entry point of the workflow in the Snakefile and the main
   configuration in the config.yaml file.
4. Commit changes to your fork/clone.
5. Create a pull request (PR) with some motivation behind the work you have
   done and possibly some explanations for tricky bits.


## Code organization
The workflow code is organized in the following folder structure inside the
git repository:

    cluster_configs     # Dir containing Snakemake profiles for clusters
	docs                # Sphinx documentation, only sources!
	envs                # Conda environment files
	rules               # The actual Snakemake rules that make up the workflow
	scripts             # Scripts used in workflow rules
	utils               # Auxiliary scripts, not used in workflow
	config.yaml         # The default config file
	CONTRIBUTING.md     # This document
	LICENSE             # LICENSE
	README.md           # The README shown in the github repo
	Snakefile           # The main workflow script

### cluster_configs
The `cluster_configs` directory should contain either:

1. Folders representing entire [Snakemake cluster profiles](https://snakemake.readthedocs.io/en/stable/executable.html#profiles).
2. Single `yaml` or `json` [cluster config files](http://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html?highlight=cluster-config#cluster-configuration).

### docs 
The documentation for the project is built automatically by
[readthedocs](www.readthedocs.org) upon every commit. The HTML documentation is
available at https://stag-mwc.readthedocs.org. Feel free to improve the
documentation, but avoid committing anything but source documents to the repo.
The documentation is written using Sphix, so all documentation sources are
written in [reStructuredText](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

### envs
The `envs` folder contains conda environments for the workflow. The ambition is
that all dependencies should be included in the main `stag-mwc.yaml`
environment, unless they have incompatible dependencies, to reduce the amount
of conda environments in total. It is absolutely preferable if all tools used
in the workflow are available via conda (either default channels, or bioconda,
conda-forge, etc.).

### rules
All workflow rules are organized in the `rules` folder. It contains a directory
hierarchy organized by overall function in the workflow, e.g., the subfolder
`taxonomic_profiling` contains rules for all taxonomic profiling tools. It is
recommended to keep one file per logical unit or tool. 

The overall concept of StaG-mwc is that analyses are performed on trimmed/cleaned
reads that have had human sequences removed, so rules should generally start
with the clean FASTQ files output from the `remove_human` step. This is of course
only a general recommendation, and some tools require the raw reads for their
analysis.

Each rule file should define the expected output files of that module and add
them to the `all_outputs` object, defined in the main Snakefile. This is
designed to allow some inclusion logic in the main Snakefile, so components can
be turned on or off without too much trouble. Output should typically be in a
subfolder inside the overall `outdir` folder. `outdir` is available as a string
in all rule files, as it is defined in the main Snakefile based on the value
set in `config.yaml`. 

Tools that require databases or other reference material to work can be
confusing or annyoing to users of the workflow. To minimize the amount of
effort required to use StaG-mwc with such tools, it is recommended to provide a
rule to download and prepare the required databases. See examples in e.g.
`rules/preproc/remove_human.smk` and `rules/taxonomic_profiling/metaphlan2.smk`. 
Make sure to include some code to verify the existence of required reference
files at the top of the rule file, so users get a fairly clean `WorkflowError`
output if a required input file is missing. See examples of that in e.g. 
`rules/mappers/bowtie2.smk`. 

### scripts
The `scripts` folder contains all scripts required by workflow rules. These
are typically read summarization or plotting scripts, but anything that is
used by rules that aren't specifically rules themselves should go in here.

### utils
The `utils` folder contains auxiliary scripts or tools that are useful in the
context of StaG-mwc, but are not necessarily used directly by the workflow.

### config.yaml
The configuration file is the main point of configuration of StaG-mwc. It
should include reasonable default values for all important settings for the
implemented tools, and contain brief comments for most values. The
configuration file is divided into separate sections, reflecting the overall
structure of the workflow.

The top section is `Pipeline steps included`, which should contain Boolean
values for all steps to include in the workflow. These variables are checked in
the main Snakefile and used to determine if rules should be included in the
workflow or not. 

The following sections reflect the folder structure inside the `rules` folder,
and are organized by tool name. If the same tool is used in several steps, it
is recommended to choose a more descriptive name. 

### Snakefile
`Snakefile` is the main workflow script. This is where all the different rules
defined in the `rules` folder are included into the overall Snakemake workflow. 
The Snakefile defines the set of samples to analyze, by searching for files
matching the `input_fn_pattern` in the `inputdir`, both defined in `config.yaml`. 

The Snakefile defines the following variables that are accessible in all rule
files:

* `config` - Python dict containing all configuration parameters defined in `config.yaml`.
* `outdir` - Python string with the path to the output directory
* `all_outputs` - Python list that contains strings with the path to all expected output files
* `SAMPLES` - Python set containing all sample names identified from the input file names.


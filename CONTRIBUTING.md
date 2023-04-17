# Contributing guidelines for StaG-mwc

## Issue tracker
We use the issue tracker in Github. Submit issues for things such as
bug reports, feature requests, or general improvement discussion topics.

## Submitting changes
The main branch of StaG-mwc should always be stable and reliable. All
development should be based on the develop branch. Please create new feature
branches from the develop branch. The develop branch is then merged into the
master branch when enough improvements have accumulated. The typical procedure
to develop new features or fix bugs in StaG-mwc looks something like this:

1. Fork or clone the repository.
2. Checkout the develop branch and create a new feature branch from there.
   Use a descriptive name and use dashes to separate words:
   ```
   git checkout develop
   git checkout -b add-megahit-assembly-step
   ```
3. Write or modify code in the scripts, rules and envs folders. Define the
   entry point of the workflow in the Snakefile and the main configuration in the
   config.yaml file.
4. If a new feature has been added, document it in the Sphinx documentation.
4. Commit changes to your fork/clone.
5. Create a pull request (PR) to the develop branch  with some descriptions of
   the work you have done and possibly some explanations for potentially tricky
   bits.
6. When the feature is considered complete, we bump the version number depending
   on the size and impact of the PR before merging the PR to the develop branch.


### Releases
New releases are made whenever enough new features have accumulated on the
develop branch. Before creating a new release, create a staging branch off of
the develop branch, and ensure the following things have been taken care of:

* All pending features that should be included in the upcoming release are
  merged.
* Double check that documentation is available and up-to-date for implemented
  features.
* Check that the version number in the documentation matches the Snakefile.

Then, merge the staging branch into master, squashing all commits, and tag
the new release. Afterwards, merge the staging branch back to develop so all
changes in the staging branch are present in develop.


## Code organization
The repo aims to follow the established Snakemake "best practices" for [code
organization](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility).


### Profiles
The `profiles` directory should contain folders representing [Snakemake
profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
suitable for use with commonly encountered cluster scheduling systems (e.g.
Slurm).


### Docs 
The documentation for the project is built automatically by
[readthedocs](www.readthedocs.org) upon every commit. The HTML documentation is
available at https://stag-mwc.readthedocs.org. Feel free to improve the
documentation, but avoid committing anything but source documents to the repo.
The documentation is written using Sphinx, so all documentation sources are
written in [reStructuredText](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).


### Envs
The `envs` folder contains conda environments for the workflow. The ambition is
that all dependencies should be included in the main `stag-mwc.yaml`
environment, unless they have incompatible dependencies, to reduce the amount
of conda environments in total. It is absolutely preferable if all tools used
in the workflow are available via conda (either default channels, or bioconda,
conda-forge, etc.), and Singularity. In cases where no official container images
are available, we build custom Singularity images using Singularity definition
files stored in the same folder. In most cases the same conda environment
specification can simply be installed in the container. Containers are
automatically built using a Github Action.


### Rules
All workflow rules are organized in the `rules` folder. It contains a directory
hierarchy organized by overall function in the workflow, e.g., the subfolder
`taxonomic_profiling` contains rules for all taxonomic profiling tools. It is
recommended to keep one file per logical unit or tool, so they can be easily
added by a single ``include:`` in the main Snakefile.

The overall concept of StaG-mwc is that analyses are performed on trimmed/cleaned
reads that have had human sequences removed, so rules should generally start
with the clean FASTQ files output from the `remove_human` step. This is of
course only a general recommendation, and some tools naturally require the raw
reads for their analysis.

Each rule file should define the expected output files of that module and
conditionally add them to the `all_outputs` object defined in the main
Snakefile. Wrap adding of the files to the ``all_outputs`` list in an
if-statement conditioned on the booleans defined in ``config.yaml`` under the
``Pipeline steps included`` section. This is the preferred way, as it makes
Snakemake aware of all rules, and uses its own dependency resolution engine to
figure out the rule graph to produce the desired output files. This way, users
can easily change which output files they want in ``config.yaml`` in an easy
way, and Snakemake figures out the rest.  Output should typically be in a
subfolder inside the overall `outdir` folder. `outdir` is available as a string
in all rule files, as it is defined in the main Snakefile based on the value
set in `config.yaml`.

Declare paths to input, output and log files using the pathlib Path objects
INPUTDIR, OUTDIR, and LOGDIR. Note that Snakemake is not yet fully pathlib
compatible so convert Path objects to strings inside `expand` statements and
log file declarations. In future versions of Snakemake this will not be necessary.

Tools that require databases or other reference material to work can be
confusing or annyoing to users of the workflow. To minimize the amount of
effort required to use StaG-mwc with such tools, it is recommended to provide a
rule to download and prepare the required databases. See examples in e.g.
`rules/preproc/remove_human.smk` and `rules/taxonomic_profiling/metaphlan2.smk`. 
Make sure to include some code to verify the existence of required reference
files at the top of the rule file, so users get a fairly clean `WorkflowError`
output if a required input file is missing. See examples of that in e.g. 
`rules/mappers/bowtie2.smk`. 

Other types of messages that users might need to see should be added to the
list of deferred user messages by calling either
`user_messages.info("message")` or `user_messages.warn("message")`. Messages
added this way will be printed upon pipeline failure or completion.

### Scripts
The `scripts` folder contains all scripts required by workflow rules. These
are typically read summarization or plotting scripts, but anything that is
used by rules that aren't specifically rules themselves should go in here.


### Utils
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

* `config` - Python `dict` containing all configuration parameters defined in
  `config.yaml`.
* `citations` - Python `set` of publication citation strings loaded from
  `rules/citatons.py`. 
* `INPUTDIR` - Python `pathlib.Path` object for the input directory
* `OUTDIR` - Python `pathlib.Path` object for the output directory
* `LOGDIR` - Python `pathlib.Path` object for the log directory
* `TMPDIR` - Python `pathlib.Path` object for the temporary directory directory specified
  in `config.yaml`.
* `DBDIR` - Python `pathlib.Path` object for the database directory specified in
  `config.yaml`.
* `user_messages` - A `UserMessages` instance with methods `info` and `warn` to
  show messages to the user.
* `all_outputs` - Python `list` that contains strings with the path to all
  expected output files
* `SAMPLES` - Python `set` containing all sample names identified from the input
  file names.


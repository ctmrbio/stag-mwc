# StaG Metagenomic Workflow Collaboration (mwc)

[![DOI](https://zenodo.org/badge/125840716.svg)](https://zenodo.org/badge/latestdoi/125840716)
[![Snakemake](https://img.shields.io/badge/snakemake-≥4.8.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![CircleCI](https://circleci.com/gh/ctmrbio/stag-mwc.svg?style=svg)](https://circleci.com/gh/ctmrbio/stag-mwc)

![StaG mwc logo](docs/source/img/stag_head_text.png "StaG mwc")

This repo contains the code for a Snakemake workflow of the StaG Metagenomic
Workflow Collaboration (mwc). Currently, the project focus is a barebones
metagenomics analysis workflow to produce primary output files from several
different metagenomic analysis tools. 

Go to https://stag-mwc.readthedocs.org for the full documentation.

## Usage

### Step 0: Install conda and Snakemake
[Conda](https://conda.io/docs/) and
[Snakemake](https://snakemake.readthedocs.io) are required to be able to use
StaG-mwc. Most people would probably want to install
[Miniconda](https://conda.io/miniconda.html) and install Snakemake into their
base environment. Conda will automatically install the required versions of 
all tools required to run StaG-mwc.

### Step 1: Clone workflow
To use StaG-mwc, you need a local copy of the workflow repository. Start by
making a clone of the repository: 

    git clone git@github.com:ctmrbio/stag-mwc

If you use StaG-mwc in a publication, please credit the authors by citing
either the URL of this repository or the project's DOI. Also, don't forget to
cite the publications of the other tools used in your workflow.

### Step 2: Configure workflow
Configure the workflow according to your needs by editing the file
`config.yaml`. The most common changes include setting the paths to input and
output folders, and configuring what steps of the workflow should be included
when running the workflow.

### Step 3: Execute workflow
Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores N

This will run the workflow locally using `N` cores. It is also possible to run
it in a Slurm-managed cluster environment, e.g. on UPPMAX Rackham:

    snakemake --use-conda --profile cluster_configs/rackham

Make sure you edit the Slurm project account in
`cluster_configs/rackham/rackham.yaml`. Refer to the official [Snakemake
documentation](https://snakemake.readthedocs.io) for further details on how to
run Snakemake workflows on other types of cluster resources.


## Testing
A very basic continuous integration test is currently in place. It merely
validates the syntax by trying to let Snakemake build the dependency graph if
all outputs are activated.


## Contributing
Refer to the contributing guidelines in `CONTRIBUTING.md` for instructions on
how to contribute to StaG-mwc.

If you intend to modify or further develop this workflow, you are welcome to
fork this reposity. Please consider sharing potential improvements via a pull
request.

## Citing
If you find StaG-mwc useful in your research, please cite the Zenodo DOI:
https://zenodo.org/record/1483891


# Logo attribution
<a href="https://www.freepik.com/free-photos-vectors/animal">Animal vector created by Patrickss - Freepik.com</a>

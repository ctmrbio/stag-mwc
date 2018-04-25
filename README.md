# StaG Metagenomic Workflow Collaboration (mwc)

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/mwc.svg?branch=master)](https://travis-ci.org/snakemake-workflows/mwc)
![StaG mwc logo](docs/source/img/stag_head_text.png "StaG mwc")

This repo contains the code for a Snakemake workflow of the StaG Metagenomic
Workflow Collaboration (mwc). Currently, the project focus is a barebones
metagenomics analysis workflow to produce primary output files from several
different metagenomic analysis tools. 

## Authors

* Fredrik Boulund (@boulund)
* Lisa Olsson (@lis4matilda)
* (your name here)

## Usage

### Step 0: Install conda and Snakemake
[Conda](https://conda.io/docs/) and
[Snakemake](https://snakemake.readthedocs.io) are required to be able to use
StaG-mwc. Most people would probably want to install
[Miniconda](https://conda.io/miniconda.html) and install Snakemake into their
base environment. Conda will automatically install the required versions of 
all tools required to run StaG-mwc.

### Step 1: Install workflow
<!--download and extract the [latest release](https://github.com/snakemake-workflows/mwc/releases). -->

If you simply want to use this workflow, clone the repository: `git clone
git@github.com:boulund/mwc`. If you intend to modify or further develop this
workflow, you are welcome to fork this reposity. Please consider sharing
potential improvements via a pull request.

If you use StaG-mwc in a publication, please credit the authors by citing
the URL of this repository and, when available, its DOI. Also, don't forget to
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
it in a Slurm-managed cluster environment using e.g. 

    snakemake --use-conda --jobs 999 --cluster-config cluster_configs/rackham.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time}"

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further
details on how to run Snakemake workflows on other types of cluster resources.


### Automatic database download
The workflow offers steps that can automatically download the required
reference databases. Note that this step is normally only required once, as
previously downloaded databases are reused. See the 
[official documentation](https://stag-mwc.readthedocs.org) for more information.


## Testing
Tests are currently not implemented. The ambition is that mwc will contain
extensive tests to verify functionality. They should be executed via continuous
integration with Travis CI. 


## Contributing
Refer to the contributing guidelines in `CONTRIBUTING.md` for instructions on how to
contribute to StaG-mwc.

# Logo attribution
<a href="https://www.freepik.com/free-photos-vectors/animal">Animal vector created by Patrickss - Freepik.com</a>

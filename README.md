# StaG Metagenomic Workflow Collaboration (mwc)

[![DOI](https://zenodo.org/badge/125840716.svg)](https://zenodo.org/badge/latestdoi/125840716)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.24.0-brightgreen.svg)](https://snakemake.bitbucket.io)

![StaG mwc logo](docs/source/img/stag_head_text.png "StaG mwc")

The StaG Metagenomic Workflow Collaboration (mwc) project focuses on providing
a metagenomics analysis workflow suitable for microbiome research and general
metagenomics analyses. 

Please visit https://stag-mwc.readthedocs.org for the full documentation.


## Usage

### Step 0: Install conda and Snakemake
[Conda](https://conda.io/docs/) and
[Snakemake](https://snakemake.readthedocs.io) are required to be able to use
StaG-mwc. Most people would probably want to install
[Miniconda](https://conda.io/miniconda.html) and install Snakemake into their
base environment. When running StaG with the `--use-conda` or
`--use-singularity` flags, all dependencies are managed automatically. If
using conda it will automatically install the required versions of all tools
required to run StaG-mwc. There is no need to combine the conda and singularity
flags: the Singularity images used by the workflow already contain all required
dependencies.

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

Note that in all examples above, `--use-conda` can be replaced with
`--use-singularity` to run in Singularity containers instead of using a locally
installed conda. Read more about it under the Running section in the docs.


## Testing
A very basic continuous integration test is currently in place. It merely
validates the syntax by trying to let Snakemake build the dependency graph if
all outputs are activated. Suggestions for how to improve the automated
testing of StaG-mwc are very welcome!


## Contributing
Refer to the contributing guidelines in `CONTRIBUTING.md` for instructions on
how to contribute to StaG-mwc.

If you intend to modify or further develop this workflow, you are welcome to
fork this reposity. Please consider sharing potential improvements via a pull
request.


## Citing
If you find StaG-mwc useful in your research, please cite the Zenodo DOI:
https://zenodo.org/badge/latestdoi/125840716

# Logo attribution
<a href="https://www.freepik.com/free-photos-vectors/animal">Animal vector created by Patrickss - Freepik.com</a>

Installing |full_name|
=======================
|full_name| depends on a number of individual tools, most prominently it
requires `Snakemake`_, which is used to manage the overall analysis workflow.
In addition, another important component for the |full_name| framework is
`Conda`_, which is used to manage additional dependencies. This documentation
assumes you already have `git`_ installed. 

.. _Snakemake: https://snakemake.readthedocs.io
.. _Conda: https://conda.io/docs/
.. _miniconda3: https://conda.io/miniconda.html
.. _git: https://git-scm.com/


Install conda and Snakemake
***************************
The first two things you need to install are:

1. `Conda`_
2. `Snakemake`_

The recommended way to get started using |full_name| is to download and install
`Conda`_. A good starting point is a clean `miniconda3`_ installation.
Miniconda3 is quick to install and does not require administrator permissions.

After installing `Conda`_ and activating the base environment, install
snakemake into your base environment::

    (base)$ conda install snakemake

These are the only external dependencies you need to install manually. The
correct versions of any remaining dependencies will be automatically downloaded
and installed when you run the workflow the first time.


Download the workflow code
**************************
Clone the |full_name| repository using git to get a copy of the workflow::

    (base)$ git clone https://www.github.com/boulund/stag-mwc

This will clone the repo into a folder called ``stag-mwc`` inside your current
working directory. You will manage all workflow-related business from inside this
folder (i.e. configuring and running the workflow).


Congratulations
***************
You have now installed |full_name|. 

The intended way of working with |full_name|, or really any Snakemake
workflow, is that you download/clone a complete copy of the repository for
each analysis you intend to make. That way, the local copy you have will
remain after the analysis has been run, so you can go back and see exactly
what was run, and how. This forms the basis of how Snakemake enables
traceability and reproducibility.
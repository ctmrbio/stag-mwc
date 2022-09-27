Running |full_name|
===================
You need to configure a workflow before you can run |full_name|. The code 
you downloaded in the previous ``git clone`` step includes a file called 
``config.yaml``, which is used to configure the workflow. 

Selecting input files
*********************
There are two ways to define which files |full_name| should run on: either
by specifying an input directory and a filename pattern, or by providing
a sample sheet. The two ways are exclusive and cannot be combined, so you have
to pick the one that suits you best. 

Input directory
---------------
If your input FASTQ files are all in the same folder and they all follow the
same filename pattern, the input directory option is often the most convenient.

Open ``config.yaml`` in your favorite editor and change input file settings
under the ``Run configuration`` heading: the input directory, the input
filename pattern. They can be declared using absolute or relative filenames
(relative to the |full_name| repository directory). Input and output
directories can technically be located anywhere, i.e. their locations are not
restricted to the repository folder, but it is recommended to keep them in the
repository directory. A common practice is to put symlinks to the files you
want to analyze in a folder called ``input`` in the repository folder.

Samplesheet
-----------
If your input FASTQ files are spread across several filesystem locations or
potentially exist in remote locations (e.g. S3), or your input FASTQ filenames
do not follow a common filename pattern, the samplesheet option is the most
convenient. The samplesheet input option also allows you to specify custom
sample names that are not derived from a substring of the input filenames.

The format of the samplesheet is tab-separated text and it must contain a
header line with at least the following three columns: ``sample_id``,
``fastq_1``, and ``fastq_2``. An example file could look like this (columns are
separated by TAB characters)::

   sample_id  fastq_1                           fastq_2
   ABC123     /path/to/sample1_1.fq.gz          /path/to/sample1_2.fq.gz
   DEF456     s3://bucketname/sample_R1.fq.gz   s3://bucketname/sample_R2.fq.gz

Open ``config.yaml`` in your favorite editor and enter the path to a
samplesheet TSV file that you have prepared in advance in the ``samplesheet``
filed under the ``Run configuration`` heading: They paths can be declared using
absolute or relative filenames (relative to the |full_name| repository
directory). Input files can be located anywhere, i.e. their locations are not
restricted to the repository folder and they can even be located in remote
storage systems like S3.


Configuring which tools to run
******************************
Next, configure the settings under the ``Pipeline steps included`` heading.
This is where you define what steps should be included in your workflow. Simply
assign ``True`` or ``False`` to the steps you want to include. Note that the
default configuration file already includes ``qc_reads`` and ``host_removal``.
These two steps are the primary read processing steps and most other steps
depends on host filtered reads (i.e. the output of the ``host_removal`` step).
Note that these two steps will pretty much always run, regardless of their
setting in the config file, because they produce output files that almost all
other workflow steps depend on. 

.. note:: 

    You can create several copies of ``config.yaml``, named whatever you want,
    in order to manage several analyses from the same |full_name| directory.
    If you create a copy called e.g. ``microbime_analysis.yaml``, you can easily
    run the workflow with this configuration file by using the ``--configfile``
    commandline argument when running the workflow.

A reference database is required in order to run the ``host_removal`` step. If
you already have it downloaded somewhere, point |full_name| to the location
using the ``db_path`` parameter under the ``remove_host`` section of ``config.yaml``.

The config file contains a parameter called ``email``. This can be used to have
the workflow send an email after a successful or failed run. Note that this 
requires that the Linux system your workflow is running on has a working email
configuration. It is also quite common that most email clients will mark email sent
from unknown random computers as spam, so don't forget to check your spam folder.


Running
*******
It is recommended to run Snakemake with the ``-n``/``--dryrun`` argument before
starting an analysis for real. Executing a dryrun will let Snakemake check that
all the requirements are available and it will then print a summary of what it
intends to do, without actually doing anything. After finishing the
configuration by editing ``config.yaml``, test your configuration with::

    snakemake --dryrun

If you are satisfied with the workflow plan output by the dryrun, you can run
the workflow. The typical command to run |full_name| on your local computer
is::

    snakemake --use-conda --cores N

where ``N`` is the maximum number of cores you want to allow for the
workflow. Snakemake will automatically reduce the number of cores available
to individual steps to this limit. Another variant of ``--cores`` is called
``--jobs``, which you might encounter occassionally. The two commands are
equivalent.

.. note::

    If several people are running StaG-mwc on a shared server or on a shared
    file system, it can be useful to use the
    ``--singularity-prefix``/``--conda-prefix`` parameter to use a common
    folder to store the conda environments created by StaG-mwc, so they can be
    re-used between different people or analyses. This reduces the risk of
    producing several copies of the same conda environment in different
    folders. This is also often necessary when running on cluster systems where
    paths are usually very deep. Then for example create a folder in your home
    directory and use that with the
    ``--singularity-prefix``/``--conda-prefix`` option.

If you want to keep your customized ``config.yaml`` in a separate file, let's 
say ``my_config.yaml``, then you can run snakemake using that custom configuration 
file with the ``--configfile my_config.yaml`` command line argument.

Another useful command line argument to snakemake is ``--keep-going``. This will 
instruct snakemake to keep going even if a job should fail, e.g. maybe the
taxonomic profiling step will fail for a sample if the sample contains no assignable
reads after quality filtering (extreme example).

If you are having trouble running |full_name| with conda, try with Singularity
(assuming you have Singularity installed on your system). There are pre-built
Singularity images that are ready to use with |full_name|. Consider using
``--singularity-prefix`` to specify a folder where Snakemake can download and
re-use the downloaded Singularity images for future invocations. The command to
run |full_name| with Singularity instead of conda is::

    snakemake --use-singularity --singularity-prefix /path/to/prefix/folder --dryrun

There are some additional details that need to be considered when using
Singularity instead of conda, most notably that you will have to specify bind
paths (specifying-bind-paths_) so that your reference databases are
accessible from within the containers when running |full_name|. It might look
something like this::

    snakemake --use-singularity --singularity-prefix /path/to/prefix/folder --singularity-args "-B /home/username/databases"

The above example assumes you have entered paths to your databases in
``config.yaml`` with a base path like the one shown in the above command
(e.g. ``/home/username/databases/kraken2/kraken2_human/``).


Running on cluster resources
****************************
In order to run |full_name| on a cluster, you need a special cluster
configuration file.  |full_name| ships with a pre-made configuration profile
for use on CTMR's Gandalf cluster and UPPMAX's Rackham cluster.  Find all
available cluster configuration profiles in the ``cluster_configs`` directory
in the repository. The cluster configuration profiles specify which cluster
scheduler account to use (e.g.  Slurm project account), as well as the number
of CPUs, time, and memory requirements for each individual step. Snakemake uses
this information when submitting jobs to the cluster scheduler.

When running on a cluster it will likely work best if you run StaG using
Singularity. The workflow comes preconfigured to download and use containers
from Singularity hub. To use Singularity launch Snakemake with the
``--use-singularity`` argument. 

.. _specifying-bind-paths: https://sylabs.io/guides/3.5/user-guide/bind_paths_and_mounts.html#specifying-bind-paths

.. note:: 

    Do not combine ``--use-conda`` with ``--use-singularity``.

    To prevent |full_name| from unnecessarily downloading the Singularity
    container images again between several projects you can use the
    ``--singularity-prefix`` to specify a directory where Snakemake can store
    the downloaded images for reuse between projects.

    Paths to databases need to be located so that they are accessible from
    inside the Singularity containers. It's easiest if they are all available
    from the same folder, so you can bind the main database folder into the
    Singularity container with e.g. ``--singularity-args "-B /db"``. Note that
    database paths need to specified in the config file so that the paths are
    correct from inside the Singularity container. Read more about specifying
    bind paths in the official Singularity docs: specifying-bind-paths_. 

To run |full_name| on e.g. CTMR's Gandalf, run the following command from
inside the workflow repository directory::

    snakemake --use-singularity --singularity-prefix /ceph/db/sing --singularity-args "-B /ceph" --profile cluster_configs/ctmr_gandalf

This will make Snakemake submit each workflow step as a separate cluster job
using the CPU and time requirements specified in ``ctmr_gandalf.yaml`` inside the
Rackham profile folder. The above command assumes you are using the default
``config.yaml`` configuration file. If you are using a custom configuration
file, just add ``--configfile <name_of_your_config_file>`` to the command line.

.. note::

    Make sure you edit ``cluster_configs/ctmr_gandalf/ctmr_gandalf.yaml`` to
    specify the Slurm project name to use for Slurm job submissions.

Some very lightweight rules will run on the submitting node (typically directly
on the login node), but the number of concurrent local jobs is limited to 1 in
the default profiles.


Execution report
****************
Snakemake provides facilites to produce an HTML report of the execution of the
workflow. An HTML report is automatically created when the workflow finishes.



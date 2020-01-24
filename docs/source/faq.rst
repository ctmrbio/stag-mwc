Frequently asked questions (FAQ)
================================
|full_name| is designed to be fairly flexible. This section answers common
questions on how to best utilize |full_name|'s flexibility. Please help expand
this section by making a Pull Request with suggestions.


Skip host removal
*****************
It is possible to skip host removal. This might be appropriate for example when
processing environmental samples that do not need host removal. It is
implemented in |full_name| in a special way: if the user sets ``host_removal:
False`` in the configuration file then |full_name| will put symlinks to the
``qc_reads`` output files in the ``host_removal`` output directory. This
"tricks" Snakemake into thinking that host removal actually occured, enabling
it to complete the dependency graph to process the data in downstream steps.


Using already trimmed and quality controlled sequencing data
************************************************************
In some cases your data has already been subjected to adapter and quality
trimming so you want to skip that step. In |full_name| it is possible to bypass
both the read QC step and the host removal steps if you need to. In order to do
so, you must "trick" Snakemake into thinking that those rules have already been
performed. The easiest way is simply to put symlinks to your input data files
inside the ``qc_reads`` output directory (or the ``host_removal`` output
directory), making sure that the fake "output" filenames match the hardcoded
output filenames from those steps: ``{sample}_{readpair}.fq.gz``.

For example, if your data has already been QC'd. Let's assume your input files
are located in a folder called ``input``::

   $ ls input
   sample1_1.fq.gz sample1_2.fq.gz

Now, you want to skip QC but not host removal. Then we need to create the output files
that Snakemake expects from the QC step::

   $ mkdir -pv output_dir/fastp
   $ cd output_dir/fastp
   $ ln -s ../../input/sample1_{1,2}.fq.gz .

Then, make sure to set the pipeline configuration so that QC is not performed:
``qc_reads: False``.

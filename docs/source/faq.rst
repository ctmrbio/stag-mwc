Frequently asked questions (FAQ)
================================
|full_name| is designed to be fairly flexible. This section answers common
questions on how to best utilize |full_name|'s flexibility. Please help expand
this section by making a Pull Request with suggestions.


Skip read QC
************
In some cases your data has already been subjected to adapter and quality
trimming so you want to skip that step. In |full_name| it is possible to bypass
both the read QC step and the host removal steps if you need to. In order to do
so, we must "trick" Snakemake into thinking that those rules have already been
performed. The rules for read QC and host removal are configured with bypasses
so that if the user sets ``host_removal: False`` or ``qc_reads: False``, those
steps will be bypassed by creating symlinks directly to the input files in the
respective output directories. 


Skip host removal
*****************
It is possible to skip host removal. This might be appropriate for example when
processing environmental samples that do not need host removal. It is
implemented in |full_name| in a special way: if the user sets ``host_removal:
False`` in the configuration file then |full_name| will put symlinks to the
``qc_reads`` output files in the ``host_removal`` output directory. This
"tricks" Snakemake into thinking that host removal actually occured, enabling
it to complete the dependency graph to process the data in downstream steps.

Pipeline stopped unexpectedly
*****************************
If pipeline ends with error or if the session is locked after being
unexpectedly disconnected and the pipeline needs to be restarted, you can try
to remove slurm metadadata files before restarting pipeline using::

    (base)$ rm -rfv .snakemake/metadata


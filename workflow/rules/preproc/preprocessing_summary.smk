# vim: syntax=python expandtab
# Summarize read counts passing through preprocessing.
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

localrules:
    preprocessing_summary

if config["qc_reads"] and config["host_removal"]["kraken2"]:
    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    read_counts = str(OUTDIR/"preprocessing_read_counts.txt")
    read_counts_plot = str(OUTDIR/"preprocessing_read_counts.pdf")

    all_outputs.append(read_counts)
    all_outputs.append(read_counts_plot)


rule preprocessing_summary:
    """Summarize read counts in preprocessing steps"""
    input:
        fastp=expand(str(LOGDIR/"fastp/{sample}.fastp.json"), sample=SAMPLES),
        kraken2=expand(str(LOGDIR/"host_removal/{sample}.kraken2.log"), sample=SAMPLES),
    output:
        table=report(str(OUTDIR/"preprocessing_read_counts.txt"),
            category="Preprocessing",
            caption="../../report/preprocessing_summary.rst"),
        plot=report(str(OUTDIR/"preprocessing_read_counts.pdf"),
            category="Preprocessing",
            caption="../../report/preprocessing_summary.rst"),
    log:
        stderr=str(LOGDIR/"preprocessing_summary.log"),
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    threads: 1
    shell:
        """
        workflow/scripts/preprocessing_summary.py \
            --fastp {input.fastp} \
            --kraken2 {input.kraken2} \
            --output-table {output.table} \
            --output-plot {output.plot} \
            2>&1 > {log.stderr}
        """


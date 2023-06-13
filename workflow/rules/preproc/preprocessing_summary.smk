# vim: syntax=python expandtab
# Summarize read counts passing through preprocessing.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

localrules:
    preprocessing_summary

# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
read_counts = OUTDIR/"preprocessing_read_counts.txt"

all_outputs.append(read_counts)


rule preprocessing_summary:
    """Summarize read counts in preprocessing steps"""
    input:
        # Host removal fastq output is used as a proxy to indicate of when all
        # quality control and host removal steps are completed
        expand(OUTDIR/"host_removal/{sample}_1.fq.gz", sample=SAMPLES),
    output:
        table=report(OUTDIR/"preprocessing_read_counts.txt",
            category="Preprocessing",
            caption="../../report/preprocessing_summary.rst"),
    log:
        stdout=LOGDIR/"preprocessing_summary.log",
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    threads: 1
    params:
        fastp_arg=lambda w: f"--fastp {LOGDIR}/fastp/*.fastp.json" if config["qc_reads"] else "",
        kraken2_arg=lambda w: f"--kraken2 {LOGDIR}/host_removal/*.kraken2.log" if config["host_removal"]["kraken2"] else "",
        bowtie2_arg=lambda w: f"--bowtie2 {LOGDIR}/host_removal/*.samtools.fastq.log" if config["host_removal"]["bowtie2"] else "",
    shell:
        """
        workflow/scripts/preprocessing_summary.py \
            {params.fastp_arg} \
            {params.kraken2_arg} \
            {params.bowtie2_arg} \
            --output-table {output.table} \
            > {log.stdout}
        """


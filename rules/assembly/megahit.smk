# vim: syntax=python expandtab
# Assembly of metagenomic samples using MEGAHIT
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:

if config["assembly"]["megahit"]:
    # Add MEGAHIT output to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    megahit_output = expand(str(OUTDIR/"assembly/megahit/{sample}/{sample}.contigs.fa"), sample=SAMPLES)
    all_outputs.extend(megahit_output)

    citations.add((
       "Li, D., Liu, C-M., Luo, R., Sadakane, K., and Lam, T-W., (2015).",
        "MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.",
        "Bioinformatics, doi: 10.1093/bioinformatics/btv033 [PMID: 25609793]."
    ))

mh_config = config["megahit"]
rule megahit:
    """Metagenome assembly using MEGAHIT"""
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        contigs=OUTDIR/"assembly/megahit/{sample}/{sample}.contigs.fa",
    log:
        stdout=str(LOGDIR/"assembly/megahit/{sample}.megahit.stdout.log"),
        stderr=str(LOGDIR/"assembly/megahit/{sample}.megahit.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    threads:
        20
    params:
        presets=mh_config["presets"],
        extra=mh_config["extra"],
    shell:
        """
        megahit \
            -1 {input.read1} \
            -2 {input.read2} \
            --num-cpu-threads {threads} \
            --out-dir assembly/megahit/{wildcards.sample} \
            > {log.stdout} \
            2> {log.stderr}
        """


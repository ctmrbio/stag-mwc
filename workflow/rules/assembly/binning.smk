# vim: syntax=python expandtab
# Assembly and binning of metagenomic assemblies
from pathlib import Path

from snakemake.exceptions import WorkflowError

a_conf = config["assembly_"]
b_conf = config["binning_"]

localrules:


if config["assembly"]["binning"]:
    bins = expand(f"{OUTDIR}/binning/{a_conf['assembler']}/{{sample}}/vamb", sample=SAMPLES)

    all_outputs.extend(bins)

    citations.add(publications["Vamb"])

rule binning:
    """Metagenomic binning using Vamb."""
    input:
        contigs=f"{OUTDIR}/metawrap/assembly/{mw_config['assembler']}/{{sample}}/final_assembly.fasta",
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        concoct_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/concoct_bins",
        #maxbin2_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/maxbin2_bins",
        metabat2_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/metabat2_bins",
    log:
        stdout=f"{LOGDIR}/metawrap/{{sample}}.binning.stdout.log",
        stderr=f"{LOGDIR}/metawrap/{{sample}}.binning.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:assembly"+singularity_branch_tag
    threads: 20
    params:
        outdir=lambda w: f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{w.sample}",
        universal=mw_config["universal"],
    shell:
        """
        gunzip -c {input.read1} > reads_1.fastq
        gunzip -c {input.read2} > reads_2.fastq
        metawrap binning \
            -a {input.contigs} \
            -o {params.outdir} \
            -t {threads} \
            -m 32 \
            --concoct \
            --metabat2 \
            --run-checkm \
            {params.universal} \
           reads_1.fastq reads_2.fastq \
            2> {log.stderr} \
            > {log.stdout}
        """


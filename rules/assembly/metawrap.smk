# vim: syntax=python expandtab
# Assembly and binning of metagenomic assemblies using MetaWRAP:
#   - MEGAHIT
#   - MaxBin2, CONCOCT, MetaBAT
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

mw_config = config["metawrap"]

localrules:

if config["assembly"]:
    assemblies = expand(f"{OUTDIR}/metawrap/assembly/{mw_config['assembler']}/{{sample}}/final_assembly.fasta", sample=SAMPLES)
    all_outputs.extend(assemblies)

    citations.add(publications["MetaWRAP"])
    citations.add(publications["MEGAHIT"])

if config["binning"]:
    primary_bins = expand(f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/concoct_bins", sample=SAMPLES)
    consolidated_bins = expand(f"{OUTDIR}/metawrap/consolidated_bins/{mw_config['assembler']}/{{sample}}/metawrap", sample=SAMPLES)
    blobology = expand(f"{OUTDIR}/metawrap/blobology/{mw_config['assembler']}/{{sample}}", sample=SAMPLES)
    all_outputs.extend(primary_bins)
    all_outputs.extend(consolidated_bins)
    all_outputs.extend(blobology)

    citations.add(publications["MetaWRAP"])
    citations.add(publications["CONCOCT"])
    citations.add(publications["MetaBAT2"])
    citations.add(publications["MaxBin2"])


rule assembly:
    """Metagenomic assembly via MetaWRAP."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        contigs=f"{OUTDIR}/metawrap/assembly/{mw_config['assembler']}/{{sample}}/final_assembly.fasta",
    log:
        stdout=f"{LOGDIR}/metawrap/{mw_config['assembler']}/{{sample}}.assembly.stdout.log",
        stderr=f"{LOGDIR}/metawrap/{mw_config['assembler']}/{{sample}}.assembly.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:assembly"
    threads:
        cluster_config["assembly"]["n"] if "assembly" in cluster_config else 20
    params:
        outdir=lambda w: f"{OUTDIR}/metawrap/assembly/{mw_config['assembler']}/{w.sample}",
        assembler=mw_config["assembler"],
        memory=140,
    shell:
        """
        metawrap assembly \
            -1 {input.read1} \
            -2 {input.read2} \
            -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            --{params.assembler} \
            2> {log.stderr} \
            > {log.stdout}
        """


rule binning:
    """Metagenomic binning using CONCOCT, MaxBin2, MetaBAT via MetaWRAP."""
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
    singularity:
        "shub://ctmrbio/stag-mwc:assembly"
    threads:
        cluster_config["binning"]["n"] if "binning" in cluster_config else 20
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

rule consolidate_bins:
    """Consolidate bins using MetaWRAP."""
    input:
        concoct_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/concoct_bins",
        #maxbin2_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/maxbin2_bins",
        metabat2_bins=f"{OUTDIR}/metawrap/binning/{mw_config['assembler']}/{{sample}}/metabat2_bins",
    output:
        metawrap_bins=f"{OUTDIR}/metawrap/consolidated_bins/{mw_config['assembler']}/{{sample}}/metawrap",
        metawrap_stats=f"{OUTDIR}/metawrap/consolidated_bins/{mw_config['assembler']}/{{sample}}/metawrap.stats",
    log:
        stdout=f"{LOGDIR}/metawrap/{{sample}}.consolidate_bins.stdout.log",
        stderr=f"{LOGDIR}/metawrap/{{sample}}.consolidate_bins.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:assembly"
    threads:
        cluster_config["consolidate_bins"]["n"] if "consolidate_bins" in cluster_config else 20
    params:
        outdir=lambda w: f"{OUTDIR}/metawrap/consolidated_bins/{w.sample}",
        minimum_completion=mw_config["minimum_completion"],
        maximum_contamination=mw_config["maximum_contamination"],
    shell:
        """
        metawrap bin_refinement \
            -o {params.outdir} \
            -t {threads} \
            -c {params.minimum_completion} \
            -x {params.maximum_contamination} \
            -A {input.concoct_bins} \
            -B {input.metabat2_bins} \
            2> {log.stderr} \
            > {log.stdout}
        """


rule blobology:
    """Visualize bins using MetaWRAP blobology module"""
    input:
        contigs=f"{OUTDIR}/metawrap/assembly/{mw_config['assembler']}/{{sample}}/final_assembly.fasta",
        metawrap_bins=f"{OUTDIR}/metawrap/consolidated_bins/{mw_config['assembler']}/{{sample}}/metawrap",
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        plots=f"{OUTDIR}/metawrap/blobology/{mw_config['assembler']}/{{sample}}",
    log:
        stdout=f"{LOGDIR}/metawrap/{{sample}}.blobology.stdout.log",
        stderr=f"{LOGDIR}/metawrap/{{sample}}.blobology.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/assembly.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:assembly"
    threads:
        cluster_config["blobology"]["n"] if "blobology" in cluster_config else 20
    params:
        outdir=lambda w: f"{OUTDIR}/metawrap/blobology/{w.sample}",
    shell:
        """
        metawrap blobology \
            -a {input.contigs} \
            -t {threads} \
            -o {params.outdir} \
            --bins {input.metawrap_bins} \
            {input.read1} {input.read2} \
            2> {log.stderr} \
            > {log.stdout}
        """


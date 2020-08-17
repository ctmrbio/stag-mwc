# vim: syntax=python expandtab tabstop=4
# AMRPlusPlus v2.0 snakemake integration
# Aron Arzoomand 12/8 2020

from pathlib import Path, PurePath
from snakemake.exceptions import WorkflowError



localrules:
    build_amr_index,


if config["antibiotic_resistance"]["amrplusplus"]:
    amrplusplus_outputs = f"{OUTDIR}/amrplusplus/ResistomeResults/AMR_analytics_matrix.csv"
    all_outputs.append(amrplusplus_outputs)

    citations.add(publications["AMRPlusPlus2"])


amrplusplus_config=config["amrplusplus"]


rule build_amr_index:
    """Get the Megares_2.0 database and annotations"""
    input:
        R1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        R2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
    output:
        megares_db=amrplusplus_config["databases"]["db"],
        megares_annot=amrplusplus_config["databases"]["annotation"],
    log:
        stdout=f"{LOGDIR}/amrplusplus/build_amr.index.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/build_amr.index.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        db_dir="scripts/amrplusplus/db",
    shell:
        """
        wget -O {output.megares_db} https://raw.githubusercontent.com/meglab-metagenomics/amrplusplus_v2/master/data/amr/megares_modified_database_v2.00.fasta 
        wget -O {output.megares_annot} https://raw.githubusercontent.com/meglab-metagenomics/amrplusplus_v2/master/data/amr/megares_modified_annotations_v2.00.csv
        bwa index {output.megares_db} \
        2> {log.stderr \
        > {log.stdout}
        """

rule align_to_amr:
    """BWA to align host filtered reads to Megares_2.0 database"""
    input:
        R1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        R2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
        megares_db=amrplusplus_config["databases"]["db"],
    output:
        alignment=f"{OUTDIR}/amrplusplus/AlignToAMR/{{sample}}.amr.alignment.sam"
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.align_to_amr.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.align_to_amr.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        megares_db=amrplusplus_config["databases"]["db"],
    shell:
        """
        bwa mem \
        {params.megares_db} \
        {input.R1} \
        {input.R2} \
        -t 10 \
        -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' > {output.alignment} 
        """

rule run_resistome:
    """Annotate alignments"""
    input:
        alignment=f"{OUTDIR}/amrplusplus/AlignToAMR/{{sample}}.amr.alignment.sam",
    output:
        gene=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.gene.tsv",
        group=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.group.tsv",
        mech=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.mechanism.tsv",
        klass=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.class.tsv",
        typ=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.type.tsv",
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.run_resistome.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.run_resistome.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        script=amrplusplus_config["resistome"]["script"],
        threshold=amrplusplus_config["resistome"]["threshold"],
        megares_db=amrplusplus_config["databases"]["db"],
        megares_annot=amrplusplus_config["databases"]["annotation"],
    shell:
        """
        {params.script} \
        -ref_fp {params.megares_db} \
        -annot_fp {params.megares_annot} \
        -sam_fp {input.alignment} \
        -gene_fp {output.gene} \
        -group_fp {output.group} \
        -mech_fp {output.mech} \
        -class_fp {output.klass} \
        -type_fp {output.typ} \
        -t {params.threshold} \
        2> {log.stderr \
        > {log.stdout}
        """

rule run_rarefaction:
    """Running rarefaction"""
    input:
        alignment=f"{OUTDIR}/amrplusplus/AlignToAMR/{{sample}}.amr.alignment.sam",
    output:
        gene=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.gene.tsv",
        group=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.group.tsv",
        mech=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.mechanism.tsv",
        klass=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.class.tsv",
        typ=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.type.tsv",
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.run_rarefaction.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.run_rarefaction.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        megares_db=amrplusplus_config["databases"]["db"],
        megares_annot=amrplusplus_config["databases"]["annotation"],
        script=amrplusplus_config["rarefaction"]["script"],
        min=amrplusplus_config["rarefaction"]["min"],
        max=amrplusplus_config["rarefaction"]["max"],
        skip=amrplusplus_config["rarefaction"]["skip"],
        samples=amrplusplus_config["rarefaction"]["samples"],
        threshold=amrplusplus_config["resistome"]["threshold"],
    shell:
        """
        {params.script} \
        -ref_fp {params.megares_db} \
        -sam_fp {input.alignment} \
        -annot_fp {params.megares_annot} \
        -gene_fp {output.gene} \
        -group_fp {output.group} \
        -mech_fp {output.mech} \
        -class_fp {output.klass} \
        -type_fp {output.typ} \
        -min {params.min} \
        -max {params.max} \
        -skip {params.skip} \
        -samples {params.samples} \
        -t {params.threshold} \
        2> {log.stderr \
        > {log.stdout}
        """

rule resistome_results:
    """Creating matrix of all gene hits"""
    input:
        gene_resistome=expand(f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.gene.tsv", sample=SAMPLES),
        gene_rarefaction=expand(f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.gene.tsv", sample=SAMPLES)
    output:
        AMR_matrix=f"{OUTDIR}/amrplusplus/ResistomeResults/AMR_analytics_matrix.csv" 
    log:
        stdout=f"{LOGDIR}/amrplusplus/resistome_results.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/resistome_results.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        script=amrplusplus_config["amr_long_to_wide"]
    shell:
        """
        python3 {params.script} \
        -i {input.gene_resistome} \
        -o {output.AMR_matrix} \
        2> {log.stderr \
        > {log.stdout}
        """

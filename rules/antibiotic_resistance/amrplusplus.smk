# vim: syntax=python expandtab tabstop=4
# AMRPlusPlus v2.0 snakemake integration
# Aron Arzoomand 2020-08-12

from pathlib import Path, PurePath
from snakemake.exceptions import WorkflowError



localrules:
    get_local_megares,
    build_amr_index,


amrplusplus_config=config["amrplusplus"]
if config["antibiotic_resistance"]["amrplusplus"]:
    if not Path(amrplusplus_config["megares"]["fasta"]).exists():
        err_message = "No database exists at: '{}' !\n".format(amrplusplus_config["megares"]["fasta"])
        err_message += "Specify the path in the amrplusplus section of config.yaml.\n"
        err_message += "If you want to run with the default megares_modified_database_v2.00.fasta, leave it as \"\""
        raise WorkflowError(err_message)
    if not Path(amrplusplus_config["megares"]["annotation"]).exists():
        err_message = "No annotations exists at: '{}' !\n".format(amrplusplus_config["megares"]["annotation"])
        err_message += "Specify the path in the amrplusplus section of config.yaml.\n"
        err_message += "If you want to run with the default megares_modified_annotations_v2.00.csv, leave it as \"\""
        raise WorkflowError(err_message)

    amrplusplus_outputs = f"{OUTDIR}/amrplusplus/ResistomeResults/AMR_analytics_matrix.csv"
    all_outputs.append(amrplusplus_outputs)

    citations.add(publications["AMRPlusPlus2"])


MEGARES_FASTA = "databases/amrplusplus/megares_modified_database_v2.00.fasta"
MEGARES_ANNOT = "databases/amrplusplus/megares_modified_annotations_v2.00.csv"


rule get_local_megares:
    """Get the Megares_2.0 database and annotations"""
    output:
        megares_db_local=f"{MEGARES_FASTA}",
        megares_annot_local=f"{MEGARES_ANNOT}",
    log:
        stdout=f"{LOGDIR}/amrplusplus/get_local_megares.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/get_local_megares.stderr.log",
    shadow:
        "shallow"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        megares_db=amrplusplus_config["megares"]["fasta"] or f"{MEGARES_FASTA}"
    shell:
        """
        wget -O {output.megares_db_local} https://raw.githubusercontent.com/meglab-metagenomics/amrplusplus_v2/master/data/amr/megares_modified_database_v2.00.fasta 
        wget -O {output.megares_annot_local} https://raw.githubusercontent.com/meglab-metagenomics/amrplusplus_v2/master/data/amr/megares_modified_annotations_v2.00.csv
        """

rule build_amr_index:
    """Build megares index"""
    input:
        megares_db=amrplusplus_config["megares"]["fasta"] or f"{MEGARES_FASTA}"
    output:
        amb=f"{amrplusplus_config['megares']['fasta']}.amb" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.amb",
        ann=f"{amrplusplus_config['megares']['fasta']}.ann" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.ann",
        bwt=f"{amrplusplus_config['megares']['fasta']}.bwt" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.bwt",
        pac=f"{amrplusplus_config['megares']['fasta']}.pac" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.pac",
        sa=f"{amrplusplus_config['megares']['fasta']}.sa" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.sa",
    log:
        stdout=f"{LOGDIR}/amrplusplus/build_amr.index.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/build_amr.index.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/amrplusplus.yaml"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    shell:
        """
        bwa index {input.megares_db} \
            2> {log.stderr} \
            > {log.stdout}
        """

rule align_to_amr:
    """BWA to align host filtered reads to Megares_2.0 database"""
    input:
        R1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        R2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
        amb=f"{amrplusplus_config['megares']['fasta']}.amb" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.amb",
        ann=f"{amrplusplus_config['megares']['fasta']}.ann" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.ann",
        bwt=f"{amrplusplus_config['megares']['fasta']}.bwt" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.bwt",
        pac=f"{amrplusplus_config['megares']['fasta']}.pac" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.pac",
        sa=f"{amrplusplus_config['megares']['fasta']}.sa" if amrplusplus_config["megares"]["fasta"] else f"{MEGARES_FASTA}.sa",
    output:
        alignment=f"{OUTDIR}/amrplusplus/AlignToAMR/{{sample}}.amr.alignment.sam"
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.align_to_amr.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.align_to_amr.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/amrplusplus.yaml"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        megares_db=amrplusplus_config["megares"]["fasta"] or f"{MEGARES_FASTA}"
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
        type=f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.type.tsv",
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.run_resistome.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.run_resistome.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/amrplusplus.yaml"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        script="scripts/amrplusplus/resistome",
        threshold=amrplusplus_config["resistome"]["threshold"],
        megares_db=amrplusplus_config["megares"]["fasta"] or f"{MEGARES_FASTA}",
        megares_annot=amrplusplus_config["megares"]["annotation"] or f"{MEGARES_ANNOT}",
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
            -type_fp {output.type} \
            -t {params.threshold} \
            2> {log.stderr} \
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
        type=f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.type.tsv",
    log:
        stdout=f"{LOGDIR}/amrplusplus/{{sample}}.run_rarefaction.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/{{sample}}.run_rarefaction.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/amrplusplus.yaml"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        megares_db=amrplusplus_config["megares"]["fasta"] or f"{MEGARES_FASTA}",
        megares_annot=amrplusplus_config["megares"]["annotation"] or f"{MEGARES_ANNOT}",
        script="scripts/amrplusplus/rarefaction",
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
            -type_fp {output.type} \
            -min {params.min} \
            -max {params.max} \
            -skip {params.skip} \
            -samples {params.samples} \
            -t {params.threshold} \
            2> {log.stderr} \
            > {log.stdout}
        """

rule resistome_results:
    """Creating matrix of all gene hits"""
    input:
        gene_resistome=expand(f"{OUTDIR}/amrplusplus/RunResistome/{{sample}}.gene.tsv", sample=SAMPLES),
        gene_rarefaction=expand(f"{OUTDIR}/amrplusplus/RunRarefaction/{{sample}}.gene.tsv", sample=SAMPLES),
    output:
        AMR_matrix=f"{OUTDIR}/amrplusplus/ResistomeResults/AMR_analytics_matrix.csv" 
    log:
        stdout=f"{LOGDIR}/amrplusplus/resistome_results.stdout.log",
        stderr=f"{LOGDIR}/amrplusplus/resistome_results.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/amrplusplus.yaml"
    singularity:
        "shub://meglab-metagenomics/amrplusplus_v2"
    params:
        script="scripts/amrplusplus/amr_long_to_wide.py"
    shell:
        """
        python3 {params.script} \
            -i {input.gene_resistome} \
            -o {output.AMR_matrix} \
            2> {log.stderr} \
            > {log.stdout}
        """

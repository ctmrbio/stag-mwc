# vim: syntax=python expandtab
# Functional profiling of metagenomic reads using HUMAnN

from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    download_humann_databases,
    normalize_humann_tables,
    join_humann_tables,
    regroup_humann_tables,
    rename_humann_tables,

h_config = config["humann"]
HMN_METHOD=h_config["norm_method"]
HMN_MODE=h_config["norm_mode"]
HMN_GRP=h_config["groups"]
HMN_RENAME=h_config["rename"]


if config["functional_profile"]["humann"]:
    if (not all([h_config["nucleotide_db"], h_config["protein_db"]]) 
        or not any([Path(h_config["nucleotide_db"]).is_dir(), Path(h_config["protein_db"]).is_dir()])):
        err_message = "Could not find HUMAnN nucleotide and protein databases at: '{}', '{}'!\n".format(h_config["nucleotide_db"], h_config["protein_db"])
        err_message += "Specify relevant paths in the humann section of config.yaml.\n"
        err_message += "If you do not want to run HUMAnN for functional profiling, set functional_profile:humann: False in config.yaml"
        raise WorkflowError(err_message)
    bt2_db_ext = ".1.bt2"

    merged_humann_tables = expand(f"{OUTDIR}/humann/all_samples.humann_{{output_type}}.tsv",
            output_type=("genefamilies", "pathcoverage", "pathabundance"))
    all_outputs.extend(merged_humann_tables)

    citations.add(publications["HUMAnN"])



rule download_humann_databases:
    """Download ChocoPhlAn and UniRef90 (diamond)"""
    output:
        f"{DBDIR}/humann/chocophlan",
        f"{DBDIR}/humann/uniref",
    log:
        f"{LOGDIR}/humann/database_download.log"
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    params:
        dbdir=config["dbdir"]+"/humann"
    shell:
        """
        humann_databases --download chocophlan full {params.dbdir} > {log}
        humann_databases --download uniref uniref90_diamond {params.dbdir} >> {log}
        """


rule humann:
    """Functional profiling using HUMAnN."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
        taxonomic_profile=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt",
    output:
        f"{OUTDIR}/humann/{{sample}}_genefamilies.tsv",
        f"{OUTDIR}/humann/{{sample}}_pathcoverage.tsv",
        f"{OUTDIR}/humann/{{sample}}_pathabundance.tsv",
    log:
        stdout=f"{LOGDIR}/humann/{{sample}}.humann.stdout.log",
        stderr=f"{LOGDIR}/humann/{{sample}}.humann.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads:
        cluster_config["humann"]["n"] if "humann" in cluster_config else 20
    resources:
        humann=1
    params:
        outdir=OUTDIR/"humann",
        nucleotide_db=h_config["nucleotide_db"],
        protein_db=h_config["protein_db"],
    shell:
        """
        cat {input.read1} {input.read2} > concat_input_reads.fq.gz

        humann \
            --input concat_input_reads.fq.gz \
            --output {params.outdir} \
            --nucleotide-database {params.nucleotide_db} \
            --protein-database {params.protein_db} \
            --output-basename {wildcards.sample} \
            --threads {threads} \
            --taxonomic-profile {input.taxonomic_profile} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule normalize_humann_tables:
    """Normalize abundance tables from HUMAnN."""
    input:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies.tsv",
        pathabundance=f"{OUTDIR}/humann/{{sample}}_pathabundance.tsv",
    output:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}.tsv",
        pathabundance=f"{OUTDIR}/humann/{{sample}}_pathabundance_{HMN_METHOD}.tsv",
    log:
        stdout=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}.stdout.log",
        stderr=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads: 
        1
    params:
        method=h_config["norm_method"],
        mode=h_config["norm_mode"],
    shell:
        """
        humann_renorm_table \
            --input {input.genefamilies} \
            --output {output.genefamilies} \
            --units {params.method} \
            --mode {params.mode} \
            > {log.stdout} \
            2> {log.stderr}
        humann_renorm_table \
            --input {input.pathabundance} \
            --output {output.pathabundance} \
            --units {params.method} \
            --mode {params.mode} \
            >> {log.stdout} \
            2>> {log.stderr}
        """


rule regroup_humann_tables:
    """regroup gene families"""
    input:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}.tsv",
    output:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}_{HMN_GRP}.tsv",
    log:
        stdout=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}_{HMN_GRP}.stdout.log",
        stderr=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}_{HMN_GRP}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads: 
        1
    params:
        groups=h_config["groups"]
    shell:
        """
        humann_regroup_table \
            --input {input.genefamilies} \
            --output {output.genefamilies} \
            --groups {params.groups} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule rename_humann_tables:
    """rename gene families"""
    input:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}_{HMN_GRP}.tsv",
    output:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}_{HMN_GRP}_{HMN_RENAME}.tsv",
    log:
        stdout=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}_{HMN_GRP}_{HMN_RENAME}.stdout.log",
        stderr=f"{LOGDIR}/humann/{{sample}}.humann_sample_{HMN_METHOD}_{HMN_GRP}_{HMN_RENAME}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads: 
        1
    params:
        rename=h_config["rename"]
    shell:
        """
        humann_rename_table \
            --input {input.genefamilies} \
            --output {output.genefamilies} \
            --names {params.rename} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule join_humann_tables:
    """Join abundance tables from HUMAnN."""
    input:
        genefamilies=expand(f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}_{HMN_GRP}_{HMN_RENAME}.tsv", sample=SAMPLES),
        pathabundance=expand(f"{OUTDIR}/humann/{{sample}}_pathabundance_{HMN_METHOD}.tsv", sample=SAMPLES),
        pathcoverage=expand(f"{OUTDIR}/humann/{{sample}}_pathcoverage.tsv", sample=SAMPLES),
    output:
        genefamilies=report(OUTDIR/"humann/all_samples.humann_genefamilies.tsv",
                category="Functional profiling",
                caption="../../report/humann_table.rst"),
        pathabundance=report(OUTDIR/"humann/all_samples.humann_pathabundance.tsv",
                category="Functional profiling",
                caption="../../report/humann_table.rst"),
        pathcoverage=report(OUTDIR/"humann/all_samples.humann_pathcoverage.tsv",
                category="Functional profiling",
                caption="../../report/humann_table.rst"),
    log:
        stdout=f"{LOGDIR}/humann/humann_join_tables.stdout.log",
        stderr=f"{LOGDIR}/humann/humann_join_tables.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    singularity:
        "shub://AroArz/singularity_playground:biobakery"
    threads: 
        1
    params:
        output_dir=OUTDIR/"humann",
        genefamilies=f"genefamilies_{HMN_METHOD}_{HMN_GRP}_{HMN_RENAME}",
        pathabundance=f"pathabundance_{HMN_METHOD}",
    shell:
        """
        humann_join_tables \
            --input {params.output_dir} \
            --output {output.genefamilies} \
            --file_name {params.genefamilies} \
            > {log.stdout} \
            2> {log.stderr}
        humann_join_tables \
            --input {params.output_dir} \
            --output {output.pathcoverage} \
            --file_name pathcoverage \
            >> {log.stdout} \
            2>> {log.stderr}
        humann_join_tables \
            --input {params.output_dir} \
            --output {output.pathabundance} \
            --file_name {params.pathabundance} \
            >> {log.stdout} \
            2>> {log.stderr}
        """

# vim: syntax=python expandtab
# Functional profiling of metagenomic reads using HUMAnN

from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    normalize_humann_tables,
    humann_join_tables,

h_config = config["humann"]
HMN_METHOD=h_config["norm_method"]
HMN_MODE=h_config["norm_mode"]
HMN_UTILITY=h_config["utility_db"]

if config["functional_profile"]["humann"]:
    if (not all([h_config["nucleotide_db"], h_config["protein_db"]]) 
        or not all([Path(h_config["nucleotide_db"]).is_dir(), Path(h_config["protein_db"]).is_dir()])):
        err_message = "Could not find HUMAnN nucleotide and protein databases at: '{}', '{}'!\n".format(h_config["nucleotide_db"], h_config["protein_db"])
        err_message += "Specify relevant paths in the humann section of config.yaml.\n"
        err_message += "If you do not want to run HUMAnN for functional profiling, set functional_profile:humann: False in config.yaml"
        err_message += "If you want to download the HUMAnN databases please see https://github.com/biobakery/humann"
        raise WorkflowError(err_message)

    merged_humann_tables = expand(f"{OUTDIR}/humann/all_samples.humann_{{output_type}}.txt",
            output_type=("genefamilies", "pathcoverage", "pathabundance"))
    all_outputs.extend(merged_humann_tables)

    citations.add(publications["HUMAnN"])


rule humann:
    """Functional profiling using HUMAnN."""
    input:
        read1=f"{OUTDIR}/host_removal/{{sample}}_1.fq.gz",
        read2=f"{OUTDIR}/host_removal/{{sample}}_2.fq.gz",
        taxonomic_profile=f"{OUTDIR}/metaphlan/{{sample}}.metaphlan.txt",
    output:
        genefamilies=f"{OUTDIR}/humann/{{sample}}_genefamilies.tsv",
        pathcoverage=f"{OUTDIR}/humann/{{sample}}_pathcoverage.tsv",
        pathabundance=f"{OUTDIR}/humann/{{sample}}_pathabundance.tsv",
    log:
        stdout=f"{LOGDIR}/humann/{{sample}}.humann.stdout.log",
        stderr=f"{LOGDIR}/humann/{{sample}}.humann.stderr.log",
        log=f"{LOGDIR}/humann/{{sample}}.humann.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    container:
        config["containers"]["humann"]
    threads: 20
    params:
        outdir=f"{OUTDIR}/humann/",
        nucleotide_db=h_config["nucleotide_db"],
        protein_db=h_config["protein_db"],
        extra=h_config["extra"],
    shell:
        """
        TEMPDIR="{resources.tmpdir}/{wildcards.sample}"
        mkdir -pv $TEMPDIR >> {log.stdout}
        mkdir -pv {params.outdir} >> {log.stdout}
        cat {input.read1} {input.read2} > $TEMPDIR/{wildcards.sample}_concat.fq.gz

        humann \
            --input $TEMPDIR/{wildcards.sample}_concat.fq.gz \
            --output $TEMPDIR \
            --nucleotide-database {params.nucleotide_db} \
            --protein-database {params.protein_db} \
            --output-basename {wildcards.sample} \
            --threads {threads} \
            --taxonomic-profile {input.taxonomic_profile} \
            {params.extra} \
            >> {log.stdout} \
            2> {log.stderr}

        cp $TEMPDIR/{wildcards.sample}*.tsv {params.outdir}
        cp $TEMPDIR/{wildcards.sample}_humann_temp/{wildcards.sample}.log {log.log}
        rm -rfv $TEMPDIR >> {log.stdout}
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
    container:
        config["containers"]["humann"]
    threads: 1
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

rule humann_join_tables:
    """Join normalized abundance tables from HUMAnN."""
    input:
        genefamilies=expand(f"{OUTDIR}/humann/{{sample}}_genefamilies_{HMN_METHOD}.tsv", sample=SAMPLES),
        pathabundance=expand(f"{OUTDIR}/humann/{{sample}}_pathabundance_{HMN_METHOD}.tsv", sample=SAMPLES),
    output:
        genefamilies=report(f"{OUTDIR}/humann/all_samples.humann_genefamilies.txt",
            category="Functional profiling",
            caption="../../report/humann_genefamilies.rst"),
        pathabundance=report(f"{OUTDIR}/humann/all_samples.humann_pathabundance.txt",
            category="Functional profiling",
            caption="../../report/humann_pathabundance.rst"),
        pathcoverage=report(f"{OUTDIR}/humann/all_samples.humann_pathcoverage.txt",
            category="Functional profiling",
            caption="../../report/humann_pathcoverage.rst"),
    log:
        stdout=f"{LOGDIR}/humann/humann_join_tables.stdout.log",
        stderr=f"{LOGDIR}/humann/humann_join_tables.stderr.log",
    shadow:
        "shallow"
    conda:
        "../../envs/humann.yaml"
    container:
        config["containers"]["humann"]
    threads: 
        1
    params:
        output_dir=f"{OUTDIR}/humann",
        genefamilies=f"genefamilies_{HMN_METHOD}",
        pathabundance=f"pathabundance_{HMN_METHOD}",
        pathcoverage=f"pathcoverage",
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
            --file_name {params.pathcoverage} \
            >> {log.stdout} \
            2>> {log.stderr}

        humann_join_tables \
            --input {params.output_dir} \
            --output {output.pathabundance} \
            --file_name {params.pathabundance} \
            >> {log.stdout} \
            2>> {log.stderr}
        """


# vim: syntax=python expandtab
# Functional profiling of metagenomic reads using HUMAnN2
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    download_humann2_databases,
    normalize_humann2_tables,
    join_humann2_tables,

h_config = config["humann2"]
if config["functional_profile"]["humann2"]:
    if not any([Path(h_config["nucleotide_db"]).exists(),
                Path(h_config["protein_db"]).exists()]):
        err_message = "Could not find HUMAnN2 nucleotide and protein databases at: '{}', '{}'!\n".format(h_config["nucleotide_db"], h_config["protein_db"])
        err_message += "Specify relevant paths in the humann2 section of config.yaml.\n"
        err_message += "Run 'snakemake download_humann2_databases' to download and build the default ChocoPhlAn and UniRef90 databases in '{dbdir}'\n".format(dbdir=DBDIR/"humann2")
        err_message += "If you do not want to run HUMAnN2 for functional profiling, set functional_profile:humann2: False in config.yaml"
        raise WorkflowError(err_message)
    bt2_db_ext = ".1.bt2"

    # Add HUMAnN2 output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    humann2_outputs = expand(str(OUTDIR/"humann2/{sample}_{output_type}.tsv"),
            sample=SAMPLES,
            output_type=("genefamilies", "pathcoverage", "pathabundance",
                        "genefamilies_{method}".format(method=h_config["norm_method"])))
    merged_humann2_tables = expand(str(OUTDIR/"humann2/all_samples.humann2_{output_type}.tsv"),
            output_type=("genefamilies", "pathcoverage", "pathabundance"))
    all_outputs.extend(humann2_outputs)
    all_outputs.extend(merged_humann2_tables)


rule download_humann2_databases:
    """Download ChocoPhlAn and UniRef90 (diamond)"""
    output:
        DBDIR/"humann2/chocophlan",
        DBDIR/"humann2/uniref",
    log:
        str(LOGDIR/"humann2/database_download.log")
    shadow:
        "shallow"
    conda:
        "../../envs/biobakery.yaml"
    params:
        dbdir=config["dbdir"]+"/humann2"
    shell:
        """
        humann2_databases --download chocophlan full {params.dbdir} > {log}
        humann2_databases --download uniref uniref90_diamond {params.dbdir} >> {log}
        """

# Ensure MPA output exists before running HUMAnN2
if not config["taxonomic_profile"]["metaphlan2"]:
    include: "../taxonomic_profiling/metaphlan2.smk"

rule humann2:
    """Functional profiling using HUMAnN2."""
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
        taxonomic_profile=OUTDIR/"metaphlan2/{sample}.metaphlan2.txt",
    output:
        OUTDIR/"humann2/{sample}_genefamilies.tsv",
        OUTDIR/"humann2/{sample}_pathcoverage.tsv",
        OUTDIR/"humann2/{sample}_pathabundance.tsv",
    log:
        stdout=str(LOGDIR/"humann2/{sample}.humann2.stdout.log"),
        stderr=str(LOGDIR/"humann2/{sample}.humann2.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/biobakery.yaml"
    threads:
        8
    params:
        outdir=OUTDIR/"humann2",
        nucleotide_db=h_config["nucleotide_db"],
        protein_db=h_config["protein_db"],
    shell:
        """
        cat {input.read1} {input.read2} > concat_input_reads.fq.gz \
        && \
        humann2 \
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


rule normalize_humann2_tables:
    """Normalize abundance tables from HUMAnN2."""
    input:
        genefamilies=OUTDIR/"humann2/{sample}_genefamilies.tsv",
        pathabundance=OUTDIR/"humann2/{sample}_pathabundance.tsv",
    output:
        genefamilies=OUTDIR/"humann2/{{sample}}_genefamilies_{method}.tsv".format(method=h_config["norm_method"]),
        pathabundance=OUTDIR/"humann2/{{sample}}_pathabundance_{method}.tsv".format(method=h_config["norm_method"]),
    log:
        stdout=str(LOGDIR/"humann2/{sample}.humann2_sample_normalize.stdout.log"),
        stderr=str(LOGDIR/"humann2/{sample}.humann2_sample_normalize.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/biobakery.yaml"
    threads: 
        1
    params:
        method=h_config["norm_method"],
        mode=h_config["norm_mode"],
    shell:
        """
        humann2_renorm_table \
            --input {input.genefamilies} \
            --output {output.genefamilies} \
            --units {params.method} \
            --mode {params.mode} \
            > {log.stdout} \
            2> {log.stderr}
        humann2_renorm_table \
            --input {input.pathabundance} \
            --output {output.pathabundance} \
            --units {params.method} \
            --mode {params.mode} \
            >> {log.stdout} \
            2>> {log.stderr}
        """


rule join_humann2_tables:
    """Join abundance tables from HUMAnN2."""
    input:
        genefamilies=expand(str(OUTDIR/"humann2/{{sample}}_genefamilies_{method}.tsv").format(method=h_config["norm_method"]),
            sample=SAMPLES),
        pathabundance=expand(str(OUTDIR/"humann2/{{sample}}_genefamilies_{method}.tsv").format(method=h_config["norm_method"]),
            sample=SAMPLES),
        pathcoverage=expand(str(OUTDIR/"humann2/{sample}_pathcoverage.tsv"), sample=SAMPLES),
    output:
        genefamilies=OUTDIR/"humann2/all_samples.humann2_genefamilies.tsv",
        pathabundance=OUTDIR/"humann2/all_samples.humann2_pathabundance.tsv",
        pathcoverage=OUTDIR/"humann2/all_samples.humann2_pathcoverage.tsv",
    log:
        stdout=str(LOGDIR/"humann2/humann2_join_tables.stdout.log"),
        stderr=str(LOGDIR/"humann2/humann2_join_tables.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/biobakery.yaml"
    threads: 
        1
    params:
        output_dir=OUTDIR/"humann2",
        genefamilies="genefamilies_{method}".format(method=h_config["norm_method"]),
        pathabundance="pathabundance_{method}".format(method=h_config["norm_method"]),
    shell:
        """
        humann2_join_tables \
            --input {params.output_dir} \
            --output {output.genefamilies} \
            --file_name {params.genefamilies} \
            > {log.stdout} \
            2> {log.stderr}
        humann2_join_tables \
            --input {params.output_dir} \
            --output {output.pathcoverage} \
            --file_name pathcoverage \
            >> {log.stdout} \
            2>> {log.stderr}
        humann2_join_tables \
            --input {params.output_dir} \
            --output {output.pathabundance} \
            --file_name {params.pathabundance} \
            >> {log.stdout} \
            2>> {log.stderr}
        """




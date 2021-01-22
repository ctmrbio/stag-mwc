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
    if (not all([h_config["nucleotide_db"], h_config["protein_db"]]) 
        or not any([Path(h_config["nucleotide_db"]).is_dir(), Path(h_config["protein_db"]).is_dir()])):
        err_message = "Could not find HUMAnN2 nucleotide and protein databases at: '{}', '{}'!\n".format(h_config["nucleotide_db"], h_config["protein_db"])
        err_message += "Specify relevant paths in the humann2 section of config.yaml.\n"
        err_message += "If you do not want to run HUMAnN2 for functional profiling, set functional_profile:humann2: False in config.yaml"
        raise WorkflowError(err_message)
    bt2_db_ext = ".1.bt2"

    # Add HUMAnN2 output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    humann2_outputs = expand(str(OUTDIR/"humann2/{sample}_{output_type}.txt"),
            sample=SAMPLES,
            output_type=("genefamilies", "pathcoverage", "pathabundance",
                        "genefamilies_{method}".format(method=h_config["norm_method"])))
    merged_humann2_tables = expand(str(OUTDIR/"humann2/all_samples.humann2_{output_type}.txt"),
            output_type=("genefamilies", "pathcoverage", "pathabundance"))
    all_outputs.extend(humann2_outputs)
    all_outputs.extend(merged_humann2_tables)

    citations.add(publications["HUMAnN2"])


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
        "../../envs/humann2.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:biobakery"
    params:
        dbdir=config["dbdir"]+"/humann2"
    shell:
        """
        humann2_databases --download chocophlan full {params.dbdir} > {log}
        humann2_databases --download uniref uniref90_diamond {params.dbdir} >> {log}
        """

rule humann2:
    """Functional profiling using HUMAnN2."""
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
        taxonomic_profile=OUTDIR/"metaphlan2/{sample}.metaphlan2.txt",
    output:
        OUTDIR/"humann2/{sample}_genefamilies.txt",
        OUTDIR/"humann2/{sample}_pathcoverage.txt",
        OUTDIR/"humann2/{sample}_pathabundance.txt",
    log:
        stdout=str(LOGDIR/"humann2/{sample}.humann2.stdout.log"),
        stderr=str(LOGDIR/"humann2/{sample}.humann2.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/humann2.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:biobakery"
    threads:
        cluster_config["humann2"]["n"] if "humann2" in cluster_config else 20
    resources:
        humann2=1
    params:
        outdir=OUTDIR/"humann2",
        nucleotide_db=h_config["nucleotide_db"],
        protein_db=h_config["protein_db"],
    shell:
        """
        # Convert MPA2 v2.96.1 output to something like MPA2 v2.7.7 output 
        # so it can be used with HUMAnN2, avoids StaG issue #138.
        # TODO: Remove this once HUMANn2 v2.9 is out.
        echo "#SampleID\t{wildcards.sample}" > mpa2_table-v2.7.7.txt
        sed '/#/d' {input.taxonomic_profile} \
            | cut -f1,3 \
            >> mpa2_table-v2.7.7.txt

        cat {input.read1} {input.read2} > concat_input_reads.fq.gz

        humann2 \
            --input concat_input_reads.fq.gz \
            --output {params.outdir} \
            --nucleotide-database {params.nucleotide_db} \
            --protein-database {params.protein_db} \
            --output-basename {wildcards.sample} \
            --threads {threads} \
            --taxonomic-profile mpa2_table-v2.7.7.txt \
            > {log.stdout} \
            2> {log.stderr}
        """


rule normalize_humann2_tables:
    """Normalize abundance tables from HUMAnN2."""
    input:
        genefamilies=OUTDIR/"humann2/{sample}_genefamilies.txt",
        pathabundance=OUTDIR/"humann2/{sample}_pathabundance.txt",
    output:
        genefamilies=OUTDIR/"humann2/{{sample}}_genefamilies_{method}.txt".format(method=h_config["norm_method"]),
        pathabundance=OUTDIR/"humann2/{{sample}}_pathabundance_{method}.txt".format(method=h_config["norm_method"]),
    log:
        stdout=str(LOGDIR/"humann2/{sample}.humann2_sample_normalize.stdout.log"),
        stderr=str(LOGDIR/"humann2/{sample}.humann2_sample_normalize.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/humann2.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:biobakery"
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
        genefamilies=expand(str(OUTDIR/"humann2/{{sample}}_genefamilies_{method}.txt").format(method=h_config["norm_method"]),
            sample=SAMPLES),
        pathabundance=expand(str(OUTDIR/"humann2/{{sample}}_genefamilies_{method}.txt").format(method=h_config["norm_method"]),
            sample=SAMPLES),
        pathcoverage=expand(str(OUTDIR/"humann2/{sample}_pathcoverage.txt"), sample=SAMPLES),
    output:
        genefamilies=report(OUTDIR/"humann2/all_samples.humann2_genefamilies.txt",
                category="Functional profiling",
                caption="../../report/humann2_table.rst"),
        pathabundance=report(OUTDIR/"humann2/all_samples.humann2_pathabundance.txt",
                category="Functional profiling",
                caption="../../report/humann2_table.rst"),
        pathcoverage=report(OUTDIR/"humann2/all_samples.humann2_pathcoverage.txt",
                category="Functional profiling",
                caption="../../report/humann2_table.rst"),
    log:
        stdout=str(LOGDIR/"humann2/humann2_join_tables.stdout.log"),
        stderr=str(LOGDIR/"humann2/humann2_join_tables.stderr.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/humann2.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:biobakery"
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




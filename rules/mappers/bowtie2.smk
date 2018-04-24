# Generic rules for alignment of reads to a reference database using Bowtie2
from snakemake.exceptions import WorkflowError
import os.path

bt2_db_extensions = (".1.bt2", ".1.bt2l")
if not any([os.path.isfile(config["bowtie2"]["db_prefix"]+ext) for ext in bt2_db_extensions]):
    err_message = "Bowtie2 index not found at: '{}'\n".format(config["bowtie2"]["db_prefix"])
    err_message += "Check path in config setting 'bowtie2:db_prefix'.\n"
    err_message += "If you want to skip mapping with bowtie2, set mappers:bowtie2:False in config.yaml."
    raise WorkflowError(err_message)
bt2_db_name = os.path.basename(config["bowtie2"]["db_prefix"])

# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
bowtie2_alignments = expand("{outdir}/bowtie2/{db_name}/{sample}.bam",
        outdir=config["outdir"],
        sample=SAMPLES,
        db_name=bt2_db_name)
bowtie2_stats = expand("{outdir}/bowtie2/{db_name}/{sample}.{stats}.txt",
        outdir=config["outdir"],
        sample=SAMPLES,
        stats=["covstats", "rpkm"],
        db_name=bt2_db_name)
counts_table = expand("{outdir}/bowtie2/{db_name}/all_samples.counts_table.tab",
        outdir=config["outdir"],
        db_name=bt2_db_name,
        sample=SAMPLES)
featureCounts = expand("{outdir}/bowtie2/{db_name}/all_samples.featureCounts{output_type}",
        outdir=config["outdir"],
        db_name=bt2_db_name,
        sample=SAMPLES,
        output_type=["", ".summary", ".table.tsv"])
all_outputs.extend(bowtie2_alignments)
all_outputs.extend(bowtie2_stats)
if config["bowtie2"]["counts_table"]["annotations"]:
    if not os.path.isfile(config["bowtie2"]["counts_table"]["annotations"]):
        err_message = "Bowtie2 counts_table annotations not found at: '{}'\n".format(config["bowtie2"]["counts_table"]["annotations"])
        err_message += "Check path in config setting 'bowtie2:counts_table:annotations'.\n"
        err_message += "If you want to skip the counts table summary for Bowtie2, set bowtie2:counts_table:annotations to '' in config.yaml."
        raise WorkflowError(err_message)
    all_outputs.extend(counts_table)
if config["bowtie2"]["featureCounts"]["annotations"]:
    if not os.path.isfile(config["bowtie2"]["featureCounts"]["annotations"]):
        err_message = "Bowtie2 featureCounts annotations not found at: '{}'\n".format(config["bowtie2"]["featureCounts"]["annotations"])
        err_message += "Check path in config setting 'bowtie2:featureCounts:annotations'.\n"
        err_message += "If you want to skip mapping with Bowtie2, set mappers:bowtie2:False in config.yaml."
        raise WorkflowError(err_message)
    all_outputs.extend(featureCounts)

rule bowtie2:
    """Align reads using Bowtie2."""
    input:
        sample=[config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
                config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz"]
    output:
        config["outdir"]+"/bowtie2/{db_name}/{{sample}}.bam".format(db_name=bt2_db_name)
    log:
        config["outdir"]+"/logs/bowtie2/{db_name}/{{sample}}.log".format(db_name=bt2_db_name)
    params:
        index=config["bowtie2"]["db_prefix"],
        extra=config["bowtie2"]["extra"],
    threads:
        8
    wrapper:
        "0.23.1/bio/bowtie2/align"


rule bowtie2_mapping_stats:
    """Summarize bowtie2 mapping statistics."""
    input:
        bam=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.bam".format(dbname=bt2_db_name)
    output:
        covstats=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.covstats.txt".format(dbname=bt2_db_name),
        rpkm=config["outdir"]+"/bowtie2/{dbname}/{{sample}}.rpkm.txt".format(dbname=bt2_db_name)
    log:
        config["outdir"]+"/logs/bowtie2/{dbname}/{{sample}}.pileup.log".format(dbname=bt2_db_name)
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
        pileup.sh \
            in={input.bam} \
            out={output.covstats} \
            rpkm={output.rpkm} \
            2> {log}
        """

rule bowtie2_counts_table:
    input:
        rpkms=expand(config["outdir"]+"/bowtie2/{dbname}/{sample}.rpkm.txt",
                dbname=bt2_db_name,
                sample=SAMPLES)
    output:
        counts=config["outdir"]+"/bowtie2/{dbname}/all_samples.counts_table.tab".format(dbname=bt2_db_name),
    log:
        config["outdir"]+"/logs/bowtie2/{dbname}/all_samples.counts_table.log".format(dbname=bt2_db_name)
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        1
    params:
        annotations=config["bowtie2"]["counts_table"]["annotations"]
    shell:
        """
        scripts/make_count_table.py \
            --annotations {params.annotations} \
            {input} \
            > {output} \
            2> {log}
        """


fc_config = config["bowtie2"]["featureCounts"]
rule bowtie2_featureCounts:
    input:
        bams=expand(config["outdir"]+"/bowtie2/{dbname}/{sample}.bam",
                dbname=bt2_db_name,
                sample=SAMPLES)
    output:
        counts=config["outdir"]+"/bowtie2/{dbname}/all_samples.featureCounts".format(dbname=bt2_db_name),
        counts_table=config["outdir"]+"/bowtie2/{dbname}/all_samples.featureCounts.table.tsv".format(dbname=bt2_db_name),
        summary=config["outdir"]+"/bowtie2/{dbname}/all_samples.featureCounts.summary".format(dbname=bt2_db_name),
    log:
        config["outdir"]+"/logs/bowtie2/{dbname}/all_samples.featureCounts.log".format(dbname=bt2_db_name)
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        4
    params:
        annotations=fc_config["annotations"],
        feature_type=lambda x: fc_config["feature_type"] if fc_config["feature_type"] else "gene",
        attribute_type=lambda x: fc_config["attribute_type"] if fc_config["attribute_type"] else "gene_id",
        extra=fc_config["extra"],
        dbname=bt2_db_name,
    shell:
        """
        featureCounts \
            -a {params.annotations} \
            -o {output.counts} \
            -t {params.feature_type} \
            -g {params.attribute_type} \
            -T {threads} \
            {params.extra} \
            {input.bams} \
            > {log} \
            2>> {log} \
        && \
        cut \
            -f1,7- \
            {output.counts}  \
            | sed '1d' \
            | sed 's|\t\w\+/bowtie2/{params.dbname}/|\t|g' \
            > {output.counts_table}
        """

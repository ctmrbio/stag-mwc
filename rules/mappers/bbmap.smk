# Rules for generic read mapping using BBMap
from snakemake.exceptions import WorkflowError
import os.path

if not os.path.isdir(os.path.join(config["bbmap"]["db_path"], "ref")):
    err_message = "BBMap index not found at: '{}'\n".format(config["bbmap"]["db_path"])
    err_message += "Check path in config setting 'bbmap:db_path'.\n"
    err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
    raise WorkflowError(err_message)
if not os.path.isfile(config["bbmap"]["featureCounts"]["annotations"]):
    err_message = "BBMap featureCounts annotations not found at: '{}'\n".format(config["bbmap"]["featureCounts"]["annotations"])
    err_message += "Check path in config setting 'bbmap:featureCounts:annotations'.\n"
    err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
    raise WorkflowError(err_message)


# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
bbmap_alignments = expand("{outdir}/bbmap/{db_name}/{sample}.{output_type}",
        outdir=config["outdir"],
        db_name=config["bbmap"]["db_name"],
        sample=SAMPLES,
        output_type=("sam.gz", "covstats.txt", "rpkm.txt"))
featureCounts = expand("{outdir}/bbmap/{db_name}/all_samples.featureCounts{output_type}",
        outdir=config["outdir"],
        db_name=config["bbmap"]["db_name"],
        sample=SAMPLES,
        output_type=["", ".summary", ".table.tsv"])
all_outputs.extend(bbmap_alignments)
all_outputs.extend(featureCounts)

bbmap_config = config["bbmap"]
bbmap_output_folder = config["outdir"]+"/bbmap/{db_name}/".format(db_name=bbmap_config["db_name"])
rule bbmap:
    """BBMap"""
    input:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        sam=bbmap_output_folder+"{sample}.sam.gz",
        covstats=bbmap_output_folder+"{sample}.covstats.txt",
        rpkm=bbmap_output_folder+"{sample}.rpkm.txt",
    log:
        stdout=config["outdir"]+"/logs/bbmap/{sample}.bbmap.stdout.log",
        stderr=config["outdir"]+"/logs/bbmap/{sample}.bbmap.statsfile.txt"
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        8
    params:
        db_path=bbmap_config["db_path"],
        min_id=bbmap_config["min_id"],
        extra=bbmap_config["extra"],
    shell:
        """
        bbmap.sh \
            threads={threads} \
            minid={params.min_id} \
            path={params.db_path} \
            in1={input.read1} \
            in2={input.read2} \
            out={output.sam} \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            {params.extra} \
            > {log.stdout} \
            2> {log.stderr}
        """


fc_config = bbmap_config["featureCounts"]
rule bbmap_featureCounts:
    input:
        bams=expand(config["outdir"]+"/bbmap/{dbname}/{sample}.sam.gz",
                dbname=bbmap_config["db_name"],
                sample=SAMPLES)
    output:
        counts=config["outdir"]+"/bbmap/{dbname}/all_samples.featureCounts".format(dbname=bbmap_config["db_name"]),
        counts_table=config["outdir"]+"/bbmap/{dbname}/all_samples.featureCounts.table.tsv".format(dbname=bbmap_config["db_name"]),
        summary=config["outdir"]+"/bbmap/{dbname}/all_samples.featureCounts.summary".format(dbname=bbmap_config["db_name"]),
    log:
        config["outdir"]+"/logs/bbmap/{dbname}/all_samples.featureCounts.log".format(dbname=bbmap_config["db_name"])
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        4
    params:
        annotations=fc_config["annotations"],
        feature_type=lambda _: fc_config["feature_type"] if fc_config["feature_type"] else "gene",
        attribute_type=lambda _: fc_config["attribute_type"] if fc_config["attribute_type"] else "gene_id",
        extra=fc_config["extra"],
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
            | sed 's|\t\w\+/bbmap/\w\+/|\t|g' \
            > {output.counts_table}
        """

# Rules for generic read mapping using BBMap
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    bbmap_counts_table,
    bbmap_featureCounts,

db_path = Path(config["bbmap"]["db_path"])
if not Path(db_path/"ref").exists():
    err_message = "BBMap index not found at: '{}'\n".format(db_path)
    err_message += "Check path in config setting 'bbmap:db_path'.\n"
    err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
    raise WorkflowError(err_message)

# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
bbmap_alignments = expand(str(OUTDIR/"bbmap/{db_name}/{sample}.{output_type}"),
        db_name=config["bbmap"]["db_name"],
        sample=SAMPLES,
        output_type=("sam.gz", "covstats.txt", "rpkm.txt"))
counts_table = expand(str(OUTDIR/"bbmap/{db_name}/all_samples.counts_table.tab"),
        db_name=config["bbmap"]["db_name"],
        sample=SAMPLES)
featureCounts = expand(str(OUTDIR/"bbmap/{db_name}/all_samples.featureCounts{output_type}"),
        db_name=config["bbmap"]["db_name"],
        sample=SAMPLES,
        output_type=["", ".summary", ".table.tsv"])
all_outputs.extend(bbmap_alignments)
if config["bbmap"]["counts_table"]["annotations"]:
    if not Path(config["bbmap"]["counts_table"]["annotations"]).exists():
        err_message = "BBMap counts table annotations not found at: '{}'\n".format(config["bbmap"]["counts_table"]["annotations"])
        err_message += "Check path in config setting 'bbmap:counts_table:annotations'.\n"
        err_message += "If you want to skip read counts summary for BBMap, set bbmap:counts_table:annotations to '' in config.yaml."
        raise WorkflowError(err_message)
    all_outputs.extend(counts_table)
if config["bbmap"]["featureCounts"]["annotations"]:
    if not Path(config["bbmap"]["featureCounts"]["annotations"]).exists():
        err_message = "BBMap featureCounts annotations not found at: '{}'\n".format(config["bbmap"]["featureCounts"]["annotations"])
        err_message += "Check path in config setting 'bbmap:featureCounts:annotations'.\n"
        err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
        raise WorkflowError(err_message)
    all_outputs.extend(featureCounts)

bbmap_config = config["bbmap"]
bbmap_output_folder = OUTDIR/"bbmap/{db_name}".format(db_name=bbmap_config["db_name"])
rule bbmap:
    """BBMap"""
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        sam=bbmap_output_folder/"{sample}.sam.gz",
        covstats=bbmap_output_folder/"{sample}.covstats.txt",
        rpkm=bbmap_output_folder/"{sample}.rpkm.txt",
    log:
        stdout=str(LOGDIR/"bbmap/{sample}.bbmap.stdout.log"),
        stderr=str(LOGDIR/"bbmap/{sample}.bbmap.statsfile.txt"),
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


rule bbmap_counts_table:
    input:
        rpkms=expand(str(OUTDIR/"bbmap/{dbname}/{sample}.rpkm.txt"),
                dbname=bbmap_config["db_name"],
                sample=SAMPLES)
    output:
        counts=OUTDIR/"bbmap/{dbname}/all_samples.counts_table.tab".format(dbname=bbmap_config["db_name"]),
    log:
        str(LOGDIR/"bbmap/{dbname}/all_samples.counts_table.log".format(dbname=bbmap_config["db_name"]))
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        1
    params:
        annotations=config["bbmap"]["counts_table"]["annotations"]
    shell:
        """
        scripts/make_count_table.py \
            --annotations {params.annotations} \
            {input} \
            > {output} \
            2> {log}
        """


fc_config = bbmap_config["featureCounts"]
rule bbmap_featureCounts:
    input:
        bams=expand(str(OUTDIR/"bbmap/{dbname}/{sample}.sam.gz"),
                dbname=bbmap_config["db_name"],
                sample=SAMPLES)
    output:
        counts=OUTDIR/"bbmap/{dbname}/all_samples.featureCounts".format(dbname=bbmap_config["db_name"]),
        counts_table=OUTDIR/"bbmap/{dbname}/all_samples.featureCounts.table.tsv".format(dbname=bbmap_config["db_name"]),
        summary=OUTDIR/"bbmap/{dbname}/all_samples.featureCounts.summary".format(dbname=bbmap_config["db_name"]),
    log:
        str(LOGDIR/"bbmap/{dbname}/all_samples.featureCounts.log".format(dbname=bbmap_config["db_name"]))
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
        dbname=bbmap_config["db_name"],
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
            | sed 's|\t\w\+/bbmap/{params.dbname}/|\t|g' \
            > {output.counts_table}
        """

# vim: syntax=python expandtab
# Rules for generic read mapping using BBMap
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

for bbmap_config in config["bbmap"]:
    db_name = bbmap_config["db_name"]
    bbmap_output_folder = OUTDIR/"bbmap/{db_name}".format(db_name=db_name)
    bbmap_logdir = LOGDIR/"bbmap/{db_name}".format(db_name=db_name)

    db_path = Path(bbmap_config["db_path"])
    if config["mappers"]["bbmap"]:
        if not Path(db_path/"ref").exists():
            err_message = "BBMap index not found at: '{}'\n".format(db_path)
            err_message += "Check path in config setting 'bbmap:db_path'.\n"
            err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
            raise WorkflowError(err_message)

        # Add final output files from this module to 'all_outputs' from the main
        # Snakefile scope. SAMPLES is also from the main Snakefile scope.
        bbmap_alignments = expand(str(OUTDIR/"bbmap/{db_name}/{sample}.{output_type}"),
                db_name=db_name,
                sample=SAMPLES,
                output_type=("sam.gz", "covstats.txt", "rpkm.txt"))
        counts_table = expand(str(OUTDIR/"bbmap/{db_name}/counts.{column}.tsv"),
                db_name=db_name,
                column=map(str.strip, bbmap_config["counts_table"]["columns"].split(",")))
        featureCounts = expand(str(OUTDIR/"bbmap/{db_name}/all_samples.featureCounts{output_type}"),
                db_name=db_name,
                sample=SAMPLES,
                output_type=["", ".summary", ".table.txt"])
        all_outputs.extend(bbmap_alignments)

        if bbmap_config["counts_table"]["annotations"]:
            if not Path(bbmap_config["counts_table"]["annotations"]).exists():
                err_message = "BBMap counts table annotations not found at: '{}'\n".format(bbmap_config["counts_table"]["annotations"])
                err_message += "Check path in config setting 'bbmap:counts_table:annotations'.\n"
                err_message += "If you want to skip read counts summary for BBMap, set bbmap:counts_table:annotations to '' in config.yaml."
                raise WorkflowError(err_message)
            all_outputs.extend(counts_table)
        if bbmap_config["featureCounts"]["annotations"]:
            if not Path(bbmap_config["featureCounts"]["annotations"]).exists():
                err_message = "BBMap featureCounts annotations not found at: '{}'\n".format(bbmap_config["featureCounts"]["annotations"])
                err_message += "Check path in config setting 'bbmap:featureCounts:annotations'.\n"
                err_message += "If you want to skip mapping with BBMap, set mappers:bbmap:False in config.yaml."
                raise WorkflowError(err_message)
            all_outputs.extend(featureCounts)
            citations.add(publications["featureCount"])

        citations.add(publications["BBMap"])

    rule:
        f"""BBMap reads to {db_name}"""
        name: f"bbmap_{db_name}"
        input:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
        output:
            sam=bbmap_output_folder/"{sample}.sam.gz",
            covstats=bbmap_output_folder/"{sample}.covstats.txt",
            rpkm=bbmap_output_folder/"{sample}.rpkm.txt",
        log:
            stdout=str(bbmap_logdir/"{sample}.bbmap.stdout.log"),
            stderr=str(bbmap_logdir/"{sample}.bbmap.statsfile.txt"),
        message:
            "Mapping {{wildcards.sample}} to {db_name} using BBMap".format(db_name=db_name)
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 8
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


    if bbmap_config["counts_table"]["annotations"] and not bbmap_config["counts_table"]["columns"]:
        raise WorkflowError("Must define annotation column(s) for count table production!")

    rule:
        f"""Summarize read counts for {db_name}"""
        name: f"bbmap_counts_{db_name}"
        input:
            rpkms=expand(str(OUTDIR/"bbmap/{db_name}/{sample}.rpkm.txt"),
                    db_name=db_name,
                    sample=SAMPLES)
        output:
            expand(str(OUTDIR/"bbmap/{db_name}/counts.{column}.tsv"),
                    db_name=db_name,
                    column=map(str.strip, bbmap_config["counts_table"]["columns"].split(","))
            )
        log:
            str(bbmap_logdir/"counts.log")
        message:
            "Summarizing read counts for {db_name}".format(db_name=db_name)
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 1
        params:
            annotations=bbmap_config["counts_table"]["annotations"],
            columns=bbmap_config["counts_table"]["columns"],
            outdir=OUTDIR/"bbmap/{db_name}/".format(db_name=db_name),
        shell:
            """
            workflow/scripts/make_count_table.py \
                --annotation-file {params.annotations} \
                --columns {params.columns} \
                --outdir {params.outdir} \
                {input} \
                2> {log}
            """


    fc_config = bbmap_config["featureCounts"]
    rule:
        f"""Summarize feature counts for {db_name}"""
        name: f"bbmap_feature_counts_{db_name}"
        input:
            bams=expand(str(OUTDIR/"bbmap/{db_name}/{sample}.sam.gz"),
                    db_name=db_name,
                    sample=SAMPLES)
        output:
            counts=OUTDIR/"bbmap/{db_name}/all_samples.featureCounts".format(db_name=db_name),
            counts_table=OUTDIR/"bbmap/{db_name}/all_samples.featureCounts.table.txt".format(db_name=db_name),
            summary=OUTDIR/"bbmap/{db_name}/all_samples.featureCounts.summary".format(db_name=db_name),
        log:
            str(bbmap_logdir/"all_samples.featureCounts.log")
        message:
            "Summarizing feature counts for {db_name}".format(db_name=db_name)
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads: 4
        params:
            annotations=fc_config["annotations"],
            feature_type=lambda _: fc_config["feature_type"] if fc_config["feature_type"] else "gene",
            attribute_type=lambda _: fc_config["attribute_type"] if fc_config["attribute_type"] else "gene_id",
            extra=fc_config["extra"],
            dbname=db_name,
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
                2>> {log}
            cut \
                -f1,7- \
                {output.counts}  \
                | sed '1d' \
                | sed 's|\t\w\+/bbmap/{params.dbname}/|\t|g' \
                > {output.counts_table}
            """

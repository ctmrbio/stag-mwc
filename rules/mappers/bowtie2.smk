# Generic rules for alignment of reads to a reference database using Bowtie2
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

bt2_db_extensions = (".1.bt2", ".1.bt2l")
for bt2_config in config["bowtie2"]:
    if not any([Path(bt2_config["db_prefix"]).with_suffix(ext) for ext in bt2_db_extensions]):
        err_message = "Bowtie2 index not found at: '{}'\n".format(bt2_config["db_prefix"])
        err_message += "Check path in config setting 'bowtie2:db_prefix'.\n"
        err_message += "If you want to skip mapping with bowtie2, set mappers:bowtie2:False in config.yaml."
        raise WorkflowError(err_message)
    bt2_db_name = Path(bt2_config["db_prefix"]).name

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    bowtie2_alignments = expand(str(OUTDIR/"bowtie2/{db_name}/{sample}.bam"),
            sample=SAMPLES,
            db_name=bt2_db_name)
    bowtie2_stats = expand(str(OUTDIR/"bowtie2/{db_name}/{sample}.{stats}.txt"),
            sample=SAMPLES,
            stats=["covstats", "rpkm"],
            db_name=bt2_db_name)
    counts_table = expand(str(OUTDIR/"bowtie2/{db_name}/all_samples.counts_table.tab"),
            db_name=bt2_db_name,
            sample=SAMPLES)
    featureCounts = expand(str(OUTDIR/"bowtie2/{db_name}/all_samples.featureCounts{output_type}"),
            db_name=bt2_db_name,
            sample=SAMPLES,
            output_type=["", ".summary", ".table.tsv"])
    all_outputs.extend(bowtie2_alignments)
    all_outputs.extend(bowtie2_stats)
    if bt2_config["counts_table"]["annotations"]:
        if not Path(bt2_config["counts_table"]["annotations"]).exists():
            err_message = "Bowtie2 counts_table annotations not found at: '{}'\n".format(bt2_config["counts_table"]["annotations"])
            err_message += "Check path in config setting 'bowtie2:counts_table:annotations'.\n"
            err_message += "If you want to skip the counts table summary for Bowtie2, set bowtie2:counts_table:annotations to '' in config.yaml."
            raise WorkflowError(err_message)
        all_outputs.extend(counts_table)
    if bt2_config["featureCounts"]["annotations"]:
        if not Path(bt2_config["featureCounts"]["annotations"]).exists():
            err_message = "Bowtie2 featureCounts annotations not found at: '{}'\n".format(bt2_config["featureCounts"]["annotations"])
            err_message += "Check path in config setting 'bowtie2:featureCounts:annotations'.\n"
            err_message += "If you want to skip mapping with Bowtie2, set mappers:bowtie2:False in config.yaml."
            raise WorkflowError(err_message)
        all_outputs.extend(featureCounts)

    rule:
        """Align reads using Bowtie2."""
        input:
            sample=[OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
                    OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz"]
        output:
            OUTDIR/"bowtie2/{db_name}/{{sample}}.bam".format(db_name=bt2_db_name)
        log:
            str(LOGDIR/"bowtie2/{db_name}/{{sample}}.log".format(db_name=bt2_db_name))
        message:
            "Mapping {{wildcards.sample}} to {db_name} using Bowtie2".format(db_name=bt2_db_name)
        params:
            index=bt2_config["db_prefix"],
            extra=bt2_config["extra"],
        threads:
            8
        wrapper:
            "0.23.1/bio/bowtie2/align"


    rule:
        """Summarize bowtie2 mapping statistics."""
        input:
            bam=OUTDIR/"bowtie2/{db_name}/{{sample}}.bam".format(db_name=bt2_db_name)
        output:
            covstats=OUTDIR/"bowtie2/{db_name}/{{sample}}.covstats.txt".format(db_name=bt2_db_name),
            rpkm=OUTDIR/"bowtie2/{db_name}/{{sample}}.rpkm.txt".format(db_name=bt2_db_name)
        log:
            str(LOGDIR/"bowtie2/{db_name}/{{sample}}.pileup.log".format(db_name=bt2_db_name))
        message:
            "Summarizing bowtie2 mapping statistics against {db_name}.".format(db_name=bt2_db_name)
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

    rule:
        """Create count table for MetaPhlAn2 mappings."""
        input:
            rpkms=expand(str(OUTDIR/"bowtie2/{db_name}/{sample}.rpkm.txt"),
                    db_name=bt2_db_name,
                    sample=SAMPLES)
        output:
            counts=OUTDIR/"bowtie2/{db_name}/all_samples.counts_table.tab".format(db_name=bt2_db_name),
        log:
            str(LOGDIR/"bowtie2/{db_name}/all_samples.counts_table.log".format(db_name=bt2_db_name))
        message:
            "Creating count table for mappings to {db_name}".format(db_name=bt2_db_name)
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            1
        params:
            annotations=bt2_config["counts_table"]["annotations"]
        shell:
            """
            scripts/make_count_table.py \
                --annotations {params.annotations} \
                {input} \
                > {output} \
                2> {log}
            """


    fc_config = bt2_config["featureCounts"]
    rule:
        """Summarize featureCounts for Bowtie2 mappings."""
        input:
            bams=expand(str(OUTDIR/"bowtie2/{db_name}/{sample}.bam"),
                    db_name=bt2_db_name,
                    sample=SAMPLES)
        output:
            counts=OUTDIR/"bowtie2/{db_name}/all_samples.featureCounts".format(db_name=bt2_db_name),
            counts_table=OUTDIR/"bowtie2/{db_name}/all_samples.featureCounts.table.tsv".format(db_name=bt2_db_name),
            summary=OUTDIR/"bowtie2/{db_name}/all_samples.featureCounts.summary".format(db_name=bt2_db_name),
        log:
            str(LOGDIR/"bowtie2/{db_name}/all_samples.featureCounts.log".format(db_name=bt2_db_name))
        message:
            "Summarizing feature counts for {db_name} mappings.".format(db_name=bt2_db_name)
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

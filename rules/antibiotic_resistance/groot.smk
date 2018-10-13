# Rules for antibiotic resistance profile groot
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    create_groot_index

groot_db_path = Path(config["groot"]["index"])
if not Path(groot_db_path).exists():
    err_message = "No groot database found at: '{}'!\n".format(groot_db_path)
    err_message += "Specify the DB path in the groot section of config.yaml.\n"
    err_message += "Run 'snakemake create_groot_index' to download and build a groot index in '{dbdir}'\n".format(dbdir=DBDIR/"groot")
    err_message += "If you do not want to run groot to create antibiotic resistance profiles, set antibiotic_resistance: False in config.yaml"
    raise WorkflowError(err_message)

groot_outputs = expand(str(OUTDIR/"groot/{sample}/{sample}.{output_type}"),
        sample=SAMPLES,
        output_type=("groot_aligned.bam", "groot_report.tsv"))
all_outputs.extend(groot_outputs)

groot_config = config["groot"]
rule create_groot_index:
    """Create groot index."""
    output:
        dbdir=DBDIR/"groot/{}".format(groot_config["db"]),
        index=DBDIR/"groot/{db}_index".format(db=groot_config["db"]),
    log:
        str(LOGDIR/"groot/groot.get_index.log")
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"groot/",
        db=groot_config["db"],
        read_length=125,
    threads:
        4
    shell:
        """
        groot get \
            --database {params.db} \
            --out {params.dbdir}/{params.db} \
            --processors {threads} \
            --logFile {log}
        groot index \
            --msaDir {params.dbdir}/{params.db} \
            --readLength {params.read_length} \
            --outDir {output.index} \
            --processors {threads} \
            --logFile {log}
        """


rule groot_align:
    """Align reads to groot index."""
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        bam=OUTDIR/"groot/{sample}/{sample}.groot_aligned.bam",
        report=OUTDIR/"groot/{sample}/{sample}.groot_report.tsv",
        plots=directory(OUTDIR/"groot/{sample}/groot-plots"),
        graphs=directory(OUTDIR/"groot/{sample}/groot-graphs"),
    log:
        align=str(LOGDIR/"groot/{sample}.groot_align.log"),
        report=str(LOGDIR/"groot/{sample}.groot_report.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        8
    params:
        index=groot_config["index"]
    shell:
        """
        groot align \
            --fastq {input.read1},{input.read2} \
            --graphDir {output.graphs} \
            --indexDir {params.index} \
            --processors {threads} \
            --logFile {log.align} \
            > {output.bam}
        groot report \
            --bamFile {output.bam} \
            --lowCov \
            --plotCov \
            --processors {threads} \
            --logFile {log.report} \
            > {output.report}
        mv groot-plots {output.plots}
        """

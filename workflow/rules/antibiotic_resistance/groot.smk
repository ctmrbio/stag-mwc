# vim: syntax=python expandtab
# Rules for antibiotic resistance profile groot
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path
from snakemake.exceptions import WorkflowError

localrules:
    groot_report


groot_db_path = Path(config["groot"]["index_dir"])
if config["antibiotic_resistance"]["groot"]:
    if not Path(groot_db_path).exists():
        err_message = "No groot database found at: '{}'!\n".format(groot_db_path)
        err_message += "Specify the DB path in the groot section of config.yaml.\n"
        err_message += "If you do not want to run groot to create antibiotic resistance profiles, set antibiotic_resistance: False in config.yaml"
        raise WorkflowError(err_message)

    groot_outputs = expand(str(OUTDIR/"groot/{sample}/{sample}.{output_type}"),
            sample=SAMPLES,
            output_type=("groot_aligned.bam", "groot_report.txt"))
    all_outputs.extend(groot_outputs)

    citations.add(publications["GROOT"])

groot_config = config["groot"]

rule groot_align:
    """Align reads to groot index."""
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        bam=OUTDIR/"groot/{sample}/{sample}.groot_aligned.bam",
        graphs=directory(OUTDIR/"groot/{sample}/groot-graphs"),
    log:
        reformat=str(LOGDIR/"groot/{sample}.reformat.log"),
        align=str(LOGDIR/"groot/{sample}.groot_align.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "docker://quay.io/biocontainers/groot:1.1.2--hc0beb16_0"
    threads: 8
    params:
        index_dir=groot_config["index_dir"],
        minlength=groot_config["minlength"],
        maxlength=groot_config["maxlength"],
    shell:
        """
        groot align \
            --fastq {input.read1},{input.read2} \
            --graphDir {output.graphs} \
            --indexDir {params.index_dir} \
            --processors {threads} \
            --log {log.align} \
            > {output.bam}
        """

rule groot_report:
    """Report and plot ARG coverage from GROOT"""
    input:
        bam=OUTDIR/"groot/{sample}/{sample}.groot_aligned.bam",
    output:
        report=OUTDIR/"groot/{sample}/{sample}.groot_report.txt",
    log:
        report=str(LOGDIR/"groot/{sample}.groot_report.log"),
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "docker://quay.io/biocontainers/groot:1.1.2--hc0beb16_0"
    threads: 1
    params:
        covcutoff=groot_config["covcutoff"],
        lowcov=lambda _: "--lowCov" if groot_config["lowcov"] else ""
    shell:
        """
        groot report \
            --bamFile {input.bam} \
            --covCutoff {params.covcutoff} \
            {params.lowcov} \
            --processors {threads} \
            --log {log.report} \
            > {output.report}
        """


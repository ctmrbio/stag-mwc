# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Centrifuge
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    download_centrifuge_database

centrifuge_db_ext = ".1.cf"
if config["taxonomic_profile"]["centrifuge"]:
    if not Path(config["centrifuge"]["db_prefix"]).with_suffix(centrifuge_db_ext).exists():
        err_message = "No Centrifuge database found at: '{}'!\n".format(config["centrifuge"]["db_prefix"])
        err_message += "Specify Centrifuge database prefix in the Centrifuge section of config.yaml.\n"
        err_message += "Run 'snakemake download_centrifuge_database' to download a copy into '{dbdir}'\n".format(dbdir=DBDIR/"centrifuge")
        err_message += "If you do not want to run Centrifuge for taxonomic profiling, set centrifuge: False in config.yaml"
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    centrifuge = expand(str(OUTDIR/"centrifuge/{sample}.{output_type}.tsv"),
            sample=SAMPLES,
            output_type=("centrifuge", "centrifuge_report"))
    all_outputs.extend(centrifuge)

    citations.append((
        "Kim, D., et al. (2016).",
        "Centrifuge: rapid and sensitive classification of metagenomic sequences.",
        "Genome research."
    ))


rule download_centrifuge_database:
    output:
        db=[DBDIR/"centrifuge/p+h+v.{n}.cf".format(n=num) for num in (1,2,3)]
    shadow:
        "shallow"
    params:
        dbdir=DBDIR/"centrifuge"
    shell:
        """
        cd {params.dbdir}
        wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/data/p+h+v.tar.gz \
        && \
        tar -xf p+h+v.tar.gz \
        && \
        rm p+h+v.tar.gz
        """


cf_config = config["centrifuge"]
rule centrifuge:
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        classifications=OUTDIR/"centrifuge/{sample}.centrifuge.tsv",
        report=OUTDIR/"centrifuge/{sample}.centrifuge_report.tsv",
    log:
        str(LOGDIR/"centrifuge/{sample}.centrifuge.log")
    shadow: 
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    threads:
        4
    params:
        db_prefix=cf_config["db_prefix"],
    shell:
        """
        centrifuge \
            -x {params.db_prefix} \
            -1 {input.read1} \
            -2 {input.read2} \
            -S {output.classifications} \
            --report-file {output.report} \
            --threads {threads} \
            2> {log}
        """


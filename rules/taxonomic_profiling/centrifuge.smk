# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Centrifuge
from snakemake.exceptions import WorkflowError
import os.path

centrifuge_db_ext = ".1.cf"
if not os.path.isfile(config["centrifuge"]["db_prefix"]+centrifuge_db_ext):
    err_message = "No Centrifuge database found at: {}!\n".format(config["centrifuge"]["db_prefix"])
    err_message += "Specify Centrifuge database prefix in the Centrifuge section of config.yaml.\n"
    err_message += "Run 'snakemake download_centrifuge_database' to download a copy into '{dbdir}/centrifuge'\n".format(dbdir=config["dbdir"])
    err_message += "If you do not want to run Centrifuge for taxonomic profiling, set centrifuge: False in config.yaml"
    raise WorkflowError(err_message)

# Add final output files from this module to 'all_outputs' from the main
# Snakefile scope. SAMPLES is also from the main Snakefile scope.
centrifuge = expand("{outdir}/centrifuge/{sample}.{output_type}.tsv",
        outdir=config["outdir"],
        sample=SAMPLES,
        output_type=("centrifuge", "centrifuge_report"))
all_outputs.extend(centrifuge)

rule download_centrifuge_database:
    output:
        db=[config["dbdir"]+"/centrifuge/p+h+v.{n}.cf".format(n=num) for num in (1,2,3)]
    shadow:
        "shallow"
    params:
        dbdir=config["dbdir"]+"/centrifuge"
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
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        classifications=config["outdir"]+"/centrifuge/{sample}.centrifuge.tsv",
        report=config["outdir"]+"/centrifuge/{sample}.centrifuge_report.tsv",
    log:
        config["outdir"]+"/logs/centrifuge/{sample}.centrifuge.log"
    shadow: 
        "shallow"
    conda:
        "../../envs/centrifuge.yaml"
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
            2> {log}
        """


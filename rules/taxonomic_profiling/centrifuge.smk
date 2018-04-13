# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Centrifuge

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


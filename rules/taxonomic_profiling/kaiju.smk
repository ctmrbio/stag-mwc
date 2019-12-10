# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    download_kaiju_database,
    create_kaiju_krona_plot,

kaiju_config = config["kaiju"]
if config["taxonomic_profile"]["kaiju"]:
    if not all([Path(kaiju_config["db"]).exists(), 
                Path(kaiju_config["nodes"]).exists(),
                Path(kaiju_config["names"]).exists()]):
        err_message = "No Kaiju database files at: '{}', '{}', '{}'!\n".format(kaiju_config["db"], kaiju_config["nodes"], kaiju_config["names"])
        err_message += "Specify relevant paths in the kaiju section of config.yaml.\n"
        err_message += "Run 'snakemake download_kaiju_database' to download a copy into '{dbdir}'\n".format(dbdir=DBDIR/"kaiju") 
        err_message += "If you do not want to run Kaiju for taxonomic profiling, set 'kaiju: False' in config.yaml"
        raise WorkflowError(err_message)

    # Add Kaiju output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    kaiju = expand(str(OUTDIR/"kaiju/{sample}.kaiju"), sample=SAMPLES)
    kaiju_reports = expand(str(OUTDIR/"kaiju/{sample}.kaiju.summary.species"), sample=SAMPLES)
    kaiju_krona = str(OUTDIR/"kaiju/all_samples.kaiju.krona.html")
    all_outputs.extend(kaiju)
    all_outputs.extend(kaiju_reports)
    all_outputs.append(kaiju_krona)

    citations.add((
        "Menzel, P., Ng, K. L., & Krogh, A. (2016).",
        "Fast and sensitive taxonomic classification for metagenomics with Kaiju.",
        "Nature communications, 7, 11257.",
        "Available online at: https://github.com/bioinformatics-centre/kaiju",
    ))
    citations.add((
        "Ondov BD, Bergman NH, and Phillippy AM.",
        "Interactive metagenomic visualization in a Web browser.",
        "BMC Bioinformatics. 2011 Sep 30; 12(1):385.",
    ))


rule download_kaiju_database:
    output:
        db=DBDIR/"kaiju/kaiju_db.fmi",
        names=DBDIR/"kaiju/names.dmp",
        nodes=DBDIR/"kaiju/nodes.dmp"
    shadow:
        "shallow"
    params:
        dbdir=DBDIR/"kaiju"
    shell:
        """
        wget http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz \
        && \
        tar -xf kaiju_index_pg.tgz \
        && \
        mv kaiju_db.fmi names.dmp nodes.dmp {params.dbdir}
        """

rule kaiju:
    input:
        read1=OUTDIR/"host_removal/{sample}_R1.host_removal.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_R2.host_removal.fq.gz",
    output:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    log:
        str(LOGDIR/"kaiju/{sample}.kaiju.log")
    shadow: 
        "shallow"
    threads:
        4
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        db=kaiju_config["db"],
        nodes=kaiju_config["nodes"],
    shell:
        """
        if [[ "{input.read1}" == *.gz ]]
        then
            kaiju \
                -z {threads} \
                -t {params.nodes} \
                -f {params.db} \
                -i <(gunzip -c {input.read1}) \
                -j <(gunzip -c {input.read2}) \
                -o {output.kaiju} > {log}
        else 
            kaiju \
                -z {threads} \
                -t {params.nodes} \
                -f {params.db} \
                -i {input.read1} \
                -j {input.read2} \
                -o {output.kaiju} > {log}
        fi
        """


rule kaiju_report:
    input:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    output:
        krona=OUTDIR/"kaiju/{sample}.krona",
        family=OUTDIR/"kaiju/{sample}.summary.family",
        genus=OUTDIR/"kaiju/{sample}.kaiju.summary.genus",
        species=OUTDIR/"kaiju/{sample}.kaiju.summary.species",
    shadow: 
        "shallow"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
		kaiju2krona \
			-t {params.nodes} \
			-n {params.names} \
			-i {input.kaiju} \
			-o {output.krona} \
            -u
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r species \
            -l superkingdom,phylum,class,order,family,genus,species \
            -o {output.species}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r genus \
            -l superkingdom,phylum,class,order,family,genus,species \
            -o {output.genus}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r family \
            -l superkingdom,phylum,class,order,family,genus,species \
            -o {output.family}
        """

rule create_kaiju_krona_plot:
    input:
        expand(str(OUTDIR/"kaiju/{sample}.krona"), sample=SAMPLES)
    output:
        krona_html=OUTDIR/"kaiju/all_samples.kaiju.krona.html",
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
		ktImportText \
			-o {output.krona_html} \
			{input}
        """

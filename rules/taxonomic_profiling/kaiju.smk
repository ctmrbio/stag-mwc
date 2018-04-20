# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju
from snakemake.exceptions import WorkflowError
import os.path

kaiju_config = config["kaiju"]
if not all([os.path.isfile(kaiju_config["db"]), 
            os.path.isfile(kaiju_config["nodes"]),
            os.path.isfile(kaiju_config["names"])]):
    err_message = "No Kaiju database files at: '{}', '{}', '{}'!\n".format(kaiju_config["db"], kaiju_config["nodes"], kaiju_config["names"])
    err_message += "Specify relevant paths in the kaiju section of config.yaml.\n"
    err_message += "Run 'snakemake download_kaiju_database' to download a copy into '{dbdir}/kaiju'\n".format(dbdir=config["dbdir"]) 
    err_message += "If you do not want to run Kaiju for taxonomic profiling, set 'kaiju: False' in config.yaml"
    raise WorkflowError(err_message)

# Add Kaiju output files to 'all_outputs' from the main Snakefile scope.
# SAMPLES is also from the main Snakefile scope.
kaiju = expand("{outdir}/kaiju/{sample}.kaiju", outdir=outdir, sample=SAMPLES)
kaiju_reports = expand("{outdir}/kaiju/{sample}.kaiju.summary.species", outdir=outdir, sample=SAMPLES)
kaiju_krona = expand("{outdir}/kaiju/all_samples.kaiju.krona.html", outdir=outdir)
all_outputs.extend(kaiju)
all_outputs.extend(kaiju_reports)
all_outputs.extend(kaiju_krona)

rule download_kaiju_database:
    output:
        db=config["dbdir"]+"/kaiju/kaiju_db.fmi",
        names=config["dbdir"]+"/kaiju/names.dmp",
        nodes=config["dbdir"]+"/kaiju/nodes.dmp"
    shadow:
        "shallow"
    params:
        dbdir=config["dbdir"]+"/kaiju"
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
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        kaiju=config["outdir"]+"/kaiju/{sample}.kaiju",
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
                -o {output.kaiju} 
        else 
            kaiju \
                -z {threads} \
                -t {params.nodes} \
                -f {params.db} \
                -i {input.read1} \
                -j {input.read2} \
                -o {output.kaiju}
        fi
        """


rule kaiju_report:
    input:
        kaiju=config["outdir"]+"/kaiju/{sample}.kaiju",
    output:
        krona=config["outdir"]+"/kaiju/{sample}.krona",
        family=config["outdir"]+"/kaiju/{sample}.summary.family",
        genus=config["outdir"]+"/kaiju/{sample}.kaiju.summary.genus",
        species=config["outdir"]+"/kaiju/{sample}.kaiju.summary.species",
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
            -p \
            -o {output.species}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r genus \
            -p \
            -o {output.genus}
        kaijuReport \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -r family \
            -p \
            -o {output.family}
        """

rule create_kaiju_krona_plot:
    input:
        expand(config["outdir"]+"/kaiju/{sample}.krona", sample=SAMPLES)
    output:
        krona_html=config["outdir"]+"/kaiju/all_samples.kaiju.krona.html",
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

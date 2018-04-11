# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju

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


kaiju_config = config["kaiju"]
rule kaiju:
    input:
        read1=config["outdir"]+"/filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=config["outdir"]+"/filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        kaiju=config["outdir"]+"/kaiju/{sample}.kaiju",
    shadow: 
        "shallow"
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
        krona_html=config["outdir"]+"/kaiju/{sample}.krona.html",
        family=config["outdir"]+"/kaiju/{sample}.summary.family",
        genus=config["outdir"]+"/kaiju/{sample}.kaiju.summary.genus",
        species=config["outdir"]+"/kaiju/{sample}.kaiju.summary.species",
    shadow: 
        "shallow"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    shell:
        """
		kaiju2krona \
			-t {params.nodes} \
			-n {params.names} \
			-i {input.kaiju} \
			-o {output.krona} \
            -u
		ktImportText \
			-o {output.krona_html} \
			{output.krona}
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


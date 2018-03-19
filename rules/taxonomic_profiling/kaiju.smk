# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju

rule download_kaiju_databases:
    output:
        db=config["dbdir"]+"/kaiju/kaiju_db.fmi",
        names=config["dbdir"]+"/kaiju/names.dmp",
        nodes=config["dbdir"]+"/kaiju/nodes.dmp"
    shadow:
        "shallow"
    shell:
        """
        wget http://kaiju.binf.ku.dk/database/kaiju_index_pg.tgz \
        && \
        tar -xf kaiju_index_pg.tgz \
        && \
        mv kaiju_db.fmi names.dmp nodes.dmp {config[dbdir]}/kaiju
        """

rule kaiju:
    input:
        read1=config["outdir"]+"/trimmed_qa/{sample}_R1.trimmed_qa.fq.gz",
        read2=config["outdir"]+"/trimmed_qa/{sample}_R2.trimmed_qa.fq.gz",
        nodes=config["dbdir"]+"/kaiju/nodes.dmp",
        db=config["dbdir"]+"/kaiju/kaiju_db.fmi",
    output:
        kaiju=config["outdir"]+"/kaiju/{sample}.kaiju",
    shadow: 
        "shallow"
    shell:
        """
        if [[ "{input.read1}" == *.gz ]]
        then
            kaiju \
                -z {threads} \
                -t {input.nodes} \
                -f {input.db} \
                -i <(gunzip -c {input.read1}) \
                -j <(gunzip -c {input.read2}) \
                -o {wildcards.sample}.kaiju
        else
            kaiju \
                -z {threads} \
                -t {input.nodes} \
                -f {input.db} \
                -i {input.read1} \
                -j {input.read2} \
                -o {wildcards.sample}.kaiju
        fi
        """


rule kaiju_report:
    input:
        kaiju=config["outdir"]+"/kaiju/{sample}.kaiju",
        nodes=config["dbdir"]+"/kaiju/nodes.dmp",
        names=config["dbdir"]+"/kaiju/names.dmp",
        db=config["dbdir"]+"/kaiju/kaiju_db.fmi",
    output:
        krona=config["outdir"]+"/kaiju/{sample}.krona",
        krona_html=config["outdir"]+"/kaiju/{sample}.krona.html",
        family=config["outdir"]+"/kaiju/{sample}.summary.family",
        genus=config["outdir"]+"/kaiju/{sample}.kaiju.summary.genus",
        species=config["outdir"]+"/kaiju/{sample}.kaiju.summary.species",
    shadow: 
        "shallow"
    shell:
        """
		kaiju2krona \
			-t {input.nodes} \
			-f {input.db} \
			-i {input.kaiju} \
			-o {output.krona}
		ktImportText \
			-o {output.krona_html} \
			{output.krona}
		for rank in family genus species; do
			kaijuReport \
				-t {input.nodes} \
                -n {input.names} \
				-f {input.db} \
				-i {input.kaiju} \
				-r {{$rank}} \
				-o {config[outdir]}/kaiju/{wildcards.sample}.kaiju.summary.{{$rank}}
		done
        """
    

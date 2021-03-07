# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    download_kaiju_database,
    create_kaiju_krona_plot,
    kaiju2krona,
    kaiju_report,
    join_kaiju_reports,
    kaiju_area_plot,

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
    kaiju_krona = str(OUTDIR/"kaiju/all_samples.kaiju.krona.html")
    kaiju_reports = expand(str(OUTDIR/"kaiju/{sample}.kaiju.{level}.txt"), sample=SAMPLES, level=kaiju_config["levels"])
    kaiju_joined_table = expand(str(OUTDIR/"kaiju/all_samples.kaiju.{level}.txt"), level=kaiju_config["levels"])
    kaiju_area_plot = expand(str(OUTDIR/"kaiju/area_plot.kaiju.pdf"))
    all_outputs.extend(kaiju)
    all_outputs.extend(kaiju_reports)
    all_outputs.append(kaiju_krona)
    all_outputs.append(kaiju_joined_table)
    #all_outputs.append(kaiju_area_plot)  # Buggy in stag 4.1

    citations.add(publications["Kaiju"])
    citations.add(publications["Krona"])


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
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    log:
        str(LOGDIR/"kaiju/{sample}.kaiju.log")
    shadow: 
        "shallow"
    threads:
        cluster_config["kaiju"]["n"] if "kaiju" in cluster_config else 4
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
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


rule kaiju2krona:
    """Convert Kaiju output to Krona input"""
    input:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    output:
        krona=OUTDIR/"kaiju/{sample}.krona",
    shadow: 
        "shallow"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        kaiju2krona \
            -t {params.nodes} \
            -n {params.names} \
            -i {input.kaiju} \
            -o {output.krona} \
            -u
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
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        ktImportText \
            -o {output.krona_html} \
            {input}
        """

rule kaiju_report:
    input:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    output:
        OUTDIR/"kaiju/{sample}.kaiju.{level}.txt",
    log:
        str(LOGDIR/"kaiju/kaiju2table.{sample}.{level}.log")
    shadow: 
        "shallow"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        kaiju2table \
            -t {params.nodes} \
            -n {params.names} \
            -r {wildcards.level} \
            -l superkingdom,phylum,class,order,family,genus,species \
            -o {output} \
            {input.kaiju} \
            2>&1 > {log}
        """


rule join_kaiju_reports:
    input:
        expand(str(OUTDIR/"kaiju/{sample}.kaiju.{{level}}.txt"), sample=SAMPLES),
    output:
        report(OUTDIR/"kaiju/all_samples.kaiju.{level}.txt",
            category="Taxonomic profiling",
            caption="../../report/kaiju_table.rst")
    log:
        str(LOGDIR/"kaiju/join_kaiju_reports.{level}.log")
    shadow: 
        "shallow"
    params:
        feature_column=kaiju_config["feature_column"],
        value_column=kaiju_config["value_column"],
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/join_tables.py \
            --feature-column {params.feature_column} \
            --value-column {params.value_column} \
            --outfile {output} \
            {input} \
            2>&1 > {log}
        """
    

rule kaiju_area_plot:
    input:
        OUTDIR/"kaiju/all_samples.kaiju.{level}.txt".format(level=kaiju_config["levels"][-1])
    output:
        report(OUTDIR/"kaiju/area_plot.kaiju.pdf",
            category="Taxonomic profiling",
            caption="../../report/area_plot.rst")
    log:
        str(LOGDIR/"kaiju/area_plot.log")
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode kaiju \
            2>&1 > {log}
        """


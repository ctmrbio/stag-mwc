# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kaiju
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
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
    kaiju = expand(OUTDIR/"kaiju/{sample}.kaiju", sample=SAMPLES)
    kaiju_krona = OUTDIR/"kaiju/all_samples.kaiju.krona.html"
    kaiju_reports = expand(OUTDIR/"kaiju/{sample}.kaiju.report_{level}.txt", sample=SAMPLES, level=kaiju_config["levels"])
    kaiju_joined_table = expand(OUTDIR/"kaiju/all_samples.kaiju.{level}.txt", level=kaiju_config["levels"])
    kaiju_area_plot = expand(OUTDIR/"kaiju/area_plot.kaiju.pdf")

    all_outputs.extend(kaiju)
    all_outputs.extend(kaiju_reports)
    all_outputs.append(kaiju_joined_table)
    #all_outputs.append(kaiju_area_plot)  # Buggy in stag 4.1

    citations.add(publications["Kaiju"])

    if kaiju_config["run_krona"]:
        all_outputs.append(kaiju_krona)
        citations.add(publications["Krona"])


rule kaiju:
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        kaiju=OUTDIR/"kaiju/{sample}.kaiju",
    log:
        LOGDIR/"kaiju/{sample}.kaiju.log",
    shadow: 
        "shallow"
    threads: 8
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    params:
        db=kaiju_config["db"],
        nodes=kaiju_config["nodes"],
    shell:
        """
        kaiju \
            -z {threads} \
            -t {params.nodes} \
            -f {params.db} \
            -i {input.read1} \
            -j {input.read2} \
            -o {output.kaiju} > {log}
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
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
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
        expand(OUTDIR/"kaiju/{sample}.krona", sample=SAMPLES)
    output:
        krona_html=report(OUTDIR/"kaiju/all_samples.kaiju.krona.html",
            category="Taxonomic profiling",
            caption="../../report/kaiju_krona.rst")
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
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
        OUTDIR/"kaiju/{sample}.kaiju.report_{level}.txt",
    log:
        LOGDIR/"kaiju/kaiju2table.{sample}.{level}.log",
    shadow: 
        "shallow"
    params:
        nodes=kaiju_config["nodes"],
        names=kaiju_config["names"],
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
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
        expand(OUTDIR/"kaiju/{sample}.kaiju.report_{{level}}.txt", sample=SAMPLES),
    output:
        report(OUTDIR/"kaiju/all_samples.kaiju.{level}.txt",
            category="Taxonomic profiling",
            caption="../../report/kaiju_table.rst")
    log:
        LOGDIR/"kaiju/join_kaiju_reports.{level}.log",
    shadow: 
        "shallow"
    params:
        feature_column=kaiju_config["feature_column"],
        value_column=kaiju_config["value_column"],
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        workflow/scripts/join_tables.py \
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
        LOGDIR/"kaiju/area_plot.log",
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        workflow/scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode kaiju \
            2>&1 > {log}
        """


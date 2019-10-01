# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kraken2
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    combine_kreports,
    create_kraken2_krona_plot,
    download_KrakenTools,
    download_minikraken2,
    kreport2krona,

kraken2_config = config["kraken2"]
if config["taxonomic_profile"]["kraken2"]:
    if not Path(kaiju_config["db"]).exists():
        err_message = "No Kraken2 database folder at: '{}'!\n".format(kraken2_config["db"])
        err_message += "Specify the path in the kraken2 section of config.yaml.\n"
        err_message += "Run 'snakemake download_minikraken2' to download a copy into '{dbdir}'\n".format(dbdir=DBDIR/"kraken2") 
        err_message += "If you do not want to run kraken2 for taxonomic profiling, set 'kraken2: False' in config.yaml"
        raise WorkflowError(err_message)

    # Add Kraken2 output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    krakens = expand(str(OUTDIR/"kraken2/{sample}.kraken"), sample=SAMPLES)
    kreports = expand(str(OUTDIR/"kraken2/{sample}.kreport"), sample=SAMPLES)
    combined_kreport = str(OUTDIR/"kraken2/all_samples.kraken2.tsv")
    kraken_krona = str(OUTDIR/"kraken2/all_samples.kraken2.krona.html")
    all_outputs.extend(krakens)
    all_outputs.extend(kreports)
    all_outputs.append(combined_kreport)
    all_outputs.append(kraken_krona)
    
    citations.add((
        "Wood DE, Salzberg SL (2014).",
        "Kraken: ultrafast metagenomic sequence classification using exact alignments.",
        "Genome Biology 2014, 15:R46.",
    ))
    citations.add((
        "Ondov BD, Bergman NH, and Phillippy AM.",
        "Interactive metagenomic visualization in a Web browser.",
        "BMC Bioinformatics. 2011 Sep 30; 12(1):385.",
    ))


rule download_minikraken2:
    output:
        db=DBDIR/"kraken2/minikraken2_v1_8GB/hash.k2d",
        names=DBDIR/"kraken2/minikraken2_v1_8GB/opts.k2d",
        nodes=DBDIR/"kraken2/minikraken2_v1_8GB/taxo.k2d"
    log:
        str(LOGDIR/"kraken2/download_minikraken2.log")
    shadow:
        "shallow"
    params:
        dbdir=DBDIR/"kraken2/minikraken2_v1_8GB"
    shell:
        """
        wget https://ccb.jhu.edu/software/kraken2/dl/minikraken2_v1_8GB.tgz > {log}
        tar -vxf minikraken2_v1_8GB.tgz  >> {log}
        mv -v *k2d {params.dbdir} >> {log}
        """


rule download_KrakenTools:
    """Download kreport2krona.py from Jennifer Lu's github."""
    output:
        "scripts/KrakenTools/kreport2krona.py"
    log:
        str(LOGDIR/"kraken2/KrakenTools.log")
    shell:
        """
        git clone \
            https://github.com/jenniferlu717/KrakenTools/tree/4db3e990db4e949b58b17a4d8bb133504818bc93 \
            scripts/KrakenTools \
            2> {log}
        """


rule kraken2:
    input:
        read1=OUTDIR/"filtered_human/{sample}_R1.filtered_human.fq.gz",
        read2=OUTDIR/"filtered_human/{sample}_R2.filtered_human.fq.gz",
    output:
        kraken=OUTDIR/"kraken2/{sample}.kraken",
        kreport=OUTDIR/"kraken2/{sample}.kreport"
    log:
        str(LOGDIR/"kraken2/{sample}.kraken2.log")
    shadow: 
        "shallow"
    threads:
        4
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        db=kraken2_config["db"],
        confidence=kraken2_config["confidence"],
        extra=kraken2_config["extra"],
    shell:
        """
        kraken2 \
            --db {params.db} \
            --confidence {params.confidence} \
            --threads {threads} \
            --output {output.kraken} \
            --report {output.kreport} \
            --use-names \
            --paired \
            {input.read1} {input.read2} \
            {params.extra} \
            2> {log}
        """
        

rule combine_kreports:
    input:
        kreports=expand(str(OUTDIR/"kraken2/{sample}.kreport"), sample=SAMPLES),
        combine="scripts/KrakenTools/combine_kreports.py",
    output:
        report(OUTDIR/"kraken2/all_samples.kraken2.tsv",
               category="Taxonomic profiling",
               caption="../../report/kraken2_table.rst"),
    log:
        str(LOGDIR/"kraken2/combined_kreport.log")
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
        {input.combine} \
            --output {output} \
            --report-files {input.kreports} \
            2>> {log} \
            >> {log}
        """


rule kreport2krona:
    input:
        kreport=OUTDIR/"kraken2/{sample}.kreport",
        kreport2krona="scripts/KrakenTools/kreport2krona.py",
    output:
        OUTDIR/"kraken2/{sample}.krona"
    log:
        str(LOGDIR/"kraken2/{sample}.kreport2krona.log")
    shadow: 
        "shallow"
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    shell:
        """
        {input.kreport2krona} \
            --report-file {input.kreport} \
            --output {output} \
            2> {log}
        """


rule create_kraken2_krona_plot:
    input:
        expand(str(OUTDIR/"kraken2/{sample}.krona"), sample=SAMPLES),
    output:
        krona_html=report(OUTDIR/"kraken2/all_samples.kraken2.krona.html",
                          category="Taxonomic profiling",
                          caption="../../report/kraken2_krona.rst"),
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

# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kraken2
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    download_minikraken2,
    create_kraken2_krona_plot,

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
    kraken_krona = str(OUTDIR/"kraken2/all_samples.kraken.krona.html")
    all_outputs.extend(krakens)
    all_outputs.extend(kreports)
    #all_outputs.append(kraken_krona) # TODO: add krona output


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
            --report {output.report} \
            --use-names \
            --paired \
            {input.read1} {input.read2} \
            {params.extra}
        """


rule create_kraken2_krona_plot:
    input:
        expand(str(OUTDIR/"kraken2/{sample}.krona"), sample=SAMPLES)
    output:
        krona_html=OUTDIR/"kraken2/all_samples.kraken2.krona.html",
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

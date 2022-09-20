# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using KrakenUniq
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    combine_krakenuniq_reports


krakenuniq_config = config["krakenuniq"]
if config["taxonomic_profile"]["krakenuniq"]:
    if not (krakenuniq_config["db"] and Path(krakenuniq_config["db"]).exists()):
        err_message = "No KrakenUniq database folder at: '{}'!\n".format(krakenuniq_config["db"])
        err_message += "Specify the path in the krakenuniq section of config.yaml.\n"
        err_message += "If you do not want to run krakenuniq for taxonomic profiling, set 'krakenuniq: False' in config.yaml"
        raise WorkflowError(err_message)

    # Add KrakenUniq output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    krakens = expand(OUTDIR/"krakenuniq/{sample}.kraken.gz", sample=SAMPLES)
    kreports = expand(OUTDIR/"krakenuniq/{sample}.kreport", sample=SAMPLES)
    combined_kreport = expand(OUTDIR/"krakenuniq/all_samples.krakenuniq.txt", sample=SAMPLES)
    all_outputs.extend(krakens)
    all_outputs.extend(kreports)
    all_outputs.append(combined_kreport)
    
    citations.add(publications["KrakenUniq"])
    citations.add(publications["Krona"])


rule krakenuniq:
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        kraken=OUTDIR/"krakenuniq/{sample}.kraken.gz" if krakenuniq_config["keep_kraken"] else temp(OUTDIR/"krakenuniq/{sample}.kraken.gz"),
        kreport=OUTDIR/"krakenuniq/{sample}.kreport" if krakenuniq_config["keep_kreport"] else temp(OUTDIR/"krakenuniq/{sample}.kreport"),
    log:
        LOGDIR/"krakenuniq/{sample}.krakenuniq.log"
    shadow: 
        "shallow"
    threads:
        cluster_config["krakenuniq"]["n"] if "krakenuniq" in cluster_config else 4
    conda:
        "../../envs/krakenuniq.yaml"
    container:
        "docker://quay.io/biocontainers/krakenuniq:0.7.3--pl5321h19e8d03_0"
    params:
        db=krakenuniq_config["db"],
        extra=krakenuniq_config["extra"],
    shell:
        """
        krakenuniq \
            --db {params.db} \
            --threads {threads} \
            --output {output.kraken} \
            --report-file {output.kreport} \
            --paired \
            {input.read1} {input.read2} \
            {params.extra} \
            2> {log}
        """


rule combine_krakenuniq_reports:
    input:
        kreports=expand(OUTDIR/"krakenuniq/{sample}.kreport", sample=SAMPLES)
    output:
        combined=OUTDIR/"krakenuniq/all_samples.krakenuniq.txt"
    log:
        LOGDIR/"krakenuniq/all_samples.krakenuniq.log"
    shadow: 
        "shallow"
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
    shell:
        """
        scripts/join_tables.py \
            --feature-column rank,taxName \
            --value-column taxReads \
            --outfile {output.combined} \
            --skiplines 3 \
            {input.kreports} \
            2> {log}
        """

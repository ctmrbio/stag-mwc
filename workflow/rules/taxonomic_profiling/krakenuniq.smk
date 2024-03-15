# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using KrakenUniq
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    krakenuniq_combine_reports,
    krakenuniq_mpa_style,
    krakenuniq_join_mpa,
    krakenuniq_krona_plot,
    krakenuniq_kreport2krona,


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
    mpa_table = OUTDIR/"krakenuniq/all_samples.krakenuniq.mpa_style.txt"
    krona_html = OUTDIR/"krakenuniq/all_samples.krakenuniq.krona.html"

    if krakenuniq_config["keep_kraken"]:
        all_outputs.extend(krakens)
    if krakenuniq_config["keep_kreport"]:
        all_outputs.extend(kreports)
    if krakenuniq_config["run_krona"]:
        all_outputs.append(krona_html)
    all_outputs.append(combined_kreport)
    all_outputs.append(mpa_table)
    
    citations.add(publications["KrakenUniq"])
    citations.add(publications["Krona"])

rule krakenuniq_merge_reads:
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        fasta=temp(OUTDIR/"krakenuniq/{sample}.tmp.fa.gz"),
    log:
        LOGDIR/"krakenuniq/{sample}.read_merge.log"
    shadow:
        "shallow"
    threads: 4
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["bbmap"]
    shell:
        """
        fuse.sh \
            in1={input.read1} \
            in2={input.read2} \
            out={output.fasta} \
            pad=1 \
            fusepairs=t \
            2> {log}
        """

rule krakenuniq:
    input:
        fasta=OUTDIR/"krakenuniq/{sample}.tmp.fa.gz",
    output:
        kraken=OUTDIR/"krakenuniq/{sample}.kraken.gz" if krakenuniq_config["keep_kraken"] else temp(OUTDIR/"krakenuniq/{sample}.kraken.gz"),
        kreport=OUTDIR/"krakenuniq/{sample}.kreport" if krakenuniq_config["keep_kreport"] else temp(OUTDIR/"krakenuniq/{sample}.kreport"),
    log:
        LOGDIR/"krakenuniq/{sample}.krakenuniq.log"
    shadow:
        "shallow"
    threads: 8
    conda:
        "../../envs/krakenuniq.yaml"
    container:
        config["containers"]["krakenuniq"]
    params:
        db=krakenuniq_config["db"],
        preload_size=krakenuniq_config["preload_size"],
        extra=krakenuniq_config["extra"],
    shell:
        """
        krakenuniq \
            --db {params.db} \
            --threads {threads} \
            --output {output.kraken} \
            --report-file {output.kreport} \
            --preload-size {params.preload_size} \
            {input.fasta} \
            {params.extra} \
            2> {log}
        """


rule krakenuniq_combine_reports:
    input:
        kreports=expand(OUTDIR/"krakenuniq/{sample}.kreport", sample=SAMPLES)
    output:
        combined=report(OUTDIR/"krakenuniq/all_samples.krakenuniq.txt",
                category="Taxonomic profiling",
                caption="../../report/krakenuniq_table.rst"),
    log:
        LOGDIR/"krakenuniq/all_samples.krakenuniq.log"
    shadow: 
        "shallow"
    threads: 1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["stag"]
    shell:
        """
        workflow/scripts/join_tables.py \
            --feature-column rank,taxName \
            --value-column taxReads \
            --outfile {output.combined} \
            --skiplines 2 \
            {input.kreports} \
            2> {log}
        """


rule krakenuniq_mpa_style:
    input:
        kreport=OUTDIR/"krakenuniq/{sample}.kreport"
    output:
        txt=OUTDIR/"krakenuniq/{sample}.mpa_style.txt",
    log:
        stdout=LOGDIR/"krakenuniq/{sample}.mpa_style.stdout",
        stderr=LOGDIR/"krakenuniq/{sample}.mpa_style.stderr",
    threads: 1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["stag"]
    shell:
        """
        workflow/scripts/KrakenTools/kreport2mpa.py \
            --report-file {input.kreport} \
            --output {output.txt} \
            --display-header \
            > {log.stdout} \
            2> {log.stderr}
        """


rule krakenuniq_join_mpa:
    input:
        txt=expand(OUTDIR/"krakenuniq/{sample}.mpa_style.txt", sample=SAMPLES),
    output:
        table=report(OUTDIR/"krakenuniq/all_samples.krakenuniq.mpa_style.txt",
               category="Taxonomic profiling",
               caption="../../report/krakenuniq_table_mpa.rst"),
    log:
        stdout=LOGDIR/"krakenuniq/join_krakenuniq_mpa_tables.stdout",
        stderr=LOGDIR/"krakenuniq/join_krakenuniq_mpa_tables.stderr",
    threads: 1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["stag"]
    params:
        value_column="reads",
        feature_column="taxon_name",
    shell:
        """
        workflow/scripts/join_tables.py \
            --outfile {output.table} \
            --value-column {params.value_column} \
            --feature-column '{params.feature_column}' \
            {input.txt} \
            > {log.stdout} \
            2> {log.stderr}
        """


rule krakenuniq_kreport2krona:
    input:
        kreport=OUTDIR/"krakenuniq/{sample}.mpa_style.txt",
    output:
        OUTDIR/"krakenuniq/{sample}.krona"
    log:
        stdout=LOGDIR/"krakenuniq/{sample}.kreport2krona.stdout",
        stderr=LOGDIR/"krakenuniq/{sample}.kreport2krona.stderr",
    shadow:
        "shallow"
    threads: 1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["stag"]
    shell:
        """
        awk -v OFS='\\t' '{{
          gsub("\\\\|","\\t",$1);
          print $2,$1;
          }}' {input.kreport} \
          > {output} \
          2> {log.stderr}
        """


rule krakenuniq_krona_plot:
    input:
        expand(OUTDIR/"krakenuniq/{sample}.krona", sample=SAMPLES),
    output:
        krona_html=report(OUTDIR/"krakenuniq/all_samples.krakenuniq.krona.html",
                          category="Taxonomic profiling",
                          caption="../../report/krakenuniq_krona.rst"),
    shadow:
        "shallow"
    threads: 1
    conda:
        "../../envs/stag-mwc.yaml"
    container:
        config["containers"]["krona"]
    shell:
        """
		ktImportText \
			-o {output.krona_html} \
			{input}
        """

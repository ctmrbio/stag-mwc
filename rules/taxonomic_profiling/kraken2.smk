# vim: syntax=python expandtab
# Taxonomic classification of metagenomic reads using Kraken2 with abundance
# estimation using Bracken
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path

from snakemake.exceptions import WorkflowError

localrules:
    bracken_mpa_style,
    bracken_kreport,
    bracken_all_levels,
    combine_kreports,
    create_kraken2_krona_plot,
    download_minikraken2,
    filter_bracken,
    join_bracken,
    join_bracken_filtered,
    join_bracken_mpa,
    join_kraken2_mpa,
    kreport2krona,
    bracken2krona,
    create_bracken_krona_plot,
    kraken2_area_plot,

kraken2_config = config["kraken2"]
if config["taxonomic_profile"]["kraken2"]:
    if not Path(kraken2_config["db"]).exists():
        err_message = "No Kraken2 database folder at: '{}'!\n".format(kraken2_config["db"])
        err_message += "Specify the path in the kraken2 section of config.yaml.\n"
        err_message += "Run 'snakemake download_minikraken2' to download a copy into '{dbdir}'\n".format(dbdir=DBDIR/"kraken2") 
        err_message += "If you do not want to run kraken2 for taxonomic profiling, set 'kraken2: False' in config.yaml"
        raise WorkflowError(err_message)

    # Add Kraken2 output files to 'all_outputs' from the main Snakefile scope.
    # SAMPLES is also from the main Snakefile scope.
    krakens = expand(str(OUTDIR/"kraken2/{sample}.kraken"), sample=SAMPLES)
    kreports = expand(str(OUTDIR/"kraken2/{sample}.kreport"), sample=SAMPLES)
    kreports_mpa_style = expand(str(OUTDIR/"kraken2/{sample}.mpa_style.txt"), sample=SAMPLES)
    joined_kreport_mpa_style = str(OUTDIR/"kraken2/all_samples.kraken2.mpa_style.txt")
    combined_kreport = str(OUTDIR/"kraken2/all_samples.kraken2.txt")
    kraken_krona = str(OUTDIR/"kraken2/all_samples.kraken2.krona.html")
    kraken_area_plot = str(OUTDIR/"kraken2/area_plot.kraken2.pdf")
    all_outputs.extend(krakens)
    all_outputs.extend(kreports)
    all_outputs.extend(kreports_mpa_style)
    all_outputs.append(joined_kreport_mpa_style)
    all_outputs.append(combined_kreport)
    all_outputs.append(kraken_krona)
    all_outputs.append(kraken_area_plot)
    
    citations.add(publications["Kraken2"])
    citations.add(publications["Krona"])


rule download_minikraken2:
    output:
        db=DBDIR/"kraken2/minikraken2_v2_8GB_201904_UPDATE/hash.k2d",
        names=DBDIR/"kraken2/minikraken2_v2_8GB_201904_UPDATE/opts.k2d",
        nodes=DBDIR/"kraken2/minikraken2_v2_8GB_201904_UPDATE/taxo.k2d",
        kmer_distrib=DBDIR/"kraken2/minikraken2_v2_8GB_201904_UPDATE/database150mers.kmer_distrib",
    log:
        str(LOGDIR/"kraken2/download_minikraken2.log")
    shadow:
        "shallow"
    params:
        dbdir=DBDIR/"kraken2/minikraken2_v2_8GB_201904_UPDATE"
    shell:
        """
        wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz > {log}
        tar -vxf minikraken2_v2_8GB_201904_UPDATE.tgz  >> {log}
        mv -v *k2d *kmer_distrib {params.dbdir} >> {log}
        """


rule kraken2:
    input:
        read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
        read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
    output:
        kraken=OUTDIR/"kraken2/{sample}.kraken",
        kreport=OUTDIR/"kraken2/{sample}.kreport",
    log:
        str(LOGDIR/"kraken2/{sample}.kraken2.log")
    shadow: 
        "shallow"
    threads:
        cluster_config["kraken2"]["n"] if "kraken2" in cluster_config else 4
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
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
        

rule kraken_mpa_style:
    input:
        kreport=OUTDIR/"kraken2/{sample}.kreport"
    output:
        txt=OUTDIR/"kraken2/{sample}.mpa_style.txt",
    log:
        str(LOGDIR/"kraken2/{sample}.mpa_style.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        kreport2mpa.py \
            --report-file {input.kreport} \
            --output {output.txt} \
            --display-header \
            2>&1 > {log}
        sed --in-place 's|{input.kreport}|taxon_name\treads|g' {output.txt}
        """


rule join_kraken2_mpa:
    input:
        txt=expand(str(OUTDIR/"kraken2/{sample}.mpa_style.txt"), sample=SAMPLES),
    output:
        table=report(OUTDIR/"kraken2/all_samples.kraken2.mpa_style.txt",
               category="Taxonomic profiling",
               caption="../../report/kraken2_table_mpa.rst"),
    log:
        str(LOGDIR/"kraken2/join_kraken2_mpa_tables.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        value_column="reads",
        feature_column="taxon_name",
    shell:
        """
        scripts/join_tables.py \
            --outfile {output.table} \
            --value-column {params.value_column} \
            --feature-column {params.feature_column} \
            {input.txt} \
            2>&1 > {log}
        """


rule kraken2_area_plot:
    input:
        OUTDIR/"kraken2/all_samples.kraken2.mpa_style.txt"
    output:
        report(OUTDIR/"kraken2/area_plot.kraken2.pdf",
            category="Taxonomic profiling",
            caption="../../report/area_plot.rst")
    log:
        str(LOGDIR/"kraken2/area_plot.kraken2.log")
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode kraken2 \
            2>&1 > {log}
        """


rule combine_kreports:
    input:
        kreports=expand(str(OUTDIR/"kraken2/{sample}.kreport"), sample=SAMPLES),
    output:
        report(OUTDIR/"kraken2/all_samples.kraken2.txt",
               category="Taxonomic profiling",
               caption="../../report/kraken2_table.rst"),
    log:
        str(LOGDIR/"kraken2/combined_kreport.log")
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/KrakenTools/combine_kreports.py \
            --output {output} \
            --report-files {input.kreports} \
            2>> {log} \
            >> {log}
        """


rule kreport2krona:
    input:
        kreport=OUTDIR/"kraken2/{sample}.kreport",
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
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/KrakenTools/kreport2krona.py \
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
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
		ktImportText \
			-o {output.krona_html} \
			{input}
        """


if config["taxonomic_profile"]["kraken2"] and kraken2_config["bracken"]["kmer_distrib"]:
    if not Path(kraken2_config["bracken"]["kmer_distrib"]).exists():
        err_message = "No Bracken kmer_distrib database file at: '{}'!\n".format(kraken2_config["bracken"]["kmer_distrib"])
        err_message += "Specify the path in the kraken2 section of config.yaml.\n"
        err_message += "Run 'snakemake download_minikraken2' to download a copy of the required files into '{dbdir}'\n".format(dbdir=DBDIR/"kraken2") 
        err_message += "If you do not want to run Bracken for abundance profiling, set 'kmer_distrib: ""' in the bracken section of config.yaml"
        raise WorkflowError(err_message)
    if kraken2_config["filter_bracken"]["include"] or kraken2_config["filter_bracken"]["exclude"]:
        filtered_brackens = expand(str(OUTDIR/"kraken2/{sample}.{level}.filtered.bracken"), sample=SAMPLES, level=kraken2_config["bracken"]["levels"].split())
        all_table = expand(str(OUTDIR/"kraken2/all_samples.{level}.bracken.txt"), level=kraken2_config["bracken"]["levels"].split())
        all_table_filtered = expand(str(OUTDIR/"kraken2/all_samples.{level}.filtered.bracken.txt"), level=kraken2_config["bracken"]["levels"].split())

        all_outputs.extend(filtered_brackens)
        all_outputs.append(all_table)
        all_outputs.append(all_table_filtered)

    citations.add(publications["Bracken"])

    brackens = expand(str(OUTDIR/"kraken2/{sample}.{level}.bracken"), sample=SAMPLES, level=kraken2_config["bracken"]["levels"].split())
    brackens_mpa_style = expand(str(OUTDIR/"kraken2/{sample}.bracken.mpa_style.txt"), sample=SAMPLES)
    bracken_area_plot = str(OUTDIR/"kraken2/area_plot.bracken.pdf")
    bracken_krona = str(OUTDIR/"kraken2/all_samples.bracken.krona.html")
    all_table_mpa = str(OUTDIR/"kraken2/all_samples.bracken.mpa_style.txt")

    all_outputs.extend(brackens)
    all_outputs.extend(brackens_mpa_style)
    all_outputs.append(bracken_krona)
    all_outputs.append(all_table_mpa)
    all_outputs.append(bracken_area_plot)


rule bracken_kreport:
    """Run Bracken summarization for Species level to get total sample
    bracken.kreport (required for mpa-conversion later)."""
    input:
        kreport=OUTDIR/"kraken2/{sample}.kreport"
    output:
        bracken_kreport=OUTDIR/"kraken2/{sample}_bracken.kreport",
    log:
        str(LOGDIR/"kraken2/{sample}.bracken.log")
    threads:
        cluster_config["bracken"]["n"] if "bracken" in cluster_profile else 2
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        kmer_distrib=kraken2_config["bracken"]["kmer_distrib"],
        thresh=kraken2_config["bracken"]["thresh"],
    shell:
        """
        est_abundance.py \
            --input {input.kreport} \
            --kmer_distr {params.kmer_distrib} \
            --output {output.bracken_kreport} \
            --level S \
            --thresh {params.thresh} \
            2>&1 > {log}
        """


rule bracken_all_levels:
    """Run Bracken summarization for all levels."""
    input:
        kreport=OUTDIR/"kraken2/{sample}.kreport"
    output:
        bracken=OUTDIR/"kraken2/{sample}.{level,[DPOCFGS]}.bracken",
    log:
        str(LOGDIR/"kraken2/{sample}.{level}.bracken.log")
    shadow:         # shadow required because est_abundance.py always creates the
        "shallow"   # sample-level output file with fixed filename: {sample}_bracken.kreport 
    threads:
        cluster_config["bracken"]["n"] if "bracken" in cluster_profile else 2
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        kmer_distrib=kraken2_config["bracken"]["kmer_distrib"],
        thresh=kraken2_config["bracken"]["thresh"],
    shell:
        """
        est_abundance.py \
            --input {input.kreport} \
            --kmer_distr {params.kmer_distrib} \
            --output {output.bracken} \
            --level {wildcards.level} \
            --thresh {params.thresh} \
            2>&1 > {log}
        """


rule bracken_mpa_style:
    input:
        kreport=OUTDIR/"kraken2/{sample}_bracken.kreport"
    output:
        txt=OUTDIR/"kraken2/{sample}.bracken.mpa_style.txt",
    log:
        str(LOGDIR/"kraken2/{sample}.bracken.mpa_style.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        kreport2mpa.py \
            --report-file {input.kreport} \
            --output {output.txt} \
            --display-header \
            2>&1 > {log}
        sed --in-place 's|{input.kreport}|taxon_name\treads|g' {output.txt}
        """


rule join_bracken_mpa:
    input:
        txt=expand(str(OUTDIR/"kraken2/{sample}.bracken.mpa_style.txt"), sample=SAMPLES),
    output:
        table=report(OUTDIR/"kraken2/all_samples.bracken.mpa_style.txt",
               category="Taxonomic profiling",
               caption="../../report/bracken_table_mpa.rst"),
    log:
        str(LOGDIR/"kraken2/join_bracken_mpa_tables.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        value_column="reads",
        feature_column="taxon_name",
    shell:
        """
        scripts/join_tables.py \
            --outfile {output.table} \
            --value-column {params.value_column} \
            --feature-column {params.feature_column} \
            {input.txt} \
            2>&1 > {log}
        """


rule bracken_area_plot:
    input:
        OUTDIR/"kraken2/all_samples.bracken.mpa_style.txt"
    output:
        report(OUTDIR/"kraken2/area_plot.bracken.pdf",
            category="Taxonomic profiling",
            caption="../../report/area_plot.rst")
    log:
        str(LOGDIR/"kraken2/area_plot.bracken.log")
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/area_plot.py \
            --table {input} \
            --output {output} \
            --mode kraken2 \
            2>&1 > {log}
        """


rule join_bracken:
    input:
        bracken=expand(str(OUTDIR/"kraken2/{sample}.{{level}}.bracken"), sample=SAMPLES),
    output:
        table=report(OUTDIR/"kraken2/all_samples.{level,[DPOCFGS]}.bracken.txt",
               category="Taxonomic profiling",
               caption="../../report/bracken_table.rst"),
    log:
        str(LOGDIR/"kraken2/join_bracken_tables.{level}.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        value_column="fraction_total_reads",
        feature_column="name",
    shell:
        """
        scripts/join_tables.py \
            --outfile {output.table} \
            --value-column {params.value_column} \
            --feature-column {params.feature_column} \
            {input.bracken} \
            2>&1 > {log}
        """
    
rule bracken2krona:
    """Convert Bracken kreport output to krona"""
    input:
        bracken_kreport=OUTDIR/"kraken2/{sample}_bracken.kreport",
    output:
        bracken_krona=OUTDIR/"kraken2/{sample}.bracken.krona",
    log:
        str(LOGDIR/"kraken2/{sample}.bracken2krona.log")
    threads:
        1
    shadow:
        "shallow"
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    shell:
        """
        scripts/KrakenTools/kreport2krona.py \
            --report-file {input.bracken_kreport} \
            --output {output.bracken_krona} \
            2>&1 > {log}
        """


rule create_bracken_krona_plot:
    input:
        expand(str(OUTDIR/"kraken2/{sample}.bracken.krona"), sample=SAMPLES),
    output:
        krona_html=report(OUTDIR/"kraken2/all_samples.bracken.krona.html",
                          category="Taxonomic profiling",
                          caption="../../report/bracken_krona.rst"),
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


rule filter_bracken:
    input:
        bracken=OUTDIR/"kraken2/{sample}.{level}.bracken",
    output:
        filtered=OUTDIR/"kraken2/{sample}.{level,[DPOCFGS]}.filtered.bracken",  
    log:
        str(LOGDIR/"kraken2/{sample}.{level}.filter_bracken.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        filter_bracken="scripts/KrakenTools/filter_bracken.out.py",
        include=kraken2_config["filter_bracken"]["include"],
        exclude=kraken2_config["filter_bracken"]["exclude"],
    shell:
        """
        {params.filter_bracken} \
            --input-file {input.bracken} \
            --output {output.filtered} \
            {params.include} \
            {params.exclude} \
            2>&1 > {log}
        """


rule join_bracken_filtered:
    input:
        bracken=expand(str(OUTDIR/"kraken2/{sample}.{{level}}.filtered.bracken"), sample=SAMPLES),
    output:
        table=report(OUTDIR/"kraken2/all_samples.{level,[DPCOFGS]}.filtered.bracken.txt",
               category="Taxonomic profiling",
               caption="../../report/bracken_table_filtered.rst"),
    log:
        str(LOGDIR/"kraken2/join_bracken_tables.{level}.log")
    threads:
        1
    conda:
        "../../envs/stag-mwc.yaml"
    singularity:
        "shub://ctmrbio/stag-mwc:stag-mwc"
    params:
        value_column="fraction_total_reads",
        feature_column="name",
    shell:
        """
        scripts/join_tables.py \
            --outfile {output.table} \
            --value-column {params.value_column} \
            --feature-column {params.feature_column} \
            {input.bracken} \
            2>&1 > {log}
        """

# vim: syntax=python expandtab
# Rules to filter host sequences from metagenomic reads
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

rh_config = config["remove_host"]
if config["host_removal"]:
    db_path = Path(config["remove_host"]["db_path"])
    if not Path(db_path/"taxo.k2d").is_file():
        err_message = "Cannot find database for host sequence removal at: '{}/*.k2d'!\n".format(db_path)
        err_message += "Specify path to folder containing Kraken2 database for host removal in config.yaml.\n"
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    host_proportions = str(OUTDIR/"host_removal/host_proportions.txt")
    if rh_config["keep_fastq"]:
        all_outputs.extend(filtered_host)
    all_outputs.append(host_proportions)

    citations.add(publications["Kraken2"])

    localrules:
        plot_proportion_host

    rule remove_host:
        """Filter reads matching host database using Kraken2."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz" if rh_config["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_1.fq.gz"),
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz" if rh_config["keep_fastq"] else temp(OUTDIR/"host_removal/{sample}_2.fq.gz"),
            host1=OUTDIR/"host_removal/{sample}.host_1.fq.gz" if rh_config["keep_host_fastq"] else temp(OUTDIR/"host_removal/{sample}.host_1.fq.gz"),
            host2=OUTDIR/"host_removal/{sample}.host_2.fq.gz" if rh_config["keep_host_fastq"] else temp(OUTDIR/"host_removal/{sample}.host_2.fq.gz"),
            kraken=OUTDIR/"host_removal/{sample}.kraken" if rh_config["keep_kraken"] else temp(OUTDIR/"host_removal/{sample}.kraken"),
            kreport=OUTDIR/"host_removal/{sample}.kreport" if rh_config["keep_kreport"] else temp(OUTDIR/"host_removal/{sample}.kreport"),
        log:
            stderr=str(LOGDIR/"host_removal/{sample}.kraken2.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads:
            cluster_config["remove_host"]["n"] if "remove_host" in cluster_config else 8
        params:
            db=rh_config["db_path"],
            confidence=rh_config["confidence"],
            extra=rh_config["extra"],
            classified=lambda w: f"{OUTDIR}/host_removal/{w.sample}.host#.fq",
            unclassified=lambda w: f"{OUTDIR}/host_removal/{w.sample}#.fq",
            fq_to_compress=lambda w: f"{OUTDIR}/host_removal/{w.sample}*.fq",
        shell:
            """
            kraken2 \
                --db {params.db} \
                --threads {threads} \
                --output {output.kraken} \
                --classified-out {params.classified} \
                --unclassified-out {params.unclassified} \
                --report  {output.kreport} \
                --paired \
                --confidence {params.confidence} \
                {params.extra} \
                {input.read1} {input.read2} \
                2> {log.stderr}
            pigz \
                --processes {threads} \
                --verbose \
                --force \
                {params.fq_to_compress} \
                2>> {log.stderr}
            """


    rule plot_proportion_host:
        """Plot proportion of reads that matched the host DB."""
        input:
            expand(str(LOGDIR/"host_removal/{sample}.kraken2.log"), sample=SAMPLES)
        output:
            histogram=report(OUTDIR/"host_removal/host_histogram.pdf",
                       category="Preprocessing",
                       caption="../../report/host_histogram.rst"),
            barplot=report(OUTDIR/"host_removal/host_barplot.pdf",
                       category="Preprocessing",
                       caption="../../report/host_barplot.rst"),
            txt=report(OUTDIR/"host_removal/host_proportions.txt",
                       category="Preprocessing",
                       caption="../../report/host_proportions.rst"),
        log:
            str(LOGDIR/"host_removal/proportion_host.log")
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        container:
            "oras://ghcr.io/ctmrbio/stag-mwc:stag-mwc"+singularity_branch_tag
        threads:
            1
        shell:
            """
            workflow/scripts/plot_proportion_kraken2.py \
                {input} \
                --histogram {output.histogram} \
                --barplot {output.barplot} \
                --table {output.txt} \
                2>&1 > {log}
            """

else:

    if not config["fastp"]["keep_output"]:
        err_message = "Set fastp keep_output in config.yaml to True in order to skip host removal.\n"
        err_message += "If you want to run host removal set remove_host in config.yaml to True"
        raise WorkflowError(err_message)
        
    
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(filtered_host)

    localrules:
        skip_remove_host,

    rule skip_remove_host:
        """Do not remove host sequences"""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
        log:
            stderr=str(LOGDIR/"host_removal/{sample}.stderr.log"),
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """

# vim: syntax=python expandtab
# Rules to filter host sequences from metagenomic reads
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

localrules:


rh_config = config["remove_host"]
if rh_config:
    db_path = Path(config["remove_host"]["db_path"])
    if not rh_config["db_path"] or not db_path.is_dir():
        err_message = "Cannot find database for host sequence removal at: '{}'!\n".format(db_path)
        err_message += "Specify path to folder containing Kraken2 database for host removal in config.yaml.\n"
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    #host_proportions = str(OUTDIR/"host_removal/host_proportions.tsv")
    all_outputs.extend(filtered_host)
    all_outputs.append(host_proportions)

    citations.add((
        "Wood, D.E., Lu, J., & Langmead, B. (2019).",
        "Improved metagenomic analysis with Kraken 2.",
        "Genome biology, 20(1), 257.",
        "https://doi.org/10.1186/s13059-019-1891-0"
    ))

    localrules:
        plot_proportion_host

    rule remove_host:
        """Filter reads matching host database using Kraken2."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
            host1=OUTDIR/"host_removal/{sample}_host_1.fq.gz",
            host2=OUTDIR/"host_removal/{sample}_host_2.fq.gz",
        log:
            statsfile=str(LOGDIR/"host_removal/{sample}.statsfile.txt"),
            stderr=str(LOGDIR/"host_removal/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            8
        params:
            db=rh_config["db"],
            confidence=rh_config["confidence",
            extra=rh_config["extra"],
            classified=lambda w: f"{OUTDIR}/host_removal/{w.sample}_#.fq.gz",
            unclassified=lambda w: f"{OUTDIR}/host_removal/{w.sample}_#.host.fq.gz",
        shell:
            """
            kraken2 \
                --db {params.db} \
                --threads {threads} \
                --output {output.kraken} \
                --classified-out {params.classified} \
                --unclassified-out {params.unclassified} \
                --paired \
                --confidence {params.confidence} \
                {params.extra} \
                2> {log.stderr}
            """


    rule plot_proportion_host:
        """Plot proportion of reads that matched the host DB."""
        input:
            expand(str(LOGDIR/"host_removal/{sample}.statsfile.txt"), sample=SAMPLES)
        output:
            histogram=report(OUTDIR/"host_removal/host_histogram.pdf",
                       category="Proportion host reads",
                       caption="../../report/host_histogram.rst"),
            barplot=report(OUTDIR/"host_removal/host_barplot.pdf",
                       category="Proportion host reads",
                       caption="../../report/host_barplot.rst"),
            tsv=report(OUTDIR/"host_removal/host_proportions.tsv",
                       category="Proportion host reads",
                       caption="../../report/host_proportions.rst"),
        log:
            str(LOGDIR/"host_removal/proportion_host.log")
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            1
        params:
            unambigous=lambda _: "--unambigous" if rh_config["plot_unambigous"] else ""
        shell:
            """
            scripts/plot_proportion_host.py \
                {input} \
                {params.unambigous} \
                --histogram {output.histogram} \
                --barplot {output.barplot} \
                --table {output.tsv} \
                2>&1 > {log}
            """

else:
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

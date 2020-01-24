# vim: syntax=python expandtab
# Rules to filter host sequences from metagenomic reads
# TODO: Remove superfluous str conversions when Snakemake is pathlib compatible.
from pathlib import Path 

from snakemake.exceptions import WorkflowError

localrules:
    download_hg19,


rule download_hg19:
    """Download masked hg19 from: 
    https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk"""
    output:
        OUTDIR/"hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"hg19"
    shell:
        """
        scripts/download_from_gdrive.py \
            -o {output}.gz \
            0B3llHR93L14wd0pSSnFULUlhcUk \
        && \
        gunzip {output}.gz
    """


rule index_hg19:
    """Create BBMap index of hg19 fasta file."""
    input:
        DBDIR/"hg19/hg19_main_mask_ribo_animal_allplant_allfungus.fa",
    output:
        DBDIR/"hg19/ref/genome/1/chr1.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr2.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr3.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr4.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr5.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr6.chrom.gz",
        DBDIR/"hg19/ref/genome/1/chr7.chrom.gz",
        DBDIR/"hg19/ref/genome/1/info.txt",
        DBDIR/"hg19/ref/genome/1/scaffolds.txt.gz",
        DBDIR/"hg19/ref/genome/1/summary.txt",
        DBDIR/"hg19/ref/index/1/chr1-3_index_k13_c2_b1.block",
        DBDIR/"hg19/ref/index/1/chr1-3_index_k13_c2_b1.block2.gz",
        DBDIR/"hg19/ref/index/1/chr4-7_index_k13_c2_b1.block",
        DBDIR/"hg19/ref/index/1/chr4-7_index_k13_c2_b1.block2.gz",
    conda:
        "../../envs/stag-mwc.yaml"
    params:
        dbdir=DBDIR/"hg19"
    shell:
        """
        bbmap.sh \
            ref={input} \
            path={params.dbdir}
        """


if config["host_removal"]:
    db_path = Path(config["remove_host"]["db_path"])
    if not Path(db_path/"ref").exists():
        err_message = "Cannot find database for host sequence removal at: '{}'!\n".format(db_path)
        err_message += "Specify path to folder containing BBMap index in config.yaml.\n"
        err_message += "Run 'snakemake index_hg19' to download and create a BBMap index in '{dbdir}'".format(dbdir=DBDIR/"hg19")
        raise WorkflowError(err_message)

    # Add final output files from this module to 'all_outputs' from the main
    # Snakefile scope. SAMPLES is also from the main Snakefile scope.
    filtered_host = expand(str(OUTDIR/"host_removal/{sample}_{readpair}.fq.gz"),
            sample=SAMPLES,
            readpair=[1,2])
    host_proportions = str(OUTDIR/"host_removal/host_proportions.tsv")
    all_outputs.extend(filtered_host)
    all_outputs.append(host_proportions)

    citations.add((
        "Bushnell, B. (2016).",
        "BBMap short read aligner.",
        "University of California, Berkeley, California.",
        "Available online at: http://sourceforge.net/projects/bbmap.",
    ))

    localrules:
        plot_proportion_host

    rh_config = config["remove_host"]
    rule remove_host:
        """Filter reads matching host database."""
        input:
            read1=OUTDIR/"fastp/{sample}_1.fq.gz",
            read2=OUTDIR/"fastp/{sample}_2.fq.gz",
        output:
            read1=OUTDIR/"host_removal/{sample}_1.fq.gz",
            read2=OUTDIR/"host_removal/{sample}_2.fq.gz",
            host=OUTDIR/"host_removal/{sample}_host.fq.gz",
        log:
            statsfile=str(LOGDIR/"host_removal/{sample}.statsfile.txt"),
            stderr=str(LOGDIR/"host_removal/{sample}.stderr.log"),
        shadow:
            "shallow"
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            16
        params:
            minid=rh_config["minid"],
            maxindel=rh_config["maxindel"],
            minhits=rh_config["minhits"],
            bandwidthratio=rh_config["bandwidthratio"],
            bandwidth=rh_config["bandwidth"],
            qtrim=rh_config["qtrim"],
            trimq=rh_config["trimq"],
            quickmatch=rh_config["quickmatch"],
            fast=rh_config["fast"],
            untrim=rh_config["untrim"],
        shell:
            """
            bbmap.sh \
                threads={threads} \
                in1={input.read1} \
                in2={input.read2} \
                path={rh_config[db_path]} \
                outu1={output.read1} \
                outu2={output.read2} \
                outm={output.host} \
                statsfile={log.statsfile} \
                minid={params.minid} \
                maxindel={params.maxindel} \
                minhits={params.minhits} \
                bandwidthratio={params.bandwidthratio} \
                bandwidth={params.bandwidth} \
                qtrim={params.qtrim} \
                trimq={params.trimq} \
                {params.quickmatch} \
                {params.fast} \
                {params.untrim} \
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
        conda:
            "../../envs/stag-mwc.yaml"
        threads:
            16
        params:
        shell:
            """
            ln -sv $(readlink -f {input.read1}) {output.read1} >> {log.stderr}
            ln -sv $(readlink -f {input.read2}) {output.read2} >> {log.stderr}
            """

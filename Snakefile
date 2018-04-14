# vim: syntax=python expandtab
#
# mwc - Metagenomic Workflow Collaboration
# Copyright (c) 2018 Authors
#
# Running snakemake -n in a clone of this repository should successfully
# execute a test dry-run of the workflow.
from snakemake.exceptions import WorkflowError

from sys import exit
import os.path

configfile: "config.yaml"
outdir = config["outdir"]
all_outputs = []

SAMPLES = set(glob_wildcards(config["inputdir"]+"/"+config["input_fn_pattern"]).sample)


if config["qc_reads"]:
    include: "rules/preproc/read_quality.smk"
    fastqc_input = expand("{outdir}/fastqc/{sample}_R{readpair}.{ext}",
            outdir=outdir,
            sample=SAMPLES,
            readpair=[1,2],
            ext=["zip", "html"])
    trimmed_qa = expand("{outdir}/trimmed_qa/{sample}_R{readpair}.trimmed_qa.fq.gz",
            outdir=outdir,
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(fastqc_input)
    all_outputs.extend(trimmed_qa)


if config["remove_human"]:
    include: "rules/preproc/remove_human.smk"
    filtered_human = expand("{outdir}/filtered_human/{sample}_R{readpair}.filtered_human.fq.gz",
            outdir=outdir,
            sample=SAMPLES,
            readpair=[1,2])
    all_outputs.extend(filtered_human)
    if not config["remove_human"]["hg19_path"]:
        print("WARNING: No hg19 database for human sequence removal specified!\n"
              "         Specify path to folder containing BBMap index of hg19 files in config.yaml.\n"
              "         Run 'snakemake index_hg19' to download and create a BBMap index in '{dbdir}/hg19'".format(dbdir=config["dbdir"]))


if config["mappers"]["bowtie2"]:
    include: "rules/mappers/bowtie2.smk"
    bowtie2_alignments = expand("{outdir}/bowtie2/{db_name}/{sample}.bam",
            outdir=outdir,
            sample=SAMPLES,
            db_name=bt2_db_name)
    bowtie2_stats = expand("{outdir}/bowtie2/{db_name}/{sample}.{stats}.txt",
            outdir=outdir,
            sample=SAMPLES,
            stats=["covstats", "rpkm"],
            db_name=bt2_db_name)
    all_outputs.extend(bowtie2_alignments)
    all_outputs.extend(bowtie2_stats)


if config["taxonomic_profile"]["centrifuge"]:
    include: "rules/taxonomic_profiling/centrifuge.smk"
    centrifuge = expand("{outdir}/centrifuge/{sample}.{output_type}.tsv",
            outdir=outdir,
            sample=SAMPLES,
            output_type=("centrifuge", "centrifuge_report"))
    all_outputs.extend(centrifuge)
    if not config["centrifuge"]["db_prefix"]:
        print("WARNING: No Centrifuge database specified!\n"
              "         Specify Centrifuge database prefix in the Centrifuge section of config.yaml.\n"
              "         Run 'snakemake download_centrifuge_database' to download a copy into '{dbdir}/centrifuge'\n".format(dbdir=config["dbdir"]) +
              "         If you do not want to run Centrifuge for taxonomic profiling, set centrifuge: False in config.yaml")


if config["taxonomic_profile"]["kaiju"]:
    include: "rules/taxonomic_profiling/kaiju.smk"
    kaiju = expand("{outdir}/kaiju/{sample}.kaiju", outdir=outdir, sample=SAMPLES)
    kaiju_reports = expand("{outdir}/kaiju/{sample}.kaiju.summary.species", outdir=outdir, sample=SAMPLES)
    all_outputs.extend(kaiju)
    all_outputs.extend(kaiju_reports)
    if not all([config["kaiju"]["db"], config["kaiju"]["nodes"], config["kaiju"]["names"]]):
        print("WARNING: No Kaiju database specified!\n"
              "         Specify relevant paths in the kaiju section of config.yaml.\n"
              "         Run 'snakemake download_kaiju_database' to download a copy into '{dbdir}/kaiju'\n".format(dbdir=config["dbdir"]) + 
              "         If you do not want to run Kaiju for taxonomic profiling, set 'kaiju: False' in config.yaml")


if config["taxonomic_profile"]["metaphlan2"]:
    include: "rules/taxonomic_profiling/metaphlan2.smk"
    mpa_outputs = expand("{outdir}/metaphlan2/{sample}.{output_type}",
            outdir=outdir,
            sample=SAMPLES,
            output_type=("bowtie2.bz2", "metaphlan2.txt"))
    mpa_combined = expand("{outdir}/metaphlan2/all_samples.metaphlan2.txt",
            outdir=outdir)
    all_outputs.extend(mpa_outputs)
    all_outputs.extend(mpa_combined)
    if not any([config["metaphlan2"]["mpa_pkl"], config["metaphlan2"]["bt2_db_prefix"]]):
        print("WARNING: No MetaPhlAn2 database specified!\n"
              "         Specify relevant paths in the metaphlan2 section of config.yaml.\n"
              "         Run 'snakemake build_metaphlan2_index' to download and build the default mpa_v20_m200 database in '{dbdir}/metaphlan2'\n".format(dbdir=config["dbdir"]) +
              "         If you do not want to run MetaPhlAn2 for taxonomic profiling, set metaphlan2: False in config.yaml")

if config["antibiotic_resistance"]:
    include: "rules/antibiotic_resistance/megares.smk"
    megares_outputs = expand("{outdir}/megares/{sample}.{output_type}",
            outdir=outdir,
            sample=SAMPLES,
            output_type=("sam.gz", "mapped_reads.fq.gz", "mhist.txt", "covstats.txt", "rpkm.txt"))
    all_outputs.extend(megares_outputs)
    if not config["megares"]["db_path"]:
        print("WARNING: No MEGARes database specified!\n"
              "         Specify the DB path in the megares section of config.yaml.\n"
              "         Run 'snakemake create_megares_index' to download and build a BBMap index in '{dbdir}/megares'\n".format(dbdir=config["dbdir"]) +
              "         If you do not want to map reads against MEGARes for antibiotic resistance gene detection, set antibiotic_resistance: False in config.yaml")



rule all:
    input:
        all_outputs

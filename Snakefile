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

#############################
# Pre-processing
#############################
if config["qc_reads"]:
    include: "rules/preproc/read_quality.smk"

if config["remove_human"]:
    include: "rules/preproc/remove_human.smk"

#############################
# Naive sample comparison
#############################
if config["sketch_compare"]:
    include: "rules/sketch_compare/sketch_compare.smk"

#############################
# Mappers
############################
if config["mappers"]["bbmap"]:
    include: "rules/mappers/bbmap.smk"

if config["mappers"]["bowtie2"]:
    include: "rules/mappers/bowtie2.smk"

#############################
# Taxonomic profiling
############################
if config["taxonomic_profile"]["centrifuge"]:
    include: "rules/taxonomic_profiling/centrifuge.smk"

if config["taxonomic_profile"]["kaiju"]:
    include: "rules/taxonomic_profiling/kaiju.smk"

if config["taxonomic_profile"]["metaphlan2"]:
    include: "rules/taxonomic_profiling/metaphlan2.smk"

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

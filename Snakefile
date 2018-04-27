# vim: syntax=python expandtab
#
#                   StaG
#   mwc - Metagenomic Workflow Collaboration
#         Copyright (c) 2018 Authors
#
# Running snakemake -n in a clone of this repository should successfully
# execute a test dry-run of the workflow.
from pathlib import Path

from snakemake.exceptions import WorkflowError

configfile: "config.yaml"
inputdir = Path(config["inputdir"])
outdir = Path(config["outdir"])
all_outputs = []

SAMPLES = set(glob_wildcards(inputdir/config["input_fn_pattern"]).sample)

#############################
# Pre-processing
#############################
if config["qc_reads"]:
    include: "rules/preproc/read_quality.smk"

if config["remove_human"]:
    include: "rules/preproc/remove_human.smk"

if config["assess_depth"]:
    include: "rules/preproc/bbcountunique.smk"

#############################
# Naive sample comparison
#############################
if config["sketch_compare"]:
    include: "rules/sketch_compare/sketch_compare.smk"

#############################
# Mappers
#############################
if config["mappers"]["bbmap"]:
    include: "rules/mappers/bbmap.smk"

if config["mappers"]["bowtie2"]:
    include: "rules/mappers/bowtie2.smk"

#############################
# Taxonomic profiling
#############################
if config["taxonomic_profile"]["centrifuge"]:
    include: "rules/taxonomic_profiling/centrifuge.smk"

if config["taxonomic_profile"]["kaiju"]:
    include: "rules/taxonomic_profiling/kaiju.smk"

if config["taxonomic_profile"]["metaphlan2"]:
    include: "rules/taxonomic_profiling/metaphlan2.smk"

#############################
# Antibiotic resistance
#############################
if config["antibiotic_resistance"]:
    include: "rules/antibiotic_resistance/megares.smk"


rule all:
    input:
        all_outputs

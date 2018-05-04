# vim: syntax=python expandtab
#
#    StaG Metagenomic Workflow Collaboration
#                 StaG-mwc
#         Copyright (c) 2018 Authors
#
# Running snakemake --use-conda -n in a clone of this repository should
# successfully execute a test dry run of the workflow.
from pathlib import Path

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
min_version("4.8.1")  # TODO: Bump version when Snakemake is pathlib compatible

stag_version = "0.1.1-dev"
print("="*60)
print("StaG Metagenomic Workflow Collaboration".center(60))
print("StaG-mwc".center(60))
print(stag_version.center(60))
print("="*60)

configfile: "config.yaml"
INPUTDIR = Path(config["inputdir"])
OUTDIR = Path(config["outdir"])
LOGDIR = Path(config["logdir"])
DBDIR = Path(config["dbdir"])
all_outputs = []

SAMPLES = set(glob_wildcards(INPUTDIR/config["input_fn_pattern"]).sample)

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

localrules: all

rule all:
    input:
        all_outputs

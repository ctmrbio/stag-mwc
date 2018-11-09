# vim: syntax=python expandtab
#
#    StaG Metagenomic Workflow Collaboration
#                 StaG-mwc
#         Copyright (c) 2018 Authors
#
# Running `snakemake --use-conda --dryrun` in a clone of this repository should
# successfully execute a test dry run of the workflow.
from pathlib import Path

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
min_version("5.2.0")  # TODO: Bump version requirement when Snakemake is pathlib compatible

stag_version = "0.3.0-dev"
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
include: "rules/preproc/read_quality.smk"
include: "rules/preproc/remove_human.smk"
include: "rules/preproc/bbcountunique.smk"

#############################
# Naive sample comparison
#############################
include: "rules/sketch_compare/sketch_compare.smk"

#############################
# Mappers
#############################
include: "rules/mappers/bbmap.smk"
include: "rules/mappers/bowtie2.smk"

#############################
# Taxonomic profiling
#############################
include: "rules/taxonomic_profiling/centrifuge.smk"
include: "rules/taxonomic_profiling/kaiju.smk"
include: "rules/taxonomic_profiling/metaphlan2.smk"

#############################
# Assembly
#############################

#############################
# Functional profiling
#############################
include: "rules/functional_profiling/humann2.smk"

#############################
# Antibiotic resistance
#############################
include: "rules/antibiotic_resistance/groot.smk"


localrules: all

rule all:
    input:
        all_outputs

# vim: syntax=python expandtab
#
#    StaG Metagenomic Workflow Collaboration
#                 StaG-mwc
#         Copyright (c) 2018-2021
#
# Full documentation available at:
# https://stag-mwc.readthedocs.org

from pathlib import Path
import copy
import subprocess
import textwrap

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
min_version("5.5.4")

from rules.publications import publications
from scripts.common import UserMessages, SampleSheet

user_messages = UserMessages()

stag_version = "0.5.2"
singularity_branch_tag = "-develop"  # Replace with "-master" before publishing new version

configfile: "config.yaml"
report: "report/workflow.rst"
citations = {publications["StaG"], publications["Snakemake"]}

INPUTDIR = Path(config["inputdir"])
OUTDIR = Path(config["outdir"])
LOGDIR = Path(config["logdir"])
TMPDIR = Path(config["tmpdir"])
DBDIR = Path(config["dbdir"])
all_outputs = []

if config["samplesheet"]:
    samplesheet = SampleSheet(config["samplesheet"], keep_local=config["keep_local"], endpoint_url=config["s3_endpoint_url"])
    SAMPLES = samplesheet.samples
    INPUT_read1 = lambda w: samplesheet.sample_info[w.sample]["read1"]
    INPUT_read2 = lambda w: samplesheet.sample_info[w.sample]["read2"]
else:
    SAMPLES = set(glob_wildcards(INPUTDIR/config["input_fn_pattern"]).sample)
    INPUT_read1 = INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="1"),
    INPUT_read2 = INPUTDIR/config["input_fn_pattern"].format(sample="{sample}", readpair="2")

onstart:
    print("\n".join([
        "="*60, 
        "StaG Metagenomic Workflow Collaboration".center(60),
        "StaG-mwc".center(60),
        stag_version.center(60),
        "="*60,
        ])
    )

    if len(SAMPLES) < 1:
        raise WorkflowError("Found no samples! Check input file options in config.yaml")
    else:
        if config["samplesheet"]:
            print(f"Found these samples in '{config['samplesheet']}':\n{SAMPLES}")
        else:
            print(f"Found these samples in '{config['inputdir']}' using input filename pattern '{config['input_fn_pattern']}':\n{SAMPLES}")


#############################
# Pre-processing
#############################
include: "rules/preproc/read_quality.smk"
include: "rules/preproc/host_removal.smk"
include: "rules/preproc/preprocessing_summary.smk"

#############################
# Naive sample analyses
#############################
include: "rules/naive/sketch_compare.smk"
include: "rules/naive/bbcountunique.smk"

#############################
# Taxonomic profiling
#############################
include: "rules/taxonomic_profiling/kaiju.smk"
include: "rules/taxonomic_profiling/kraken2.smk"
include: "rules/taxonomic_profiling/krakenuniq.smk"
include: "rules/taxonomic_profiling/metaphlan.smk"

#############################
# Strain level profiling
#############################
include: "rules/taxonomic_profiling/strainphlan.smk"

#############################
# Functional profiling
#############################
include: "rules/functional_profiling/humann.smk"

#############################
# Antibiotic resistance
#############################
include: "rules/antibiotic_resistance/groot.smk"
include: "rules/antibiotic_resistance/amrplusplus.smk"

#############################
# Mappers
#############################
include: "rules/mappers/bbmap.smk"
include: "rules/mappers/bowtie2.smk"

#############################
# Assembly
#############################
include: "rules/assembly/metawrap.smk"

#############################
# MultiQC
#############################
include: "rules/multiqc/multiqc.smk"


localrules: all

rule all:
    input:
        all_outputs
 

onerror:
    print("\n".join([
        "="*60, 
        "StaG Metagenomic Workflow Collaboration".center(60),
        "StaG-mwc".center(60),
        stag_version.center(60),
        "",
        "There was an error executing the workflow!".center(60),
        "",
        ])
    )

    user_messages.print_messages()

    print("\n".join([
        "",
        "Check log output, fix issues, then rerun workflow.".center(60),
        "="*60,
        ])
    )

    if config["email"]:
        shell("""mail -s "StaG-mwc run failed!" {config[email]} < {log}""")


onsuccess:
    print("\n".join([
        "="*60, 
        "StaG Metagenomic Workflow Collaboration".center(60),
        "StaG-mwc".center(60),
        stag_version.center(60),
        "",
        ])
    )

    user_messages.print_messages()

    print("\n".join([
        "",
        "Workflow completed successfully".center(60),
        "="*60,
        ])
    )

    if config["email"]:
        shell("""mail -s "StaG-mwc run completed" {config[email]} < {log}""")

    if config["report"]:
        from sys import argv
        from datetime import datetime

        report_datetime = datetime.now().strftime("%Y%m%d-%H%S")

        citation_filename = f"citations-{report_datetime}.rst"
        with open(citation_filename, "w") as citation_file:
            for citation in sorted(citations):
                citation_file.write("* "+citation+"\n")
        citations_link = Path("citations.rst")
        if citations_link.exists():
            Path("citations.rst").unlink()
        Path("citations.rst").symlink_to(citation_filename)

        unlock_call = copy.deepcopy(argv)
        unlock_call.append("--unlock")

        report_args = copy.deepcopy(argv)
        report_args.extend(["--report", f"{config['report']}-{report_datetime}.zip"])

        # Report generation doesn't work if --jobs 
        # or --use-singularity are specified,
        # so we strip all args related to these from argv
        # for report generation call
        skip = False
        report_call = []
        for arg in report_args:
            if arg == "--use-singularity":
                continue
            if arg == "--singularity-args":
                skip = True
                continue
            if arg == "--singularity-prefix":
                skip = True
                continue
            if arg == "--jobs":
                skip = True
                continue
            if skip:
                skip = False
                continue
            report_call.append(arg)

        subprocess.run(unlock_call)
        subprocess.run(report_call)


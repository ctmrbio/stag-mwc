# vim: syntax=python expandtab
#
#    StaG Metagenomic Workflow Collaboration
#                 StaG-mwc
#         Copyright (c) 2018-2020
#
# Full documentation available at:
# https://stag-mwc.readthedocs.org

from pathlib import Path
import textwrap

from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
min_version("5.5.4")

from rules.publications import publications
from scripts.common import UserMessages

user_messages = UserMessages()

stag_version = "0.4.3-dev"

configfile: "config.yaml"
report: "report/workflow.rst"
citations = {publications["StaG"], publications["Snakemake"]}

INPUTDIR = Path(config["inputdir"])
OUTDIR = Path(config["outdir"])
LOGDIR = Path(config["logdir"])
TMPDIR = Path(config["tmpdir"])
DBDIR = Path(config["dbdir"])
all_outputs = []

SAMPLES = set(glob_wildcards(INPUTDIR/config["input_fn_pattern"]).sample)

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
        raise WorkflowError("Found no samples! Check input file pattern and path in config.yaml")
    else:
        print(f"Found the following samples in inputdir using input filename pattern '{config['input_fn_pattern']}':\n{SAMPLES}")


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

        shell("{snakemake_call} --unlock".format(snakemake_call=argv[0]))
        shell("{snakemake_call} --report {report}-{datetime}.html".format(
            snakemake_call=argv[0],
            report=config["report"],
            datetime=report_datetime,
        ))
        shell("{snakemake_call} --report {report}-{datetime}.zip".format(
            snakemake_call=argv[0],
            report=config["report"],
            datetime=report_datetime,
        ))

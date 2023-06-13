#!/usr/bin/env python3
"""Summarize read counts from QC and host removal steps."""
__author__ = "CTMR, Fredrik Boulund"
__date__ = "2021-2023"
__version__ = "0.3"

from sys import argv, exit
from pathlib import Path
import json
import argparse

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

import pandas as pd

def parse_args():
    desc = f"{__doc__} Copyright (c) {__author__} {__date__}. Version v{__version__}"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--fastp", metavar="sample.json", nargs="+",
            help="fastp JSON output file.")
    parser.add_argument("--kraken2", metavar="sample.kraken2.log", nargs="+",
            help="Kraken2 log output.")
    parser.add_argument("--bowtie2", metavar="sample.samtools.fastq.log", nargs="+",
            help="Bowtie2 samtools fastq log output.")
    parser.add_argument("-o", "--output-table", metavar="TSV",
            default="read_processing_summary.txt",
            help="Filename of output table in tsv format [%(default)s].")
    parser.add_argument("-p", "--output-plot", metavar="PDF",
            default="",
            help="Filename of output table in PDF format [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_bowtie2_samtools_fastq_logs(logfiles):
    for logfile in logfiles:
        with open(logfile) as f:
            sample_name = Path(logfile).stem.split(".")[0]
            for line in f:
                if not line.startswith("[M::bam2fq_mainloop]"):
                    raise ValueError
                if "bam2fq_mainloop] processed" in line:
                    yield {
                        "Sample": sample_name,
                        "after_bowtie2_host_removal": int(int(line.split()[2])/2),  # /2 because bowtie2 counts both pairs
                    }


def parse_kraken2_logs(logfiles):
    for logfile in logfiles:
        with open(logfile) as f:
            sample_name = Path(logfile).stem.split(".")[0]
            for line in f:
                if " unclassified" in line:
                    yield {
                        "Sample": sample_name,
                        "after_kraken2_host_removal": int(line.strip().split()[0]),
                    }


def parse_fastp_logs(logfiles):
    for logfile in logfiles:
        sample_name = Path(logfile).stem.split(".")[0]
        with open(logfile) as f:
            fastp_data = json.load(f)
            yield {
                "Sample": sample_name, 
                "before_fastp": int(fastp_data["summary"]["before_filtering"]["total_reads"]/2),  # /2 because fastp counts both pairs
                "after_fastp": int(fastp_data["summary"]["after_filtering"]["total_reads"]/2), 
                "duplication": float(fastp_data["duplication"]["rate"]),
            }


if __name__ == "__main__":
    args = parse_args()

    dfs = {
        "fastp": pd.DataFrame(),
        "kraken2": pd.DataFrame(),
        "bowtie2": pd.DataFrame(),
    }

    if args.fastp:
        data_fastp = list(parse_fastp_logs(args.fastp))
        dfs["fastp"] = pd.DataFrame(data_fastp).set_index("Sample")
    if args.kraken2:
        data_kraken2 = list(parse_kraken2_logs(args.kraken2))
        dfs["kraken2"] = pd.DataFrame(data_kraken2).set_index("Sample")
    if args.bowtie2:
        data_bowtie2 = list(parse_bowtie2_samtools_fastq_logs(args.bowtie2))
        dfs["bowtie2"] = pd.DataFrame(data_bowtie2).set_index("Sample")

    df = pd.concat(dfs.values(), axis="columns")

    column_order = [
        "duplication",
        "before_fastp",
        "after_fastp",
        "after_kraken2_host_removal",
        "after_bowtie2_host_removal",
    ]
    final_columns = [c for c in column_order if c in df.columns]
    df = df[final_columns]

    df.to_csv(args.output_table, sep="\t")

    if args.output_plot:
        fig, ax = plt.subplots(figsize=(6, 5))
        df[final_columns[1:]]\
            .transpose()\
            .plot(kind="line", style=".-", ax=ax)
        ax.set_title("Reads passing through QC and host removal")
        ax.set_xlabel("Stage")
        ax.set_ylabel("Reads")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, -0.1))
        fig.savefig(args.output_plot, bbox_inches="tight")


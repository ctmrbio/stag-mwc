#!/usr/bin/env python3
"""Summarize read counts from all QC steps"""
__author__ = "CTMR, Fredrik Boulund"
__date__ = "2021"
__version__ = "0.2"

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
    parser.add_argument("--fastp", metavar="sample.stderr.txt", nargs="+",
            required=True,
            help="fastp stderr.txt")
    parser.add_argument("--kraken2", metavar="sample.kraken2.log", nargs="+",
            required=True,
            help="Kraken2 log output")
    parser.add_argument("-o", "--output-table", metavar="TSV",
            default="read_processing_summary.txt",
            help="filename of output table in tsv format [%(default)s].")
    parser.add_argument("-p", "--output-plot", metavar="PDF",
            default="read_processing_summary.pdf",
            help="filename of output table in PDF format [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_kraken2_logs(logfiles):
    for logfile in logfiles:
        with open(logfile) as f:
            sample_name = Path(logfile).stem.split(".")[0]
            for line in f:
                if " unclassified" in line:
                    yield {"Sample": sample_name, "after_host_removal": int(line.strip().split()[0])}


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

    read_qc = list(parse_fastp_logs(args.fastp))
    host_removal = list(parse_kraken2_logs(args.kraken2))

    df_read_qc = pd.DataFrame(read_qc).set_index("Sample")
    df_host_removal = pd.DataFrame(host_removal).set_index("Sample")

    final_column_order = [
        "duplication",
        "before_fastp",
        "after_fastp",
        "after_host_removal",
    ]
    df = df_read_qc.join(df_host_removal)[final_column_order]

    fig, ax = plt.subplots(figsize=(6, 5))
    df[final_column_order[1:]]\
        .div(df["before_fastp"], axis=0)\
        .transpose()\
        .plot(kind="line", style=".-", ax=ax)
    ax.set_title("Proportion reads passing through QC and host removal")
    ax.set_xlabel("Stage")
    ax.set_ylabel("% reads")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(0, -0.1))
    fig.savefig(args.output_plot, bbox_inches="tight")

    df.to_csv(args.output_table, sep="\t")



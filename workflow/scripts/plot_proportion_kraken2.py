#!/usr/bin/env python3
"""Plot proportion of reads classified by Kraken2."""
__author__ = "CTMR, Fredrik Boulund"
__date__ = "2020"
__version__ = "0.1"

from sys import argv, exit
from pathlib import Path
import argparse

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

import pandas as pd

def parse_args():
    desc = f"{__doc__} Copyright (c) {__author__} {__date__}. Version v{__version__}"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("log_output", metavar="LOG", nargs="+",
            help="Kraken2 log output (txt).")
    parser.add_argument("-H", "--histogram", dest="histogram", metavar="FILE",
            default="histogram.pdf",
            help="Filename of output histogram plot [%(default)s].")
    parser.add_argument("-b", "--barplot", dest="barplot", metavar="FILE",
            default="barplot.pdf",
            help="Filename of output barplot [%(default)s].")
    parser.add_argument("-t", "--table", dest="table", metavar="FILE",
            default="proportions.tsv",
            help="Filename of histogram data in TSV format [%(default)s].")
    parser.add_argument("-u", "--unclassified", dest="unclassified", action="store_true",
            default=False,
            help="Plot proportion unclassified reads instead of classified reads [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_kraken2_logs(logfiles, unclassified):
    search_string = "unclassified" if unclassified else " classified"
    for logfile in logfiles:
        with open(logfile) as f:
            sample_name = Path(logfile).stem.split(".")[0]
            for line in f:
                if search_string in line:
                    yield sample_name, float(line.split("(")[1].split(")")[0].strip("%"))


if __name__ == "__main__":
    options = parse_args()

    proportions = list(parse_kraken2_logs(options.log_output, options.unclassified))
    action = "unclassified" if options.unclassified else "classified"
    
    df = pd.DataFrame(proportions, columns=["Sample", "Proportion"]).set_index("Sample").rename(columns={"Proportion": f"% {action}"})
    print("Loaded {} proportions for {} samples.".format(df.shape[0], len(df.index.unique())))

    fig, ax = plt.subplots(figsize=(7, 5))
    df.plot(kind="hist", ax=ax, legend=None)
    ax.set_title(f"Proportion {action} reads")
    ax.set_xlabel(f"Proportion {action} reads")
    ax.set_ylabel("Frequency")
    fig.savefig(options.histogram, bbox_inches="tight")

    length_longest_sample_name = max([s for s in df.index.str.len()])
    fig2_width = max(5, length_longest_sample_name * 0.4)
    fig2_height = max(3, df.shape[0] * 0.25)

    fig2, ax2 = plt.subplots(figsize=(fig2_width, fig2_height))
    df.plot(kind="barh", ax=ax2, legend=None)
    ax2.set_title(f"Proportion {action} reads")
    ax2.set_xlabel(f"Proportion {action} reads")
    ax2.set_ylabel("Sample")
    fig2.savefig(options.barplot, bbox_inches="tight")

    df.to_csv(options.table, sep="\t")



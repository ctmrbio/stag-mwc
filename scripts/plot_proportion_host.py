#!/usr/bin/env python3
"""Plot proportion of reads mapped to host DB with BBMap."""
__author__ = "CTMR, Fredrik Boulund"
__date__ = "2020"
__version__ = "0.2"

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
    parser.add_argument("statsfile", metavar="statsfile", nargs="+",
            help="BBMap statsfile(s).")
    parser.add_argument("-H", "--histogram", dest="histogram", metavar="FILE",
            default="host_histogram.pdf",
            help="Filename of output histogram plot [%(default)s].")
    parser.add_argument("-b", "--barplot", dest="barplot", metavar="FILE",
            default="host_barplot.pdf",
            help="Filename of output barplot [%(default)s].")
    parser.add_argument("-t", "--table", dest="table", metavar="FILE",
            default="host_proportions.tsv",
            help="Filename of histogram data in TSV format [%(default)s].")
    parser.add_argument("-u", "--unambigous", dest="unambigous", action="store_true",
            default=False,
            help="Plot proportion unambigously mapped reads instead of all mapped reads [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_statsfile(statsfile, unambigous):
    search_string = "mapped"
    if unambigous: 
        search_string = "unambiguous"
    for statsfile in statsfile:
        with open(statsfile) as f:
            sample_name = Path(statsfile).stem.split(".")[0]
            for line in f:
                if line.startswith(search_string):
                    yield sample_name, float(line.split()[1].strip("%"))


if __name__ == "__main__":
    options = parse_args()

    proportions = list(parse_statsfile(options.statsfile, options.unambigous))
    
    df = pd.DataFrame(
            proportions,
            columns=["Sample", "Proportion"]).set_index("Sample").rename(columns={"Proportion": "% host"})
    print("Loaded {} proportions for {} samples.".format(
            df.shape[0],
            len(df.index.unique())))

    fig, ax = plt.subplots()
    df.plot(kind="hist", ax=ax)
    ax.set_title("Proportion host reads")
    ax.set_xlabel("Proportion host reads")
    ax.set_ylabel("Frequency")
    fig.savefig(options.histogram)


    fig2, ax2 = plt.subplots()
    df.plot(kind="barh", ax=ax2)
    ax2.set_title("Proportion host reads")
    ax2.set_xlabel("Proportion host reads")
    ax2.set_ylabel("Sample")
    fig2.savefig(options.barplot)

    df.to_csv(options.table, sep="\t")



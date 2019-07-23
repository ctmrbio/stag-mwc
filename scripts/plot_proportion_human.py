#!/usr/bin/env python
# Plot proportion of reads mapped to hg19 with BBMap.
# Fredrik Boulund 2019

from sys import argv, exit
from pathlib import Path
import argparse

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt

import pandas as pd

def parse_args():
    desc = "Plot histogram proportion human reads."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("statsfile", metavar="statsfile", nargs="+",
            help="BBMap statsfile(s).")
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE",
            default="proportion_human.pdf",
            help="Filename of output histogram plot [%(default)s].")
    parser.add_argument("-t", "--table", dest="table", metavar="FILE",
            default="proportion_human.tsv",
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
            columns=["Sample", "Proportion"]).set_index("Sample")
    print("Loaded {} proportions for {} samples.".format(
            df.shape[0],
            len(df.index.unique())))

    fig, ax = plt.subplots()
    df.plot(kind="hist", ax=ax)
    ax.set_title("Proportion human reads")
    ax.set_xlabel("Proportion human reads")
    ax.set_ylabel("Frequency")
    fig.savefig(options.outfile)

    df.to_csv(options.table, sep="\t")



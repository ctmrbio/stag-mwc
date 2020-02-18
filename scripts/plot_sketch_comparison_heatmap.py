#!/usr/bin/env python3
"""Plot heatmap of MinHash sketch comparisons by BBTool's sketch.sh/comparesketch.sh."""
__author__ = "Fredrik Boulund"
__date__ = "2018"
__version__ = "0.2.0"

from sys import argv, exit
from pathlib import Path
import argparse

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})

import pandas as pd
import seaborn as sns


def parse_args():
    desc = f"{__doc__} Version {__version__}. Copyright (c) {__author__} {__date__}."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("alltoall", metavar="alltoall",
            help="Output table from BBMap's comparesketch.sh in format=3.")
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE",
            default="all_vs_all.pdf",
            help="Filename of heatmap plot [%(default)s].")
    parser.add_argument("-c", "--clustered", dest="clustered", metavar="FILE",
            default="all_vs_all.clustered.pdf",
            help="Filename of clustered heatmap plot [%(default)s].")
    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


if __name__ == "__main__":
    options = parse_args()
    
    df = pd.read_table(
            options.alltoall, 
            index_col=False)
    print("Loaded data for {} sample comparisons.".format(df.shape[0]))

    similarity_matrix = df.pivot(index="#Query", 
            columns="Ref", values="ANI").fillna(100)

    corr = similarity_matrix.corr().fillna(0)
    g = sns.heatmap(corr, annot=True, fmt="2.1f", annot_kws={"fontsize": 2})
    g.set_title("Sample similarity")
    #g.set_xticklabels(g.get_xticklabels(), fontsize=4)  #WIP
    #g.set_yticklabels(g.get_yticklabels(), rotation=0, fontsize=4) #WIP
    g.set_ylabel("")
    g.set_xlabel("")
    g.figure.savefig(str(Path(options.outfile)))

    g = sns.clustermap(corr, annot=True, fmt="2.1f", annot_kws={"fontsize": 2})
    g.fig.suptitle("Sample similarity (clustered)")
    g.savefig(str(Path(options.clustered)))



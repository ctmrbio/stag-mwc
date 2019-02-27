#!/usr/bin/env python
# Plot heatmap of MinHash sketch comparisons by BBTool's sketch.sh/comparesketch.sh.
# Fredrik Boulund 2018

from sys import argv, exit
from pathlib import Path
import argparse

import matplotlib as mpl
mpl.use("agg")
mpl.rcParams.update({'figure.autolayout': True})

import pandas as pd
import seaborn as sns


def parse_args():
    desc = "Plot heatmap of sketch comparisons."
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
    
    colnames = ["Query", "Ref", "ANI", "QSize", "RefSize", "QBases"]
    df = pd.read_csv(
            options.alltoall, 
            sep="\t", 
            index_col=False, 
            skiprows=1,
            names=colnames)
    print("Loaded data for {} sample comparisons.".format(df.shape[0]))

    similarity_matrix = df.pivot(index="Query", 
            columns="Ref", values="ANI").fillna(100)

    corr = similarity_matrix.corr()
    g = sns.heatmap(corr, annot=True, fmt="2.1f", annot_kws={"fontsize": 5})
    g.set_title("Sample similarity")
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_ylabel("")
    g.set_xlabel("")
    g.figure.savefig(str(Path(options.outfile)))

    g = sns.clustermap(corr, annot=True, fmt="2.1f", annot_kws={"fontsize": 5})
    g.fig.suptitle("Sample similarity (clustered)")
    g.savefig(str(Path(options.clustered)))



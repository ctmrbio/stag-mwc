#!/usr/bin/env python
# Plot heatmap of MinHash sketch comparisons by BBTool's sketch.sh/comparesketch.sh.
# Fredrik Boulund 2018

from sys import argv, exit
import argparse

import matplotlib as mpl
mpl.use("agg")
import pandas as pd
import seaborn as sns


def parse_args():
    desc = "Plot heatmap of sketch comparisons."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("alltoall", metavar="alltoall",
            help="Output table from comparesketch.sh in format=3.")
    parser.add_argument("-o", "--outfile", dest="outfile", metavar="FILE",
            default="all_vs_all.pdf",
            help="Write heatmap plot to FILE [%(default)s].")
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

    g = sns.heatmap(similarity_matrix)
    g.set_title("Sample similarity")
    g.set_yticklabels(g.get_yticklabels(), rotation=0)
    g.set_ylabel("")
    g.set_xlabel("")

    g.get_figure().savefig(options.outfile)

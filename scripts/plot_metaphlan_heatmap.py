#!/usr/bin/env python3
"""Plot MetaPhlAn heatmap"""
__author__ = "Fredrik Boulund"
__date__ = "2022"
__version__ = "0.4"

from sys import argv, exit
from collections import defaultdict
from pathlib import Path
import argparse
import logging

import numpy as np
import pandas as pd
import seaborn as sns

TAXLEVELS = [
    "Kingdom", 
    "Phylum", 
    "Class", 
    "Order", 
    "Family", 
    "Genus", 
    "Species", 
    "Strain",
]

def parse_args():
    desc = f"{__doc__} v{__version__}. {__author__} (c) {__date__}."
    parser = argparse.ArgumentParser(description=desc, epilog="Version "+__version__)
    parser.add_argument("mpa_table",
            help="MetaPhlAn TSV table to plot.")
    parser.add_argument("-o", "--outfile-prefix", dest="outfile_prefix",
            default="mpa_heatmap",
            help="Outfile name [%(default)s]. "
                 "Will be appended with <taxonomic_level>_top<N>.{png,pdf}")
    parser.add_argument("-f", "--force", action="store_true",
            default=False,
            help="Overwrite output file if it already exists [%(default)s].")
    parser.add_argument("-l", "--level", 
            default="Species",
            choices=TAXLEVELS,
            help="Taxonomic level to summarize results for [%(default)s].")
    parser.add_argument("-t", "--topN", metavar="N",
            default=50,
            type=int,
            help="Only plot the top N taxa [%(default)s].")
    parser.add_argument("-p", "--pseudocount", metavar="P",
            default=-1,
            type=float,
            help="Use custom pseudocount, a negative value means to "
                 "autocompute a pseudocount as the median of the 0.01th "
                 "quantile across all samples [%(default)s].")
    parser.add_argument("-c", "--colormap",
            default="viridis",
            help="Matplotlib colormap to use [%(default)s].")
    parser.add_argument("-M", "--method",
            default="average",
            help="Linkage method to use, "
                 "see scipy.cluster.hierarchy.linkage docs [%(default)s].")
    parser.add_argument("-m", "--metric",
            default="euclidean",
            help="Distance metric to use, "
                 "see scipy.spatial.distance.pdist docs [%(default)s].")
    parser.add_argument("-L", "--loglevel", choices=["INFO", "DEBUG"],
            default="INFO",
            help="Set logging level [%(default)s].")
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    return parser.parse_args()


def parse_mpa_table(mpa_tsv):
    """Read joined MetaPhlAn tables into a Pandas DataFrame.

    * Convert ranks from first column into hierarchical MultiIndex
    """
    with open(mpa_tsv) as f:
        for lineno, line in enumerate(f):
            if line.startswith("#"):
                continue
            elif line.startswith("clade_name"):
                skiprows = lineno
                dropcols = ["clade_name"]
                break
            elif not line.startswith("#"):
                logger.error(f"Don't know how to process table")
                exit(3)

    df = pd.read_csv(mpa_tsv, sep="\t", skiprows=skiprows)

    logger.debug(df.head())

    lineages = df[dropcols[0]].str.split("|", expand=True)
    levels_present = TAXLEVELS[:len(lineages.columns)]  # Some tables don't have strain or species assignments
    df[levels_present] = lineages\
            .rename(columns={key: level for key, level in zip(range(len(levels_present)), levels_present)})
    mpa_table = df.drop(columns=dropcols).set_index(levels_present)

    logger.debug(f"Parsed data dimensions: {mpa_table.shape}")
    logger.debug(mpa_table.sample(10))

    return mpa_table


def extract_specific_level(mpa_table, level):
    """Extract abundances for a specific taxonomic level."""

    level_pos = mpa_table.index.names.index(level)

    if level_pos+1 == len(mpa_table.index.names):
        level_only = ~mpa_table.index.get_level_values(level).isnull()
        mpa_level = mpa_table.loc[level_only]
    else:
        level_assigned = ~mpa_table.index.get_level_values(level).isnull()
        next_level_assigned = ~mpa_table.index.get_level_values(mpa_table.index.names[level_pos+1]).isnull()
        level_only = level_assigned & ~next_level_assigned  # AND NOT 
        mpa_level = mpa_table.loc[level_only]

    ranks = mpa_table.index.names.copy()
    ranks.remove(level)
    mpa_level.index = mpa_level.index.droplevel(ranks)
    logger.debug(f"Table dimensions after extracting {level}-level only: {mpa_level.shape}")
    return mpa_level


def plot_clustermap(mpa_table, topN, pseudocount, colormap, method, metric):
    """Plot Seaborn clustermap."""

    top_taxa = mpa_table.median(axis=1).nlargest(topN)
    mpa_topN = mpa_table.loc[mpa_table.index.isin(top_taxa.index)]
    logger.debug(f"Table dimensions after extracting top {topN} taxa: {mpa_topN.shape}")

    if pseudocount < 0:
        pseudocount = mpa_topN.quantile(0.05).median() / 10
        if pseudocount < 1e-10:
            logger.warning(f"Automatically generated pseudocount is very low: {pseudocount}! "
                            "Setting pseudocount to 1e-10.")
            pseudocount = 1e-10
        logger.debug(f"Automatically generated pseudocount is: {pseudocount}")

    figwidth = mpa_topN.shape[1]
    figheight = 10+topN/5

    sns.set("notebook")
    clustergrid = sns.clustermap(
            mpa_topN.apply(lambda x: np.log10(x+pseudocount)),
            figsize=(figwidth, figheight),
            method=method,
            metric=metric,
            cmap=colormap,
            cbar_kws={"label": "$log_{10}$(abundance)"},
    )
    return clustergrid


def main(mpa_table, outfile_prefix, overwrite, level, topN, pseudocount, colormap, method, metric):
    mpa_table = parse_mpa_table(mpa_table)
    mpa_level = extract_specific_level(mpa_table, level)
    clustermap = plot_clustermap(mpa_level, topN, pseudocount, colormap, method, metric)

    outfile_png = Path(f"{outfile_prefix}.{level}_top{topN}.png")
    outfile_pdf = Path(f"{outfile_prefix}.{level}_top{topN}.pdf")
    if (outfile_png.exists() or outfile_pdf.exists()) and not overwrite:
        logger.error(f"Output file {outfile_png} or {outfile_pdf} already exists and --force is not set.")
        exit(2)

    clustermap.ax_heatmap.set_xticklabels(clustermap.ax_heatmap.get_xticklabels(), rotation=90)

    clustermap.savefig(outfile_png)
    clustermap.savefig(outfile_pdf)


if __name__ == "__main__":
    args = parse_args()
    logger = logging.getLogger(__name__)
    loglevels = {"INFO": logging.INFO, "DEBUG": logging.DEBUG}
    logging.basicConfig(format='%(levelname)s: %(message)s', level=loglevels[args.loglevel])

    main(
        args.mpa_table, 
        args.outfile_prefix, 
        args.force,
        args.level,
        args.topN,
        args.pseudocount,
        args.colormap,
        args.method,
        args.metric,
    )

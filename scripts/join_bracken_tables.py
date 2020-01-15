#!/usr/bin/env python3
""" Join Bracken tables """
__author__ = "Fredrik Boulund"
__date__ = "2020"
__version__ = "0.2"

from sys import argv, exit
from functools import reduce, partial
from pathlib import Path
import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("BRACKEN", nargs="+",
            help="Bracken TSV table with columns: name, taxonomy_id, taxonomy_lvl, fraction_total_reads.")
    parser.add_argument("-c", "--value-column", dest="value_column",
            default="fraction_total_reads",
            help="Value column to use [%(default)s].")
    parser.add_argument("-o", "--outfile", dest="outfile",
            default="bracken_joined.tsv",
            help="Outfile name [%(default)s].")
    parser.add_argument("-n", "--fillna", dest="fillna", metavar="FLOAT",
            default=0.0,
            type=float,
            help="Fill NA values in merged table with FLOAT [%(default)s].")
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    return parser.parse_args()


def main(bracken_files, value_column, outfile, fillna):
    tables = []
    for bracken_file in bracken_files:
        sample_name = Path(bracken_file).name.split(".")[0]
        tables\
            .append(pd.read_csv(bracken_file, sep="\t")\
            .set_index(["name", "taxonomy_id"])\
            .rename(columns={value_column: sample_name})\
            .loc[:, [sample_name]])  # Ugly hack to get a single-column DataFrame

    df = tables[0]
    for table in tables[1:]:
        df = df.join(table, how="outer")
    df.fillna(fillna, inplace=True)

    df.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    args = parse_args()
    if len(args.BRACKEN) < 2:
        print("Need at least two tables to merge!")
        exit(1)
    main(args.BRACKEN, args.value_column, args.outfile, args.fillna)

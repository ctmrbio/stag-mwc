#!/usr/bin/env python3
"""Join per sample feature tables into a large combined table."""
__author__ = "Fredrik Boulund"
__date__ = "2020-2022"
__version__ = "1.1"

from sys import argv, exit
from functools import reduce, partial
from pathlib import Path
import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("TABLE", nargs="+",
            help="TSV table with columns headers.")
    parser.add_argument("-f", "--feature-column", dest="feature_column",
            default="name",
            help="Column header of feature column to use, "
                 "typically containing taxa names. "
                 "Select several columns by separating with comma (e.g. name,taxid) "
                 "[%(default)s].")
    parser.add_argument("-c", "--value-column", dest="value_column",
            default="fraction_total_reads",
            help="Column header of value column to use, "
                 "typically containing counts or abundances [%(default)s].")
    parser.add_argument("-o", "--outfile", dest="outfile",
            default="joined_table.tsv",
            help="Outfile name [%(default)s].")
    parser.add_argument("-n", "--fillna", dest="fillna", metavar="FLOAT",
            default=0.0,
            type=float,
            help="Fill NA values in merged table with FLOAT [%(default)s].")
    parser.add_argument("-s", "--skiplines", dest="skiplines", metavar="N",
            default=0,
            type=int,
            help="Skip N lines before parsing header (e.g. for files "
                 "containing comments before the real header) [%(default)s].")
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    return parser.parse_args()


def main(table_files, feature_column, value_column, outfile, fillna, skiplines):
    feature_columns = feature_column.split(",")

    tables = []
    for table_file in table_files:
        sample_name = Path(table_file).name.split(".")[0]
        tables\
            .append(pd.read_csv(table_file, sep="\t", skiprows=skiplines)\
            .set_index(feature_columns)\
            .rename(columns={value_column: sample_name})\
            .loc[:, [sample_name]])  # Ugly hack to get a single-column DataFrame

    df = tables[0]
    for table in tables[1:]:
        df = df.join(table, how="outer")
    df.fillna(fillna, inplace=True)

    df.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    args = parse_args()
    if len(args.TABLE) < 2:
        print("Need at least two tables to merge!")
        exit(1)
    main(args.TABLE, args.feature_column, args.value_column, args.outfile, args.fillna, args.skiplines)

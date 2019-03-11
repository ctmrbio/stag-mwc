#!/usr/bin/env python3
"""Make count table of all samples from BBMap pileup.sh rpkm tables, 
and multi-column annotation file.
Create one output file per given annotation column.
"""
__author__ = "Fredrik Boulund"
__date__ = "2019-03-07"
__version__ = "2.0.0"

from sys import argv, exit, stderr
from collections import defaultdict
from pathlib import Path
import argparse
import logging
import csv

logging.basicConfig(format="%(levelname)s: %(message)s")


def parse_args():
    desc = "{} Version v{}. Copyright (c) {}.".format(__doc__, __version__, __author__, __date__[:4])
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument("RPKM", nargs="+",
            help="RPKM file(s) from BBMap pileup.sh.")
    parser.add_argument("-c", "--columns", dest="columns",
            default="",
            help="Comma-separated list of column names to include [all columns].")
    parser.add_argument("-a", "--annotation-file", dest="annotation_file", 
            required=True,
            help="Two-column tab-separated annotation file.")
    parser.add_argument("-o", "--outdir", dest="outdir", metavar="DIR",
            default="",
            help="Directory for output files, will create one output file per selected column.")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_rpkm(rpkm_file):
    read_counts = {}
    with open(rpkm_file) as f:
        firstline = f.readline()
        if not firstline.startswith("#File"):
            logging.error("File does not look like a BBMap pileup.sh RPKM: %s", rpkm_file)
        _ = [f.readline() for l in range(4)] # Skip remaining header lines: #Reads, #Mapped, #RefSequences, Table header
        for line_no, line in enumerate(f, start=1):
            try:
                ref, length, bases, coverage, reads, RPKM, frags, FPKM = line.strip().split("\t")
            except ValueError:
                logging.error("Could not parse RPKM file line %s: %s", line_no, rpkm_file)
                continue
            if int(reads) != 0:
                ref = ref.split()[0]  # Truncate reference header on first space
                read_counts[ref] = int(reads)
    return read_counts


def parse_annotations(annotation_file):
    annotations = defaultdict(dict)
    with open(annotation_file) as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        for line in csv_reader:
            ref = list(line.values())[0].split()[0]  # Truncate reference header on first space
            for colname, value in list(line.items())[1:]:
                annotations[colname][ref] = value
    return annotations


def merge_counts(annotations, rpkms):
    output_table = {"Unknown": [0 for n in range(len(rpkms))]}
    for annotation in set(annotations.values()):
        output_table[annotation] = [0 for n in range(len(rpkms))]
    for idx, rpkm in enumerate(rpkms):
        for ref, count in rpkm.items():
            try:
                output_table[annotations[ref]][idx] += count
            except KeyError:
                logging.warning("Found no annotation for '%s', assigning to 'Unknown'", ref)
                output_table["Unknown"][idx] += count
    return output_table


def write_table(table_data, sample_names, outfile):
    with open(str(outfile), "w") as outf:
        header = "\t".join(["Annotation"] + [sample_name for sample_name in sample_names]) + "\n"
        outf.write(header)
        for ref, counts in table_data.items():
            outf.write("{}\t{}\n".format(ref, "\t".join(str(count) for count in counts)))


if __name__ == "__main__":
    args = parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    rpkms = []
    for rpkm_file in args.RPKM:
        rpkms.append(parse_rpkm(rpkm_file))

    annotations = parse_annotations(args.annotation_file)

    if args.columns:
        selected_columns = []
        for col in args.columns.split(","):
            if col in annotations:
                selected_columns.append(col)
            else:
                logging.warning("Column %s not found in annotation file!", col)
    else:
        selected_columns = list(annotations.keys())
    
    for selected_column in selected_columns:
        table_data = merge_counts(annotations[selected_column], rpkms)
        sample_names = [Path(fn).stem.split(".")[0] for fn in args.RPKM]

        table_filename = Path(args.outdir) / "counts.{}.tsv".format(selected_column)
        write_table(table_data, sample_names, table_filename)
        logging.debug("Wrote", table_filename)


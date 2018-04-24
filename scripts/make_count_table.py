#!/usr/bin/env python3
"""Make count table of multiple samples from BBMap pileup.sh rpkm tables, and two-column annotation file."""
__author__ = "Fredrik Boulund"
__date__ = "2018-04-24"
__version__ = "2.0.0"

from sys import argv, exit, stderr
from collections import defaultdict
import os.path
import argparse

import numpy as np

NP_DTYPE = np.uint16


def parse_args():
    desc = "{} Version v{}. Copyright (c) {}.".format(__doc__, __version__, __author__, __date__[:4])
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument("RPKM", nargs="+",
            help="RPKM file(s) from BBMap pileup.sh.")
    parser.add_argument("-a", "--annotations", required=True,
            help="Two-column tab-separated annotation file.")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_rpkm(rpkm_file):
    read_counts = {}
    with open(rpkm_file) as f:
        firstline = f.readline()
        if not firstline.startswith("#File"):
            print("ERROR: File does not look like a BBMap pileup.sh RPKM: {}".format(rpkm_file),
                    file=stderr)
        _ = [f.readline() for l in range(4)] # Skip remaining header lines: #Reads, #Mapped, #RefSequences, Table header
        for line_no, line in enumerate(f, start=1):
            try:
                ref, length, bases, coverage, reads, RPKM, frags, FPKM = line.strip().split("\t")
            except ValueError:
                print("ERROR: Could not parse RPKM file line {}:\n{}".format(line_no, rpkm_file),
                        file=stderr)
            yield ref, NP_DTYPE(reads)


def parse_annotations(annotation_file):
    annotations = {}
    with open(annotation_file) as f:
        for line_no, line in enumerate(f, start=1):
            try:
                ref, annotation = line.strip().split("\t")
            except ValueError:
                print("ERROR: Could not parse annotation file line {}:\n{}".format(line_no, annotation_file),
                        file=stderr)
            annotations[ref] = annotation
    return annotations


def merge_counts(annotations, rpkms):
    output_table = {"Unknown": np.zeros(len(rpkms), dtype=NP_DTYPE)}
    for annotation in set(annotations.values()):
        output_table[annotation] = np.zeros(len(rpkms), dtype=NP_DTYPE)
    for idx, rpkm_generators in enumerate(zip(*rpkms)):
        refs, counts = zip(*rpkm_generators)
        if not len(set(refs)) == 1:
            print("ERROR: RPKM files not in the same order, error on line {}:\n{}".format(idx, refs),
                    file=stderr)
        try:
            output_table[annotations[refs[0]]] += counts
        except KeyError:
            print("WARNING: Found no annotation for '{}', assigning to 'Unknown'".format(ref),
                    file=stderr)
            output_table["Unknown"][idx] += counts
    return output_table


def print_table(table_data, sample_names):
    header = "\t".join(["Annotation"] + [sample_name for sample_name in sample_names])
    print(header)
    for ref, counts in table_data.items():
        print("{}\t{}".format(ref, "\t".join(str(count) for count in counts)))


if __name__ == "__main__":
    args = parse_args()

    rpkms = []
    for rpkm_file in args.RPKM:
        rpkms.append(parse_rpkm(rpkm_file))

    annotations = parse_annotations(args.annotations)
    table_data = merge_counts(annotations, rpkms)
    sample_names = [os.path.basename(fn).split(".")[0] for fn in args.RPKM]
    print_table(table_data, sample_names)



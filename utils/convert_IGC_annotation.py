#!/usr/bin/env python3
# Fredrik Boulund
# 2018-04-24
"""Convert IGC annotation file to GTF or SAF."""

from sys import argv, exit, stderr
from functools import partial
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("IGC_ANNOTATION", 
            help="IGC annotation file. Typically 'IGC.annotation_OF.summary'. "
            "Can convert IGC_MEDUSA's CombinedNRGenesReduces.KEGG.txt to GTF "
            "format if used with '-f igc_medusa'.")
    parser.add_argument("-f", "--format",
            default="gtf",
            choices=["gtf", "saf", "igc_medusa"],
            help="Output format SAF or GTF [%(default)s].")
    parser.add_argument("-a", "--annotation-type", dest="annotation_type",
            default="eggnog",
            choices=["eggnog", "kegg", "phylum", "genus", "kegg_func", "eggnog_func"],
            help="Annotation type to put in SAF format GeneID column [%(default)s].")
    parser.add_argument("--gene-lengths", dest="gene_lengths",
            default="",
            help="Path to two-column tab-separated file containing sequence length (col 1)"
                 "and sequence name (col 2). Only used to convert igc_medusa, and if "
                 "not supplied, the script will set the gene length of all genes to "
                 "45297, which is the maximum observed gene length in igc_medusa "
                 "2018-04-24.")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def parse_igc_annotation(annotation_file):
    """
    Yield lines from IGC annotation summary file as tuples:
        (gene_id, gene_name, gene_length, gene_completeness,
            cohort, phylum, genus, kegg, eggnog, sample_occurence_freq,
            individual_occurence_freq, kegg_func_cat, eggnog_func_cat,
            cohort) = line.split("\t")
    """
    with open(annotation_file) as f:
        for line_no, line in enumerate(f, start=1):
            try:
                yield line.strip().split("\t")
            except ValueError:
                print("ERROR parsing line {}:\n{}".format(line_no, line), 
                        file=stderr)
                exit(1)

def parse_gene_lengths(gene_lengths_file):
    gene_lengths = {}
    with open(gene_lengths_file) as f:
        for line_no, line in enumerate(f, start=1):
            try:
                length, gene_name = line.strip().split()
                gene_lengths[gene_name] = length
            except ValueError:
                print("ERROR parsing gene_lengths, line {}:\n{}".format(line_no, line),
                        file=stderr)
                exit(1)
    return gene_lengths


def convert_medusa_to_gtf(annotation_lines, gene_lengths, max_seqlen=45297):
    """
    Convert IGC_medusa to very basic GTF format.
    max_seqlen is the length of longest sequence in the database,
    which will be used as the length of all seqs in the GTF.
    """
    gtf_line = "{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attribute}"
    for gene_id, kegg in annotation_lines:
        print(gtf_line.format(
                    seqname=gene_id,
                    source="IGC_medusa",
                    feature="gene",
                    start="1",
                    end=gene_lengths.get(gene_id, max_seqlen),
                    score=".",
                    strand="+",
                    frame="0",
                    attribute="; ".join(['gene_id "{}"'.format(kegg)]))
        )


def convert_to_saf(annotation_lines, annotation_type):
    saf_line = "{GeneID}\t{Chr}\t{Start}\t{End}\t{Strand}"
    print(saf_line.format(GeneID="GeneID", Chr="Chr", Start="Start", End="End", Strand="Strand"))
    for (gene_id, gene_name, gene_length, gene_completeness,
         cohort, phylum, genus, kegg, eggnog, sample_occurence_freq,
         individual_occurence_freq, kegg_func_cat, eggnog_func_cat,
         cohort) in annotation_lines:
        if annotation_type == "eggnog":
            GeneID = eggnog
        elif annotation_type == "kegg":
            GeneID = kegg
        elif annotation_type == "phylum":
            GeneID = phylum
        elif annotation_type == "genus":
            GeneID = genus
        elif annotation_type == "kegg_func":
            GeneID = kegg_func_cat
        elif annotation_type == "eggnog_func":
            GeneID = eggnog_func_cat
        else:
            GeneID = gene_id
        print(saf_line.format(
                    GeneID=GeneID,
                    Chr=gene_name,
                    Start="1",
                    End=gene_length,
                    Strand="+"))


def convert_to_gtf(annotation_lines):
    gtf_line = "{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attribute}"
    for (gene_id, gene_name, gene_length, gene_completeness,
         cohort, phylum, genus, kegg, eggnog, sample_occurence_freq,
         individual_occurence_freq, kegg_func_cat, eggnog_func_cat,
         cohort) in annotation_lines:
        print(gtf_line.format(
                    seqname=gene_name,
                    source="IGC",
                    feature="gene",
                    start="1",
                    end=gene_length,
                    score=".",
                    strand="+",
                    frame="0",
                    attribute="; ".join(['gene_id "{}"'.format(gene_name),
                                         'gene_completeness "{}"'.format(gene_completeness),
                                         'cohort "{}"'.format(cohort),
                                         'phylum "{}"'.format(phylum),
                                         'genus "{}"'.format(genus),
                                         'kegg "{}"'.format(kegg),
                                         'eggnog "{}"'.format(eggnog),
                                         'sample_occurence_freq "{}"'.format(sample_occurence_freq),
                                         'individual_occurence_freq "{}"'.format(individual_occurence_freq),
                                         'kegg_func_cat "{}"'.format(kegg_func_cat),
                                         'eggnog_func_cat "{}"'.format(eggnog_func_cat),
                                         'cohort "{}"'.format(cohort.replace(";", ","))])
            )
        )



if __name__ == "__main__":
    args = parse_args()
    if args.gene_lengths:
        gene_lengths = parse_gene_lengths(args.gene_lengths)
    else:
        gene_lengths = {}
    converters = {"gtf": convert_to_gtf, 
                  "saf": partial(convert_to_saf, annotation_type=args.annotation_type),
                  "igc_medusa": partial(convert_medusa_to_gtf, gene_lengths=gene_lengths)}
    annotation_lines = parse_igc_annotation(args.IGC_ANNOTATION)
    converters[args.format](annotation_lines)

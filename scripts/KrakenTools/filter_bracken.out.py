#!/usr/bin/env python
#################################################################
#filter_bracken_out.py allows users to filter Bracken output files 
#Copyright (C) 2019-2020 Jennifer Lu, jennifer.lu717@gmail.com
#
#This file is part of Kraken-Tools.

#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
#of the Software, and to permit persons to whom the Software is furnished to do
#so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
####################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 07/25/2019
#
#This program reads in a Bracken output file and user-provided taxids
#to either extract or exclude those taxids. Total fractions will be
#recalculated for read counts for the remaining taxids
#
#Parameters:
#   -h, --help...............show help message
#   -i X, --input-file X.....input bracken output filename
#   -o X, --output X.........output bracken-style output filename 
#   --include X..............include the specified taxids
#   --exclude X..............exclude the specified taxids 
#Input/Output file format (tab-delimited)
#   - name
#   - taxonomy ID
#   - level ID (S = species, G = genus, etc)
#   - Kraken assigned reads
#   - added reads with abundance reestimation
#   - total reads after abundance reestimation
#   - fraction of total reads
#######################################################################
import os, sys, argparse
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-file', dest='in_file', required=True,
        help='Input bracken OUTPUT file. [NOT the report file]')
    parser.add_argument('-o','--output','--output-file', dest='out_file', required=True,
        help='Output bracken OUTPUT file.')
    parser.add_argument('--include', dest='t_include',nargs='*', type=str, required=False, 
        help='List of taxonomy IDs to include in output [space-delimited] - default=All', default=[])
    parser.add_argument('--exclude', dest='t_exclude',nargs='*', type=str, required=False,
        help='List of taxonomy IDs to exclude in output [space-delimited] - default=None',default=[])
    args = parser.parse_args() 
    
    #CHECK#1: either taxonomy IDs are included or excluded
    if len(args.t_include) == 0 and len(args.t_exclude) == 0:
        sys.stderr.write("User must include at least one taxonomy ID to include or exclude\n")
        sys.stderr.write("Please specify either --include or --exclude\n")
        sys.exit(1)
    #CHECK#2: if both are specified, make sure none exists in both lists
    if len(args.t_include) > 0 and len(args.t_exclude) > 0:
        for val in args.t_include:
            if val in args.t_exclude:
                sys.stderr.write("%s cannot be in include AND exclude lists\n" % val)
                sys.exit(1)
    include = False
    exclude = False
    if len(args.t_include) > 0:
        include = True
    if len(args.t_exclude) > 0:
        exclude = True       
    #Process input file 
    sys.stdout.write(">> Reading Bracken output file: %s\n" % args.in_file)
    i_file = open(args.in_file,'r') 
    tot_reads = 0
    excl_reads = 0
    first = True
    first_line = ""
    save_taxid2all = {}
    save_taxid2reads = {} 
    for line in i_file:
        line = line.strip()
        #Check format
        if first:
            if line.split("\t") != ["name","taxonomy_id","taxonomy_lvl","kraken_assigned_reads","added_reads","new_est_reads","fraction_total_reads"]:
                sys.stderr.write("\t%s not in Bracken output format\n" % args.in_file)
                sys.exit(1)
            first = False
            firstline = line + "\n"
            continue
        #Get reads 
        l_vals = line.split("\t")
        if include:
            if l_vals[1] in args.t_include: 
                save_taxid2all[l_vals[1]] = l_vals[0:6]
                save_taxid2reads[l_vals[1]] = int(l_vals[5])
                tot_reads += int(l_vals[5])
            else:
                excl_reads += int(l_vals[5])
        elif exclude:
            if l_vals[1] not in args.t_exclude:
                save_taxid2all[l_vals[1]] = l_vals[0:6]
                save_taxid2reads[l_vals[1]] = int(l_vals[5])
                tot_reads += int(l_vals[5])
            else:
                excl_reads += int(l_vals[5])
    sys.stdout.write("\t%i reads remaining (%i reads excluded)\n" % (tot_reads, excl_reads))         
    i_file.close()

    #Write output file with updated fractions
    sys.stdout.write(">> Writing filtered Bracken output to %s\n" % args.out_file)
    o_file = open(args.out_file,'w')
    o_file.write(firstline)
    for [taxid, reads] in sorted(save_taxid2reads.items(), key=lambda kv: kv[1], reverse=True):
        new_fraction = float(reads) / float(tot_reads)
        for val in save_taxid2all[taxid]:
            o_file.write(val + "\t")
        o_file.write("%0.10f\n" % new_fraction)
    o_file.close()

#######################################################################
if __name__ == "__main__":
    main()
#######################END OF PROGRAM##################################

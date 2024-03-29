#!/usr/bin/env python
####################################################################
#fix_unmapped.py analyzes a text file of accession IDs and maps them
#to their respective taxonomy IDs given accession2taxid files. 
#Copyright 2019-2020 Jennifer Lu
#
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
####################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 06/13/2019
#
#This program takes in a text file of accession IDs and maps them
#to their respective taxonomy IDs given accession2taxid files. 
#
#Parameters:
#   -h, --help......................show help message
#   -i X, --input X, --input_file X.....input file with accessions
#   --accession2taxid X.........list all accession2taxid files 
#               :: list any number of files separated by spaces
#               :: files can be gzipped (with extension *.gz) or not
#   -o X, --output X, --output_file X...output file with accession/taxid mapping
#               :: output format: two tab-delimited columns, no header  
####################################################################
import os, sys, argparse
import gzip 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input','--input_file', 
        type=str, dest='in_file', required=True,
        help='Input file containing accession IDs to map. \
            Multi-column files accepted. Only accessions in \
            the first column will be mapped.')
    parser.add_argument('--accession2taxid', dest='ref_files', type=str, required=True,
        nargs='+', help='Accession2taxid reference mappings to search. \
            NCBI accession2taxid format required: 4 columns with accessions \
            in column 1 and taxonomy IDs in column 3.')
    parser.add_argument('-o','--output','--output_file',
        type=str, dest='out_file', required=True,
        help='Output file with 2 tab-delimited columns for accessions and taxids')
    parser.add_argument('-r','--remaining',required=False, 
        default='still_unmapped.txt',dest='rem_file', 
        help='Name of text file containing non-found accessions from input file')
    args = parser.parse_args()
    
    #STEP 1: READ IN ACCESSIONS 
    count_a = 0
    seq2taxid = {} 
    i_file = open(args.in_file,'r')
    sys.stdout.write(">> STEP 1: READING %s FOR ACCESSIONS\n" % args.in_file)
    sys.stdout.write("\t%i accessions read" % count_a)
    sys.stdout.flush()
    for line in i_file:
        line = line.strip()
        curr_a = line.split("\t")[0]  #if file has more than one column
        #Dont save duplicates
        if curr_a not in seq2taxid:
            count_a += 1
            seq2taxid[curr_a] = -1
            sys.stdout.write("\r\t%i accessions read" % count_a)
            sys.stdout.flush()
    sys.stdout.write("\r\t%i accessions read\n" % count_a)
    sys.stdout.flush()
    i_file.close() 
    
    #STEP 2: READ IN REFERENCE MAPS 
    count_found = 0
    sys.stdout.write(">> STEP 2: SEARCHING %i ACCESSION2TAXID FILES\n" % len(args.ref_files))
    sys.stdout.flush()
    for ref in args.ref_files: 
        #gzipped file
        if ref[-3:] == ".gz":
            r_file = gzip.open(ref,'r')
        else:
            r_file = open(ref,'r')
        sys.stdout.write("\tReading %s\n" % ref) 
        sys.stdout.write("\t\t%i / %i accessions found" % (count_found, count_a))
        sys.stdout.flush()
        for line in r_file: 
            line = line.strip()
            l_vals = line.split("\t")
            curr_a = l_vals[0]
            #Found accession - save
            if curr_a in seq2taxid and seq2taxid[curr_a] == -1:
                seq2taxid[curr_a] = int(l_vals[2])
                count_found += 1
                sys.stdout.write("\r\t\t%i / %i accessions found" % (count_found, count_a))
                sys.stdout.flush()
                if count_found == count_a:
                    break 
        sys.stdout.write("\r\t\t%i / %i accessions found\n" % (count_found, count_a))
        sys.stdout.flush()
        if count_found == count_a:
            break
    sys.stdout.write("\tfinished reading provided accession2taxid files\n")
    sys.stdout.flush()

    #STEP 3: OUTPUT FOUND ACCESSIONS AND REMAINING ACCESSIONS
    sys.stdout.write(">> STEP 3: PRINTING ACCESSION2TAXIDS TO %s\n" % args.out_file)
    sys.stdout.flush()
    if count_found < count_a:
        r_file = open(args.rem_file,'w')
    o_file = open(args.out_file,'w')
    for seq in seq2taxid:
        if seq2taxid[seq] == -1: 
            r_file.write("%s\n" % seq)
        else:
            o_file.write("%s\t%i\n" % (seq, seq2taxid[seq]))
    o_file.close()
    r_file.close()
if __name__ == "__main__":
    main()

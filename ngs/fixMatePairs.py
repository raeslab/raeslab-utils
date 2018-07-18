#!/usr/bin/env python
# Tested only with BioPython 1.54 and known not to work with python 1.63b
'''
A script to remove broken pairs from FASTQ files
'''
import sys
import os
import re
# Make sure BioPython version supports del dict[]
from Bio import SeqIO
from optparse import OptionParser

def main():
    usage = "%prog --in-fastq-1 file-1.fastq --in-fastq-2 file-2.fastq --out-fastq-1 output-1.fastq --out-fastq-2 output-2.fastq -b broken-pairs.fastq <trim-PE-info>"
    option_parser = OptionParser(usage=usage)
    option_parser.add_option(
        "--in-fastq-1", dest="i_fq_1",
        help="Input fastq file 1")
    option_parser.add_option(
        "--in-fastq-2", dest="i_fq_2",
        help="Input fastq file 2")
    option_parser.add_option(
        "--out-fastq-1", dest="o_fq_1",
        help="Output fastq file 1")
    option_parser.add_option(  
        "--out-fastq-2", dest="o_fq_2",
        help="Output fastq file 2")
    option_parser.add_option(
        "-b", "--out-broken-pairs", dest="o_broken_pairs",
        help="Output for broken pairs")
    option_parser.add_option(
        "--trim-PE-info", action="store_true", dest="trim_PE", default=False,
        help="If PE information /1 /2 in sequence ids should be removed in order to match sequence ids. Output names won't be changed")

    (options, args) = option_parser.parse_args()
    
    if((len(sys.argv) <= 1)):
        print __doc__
        option_parser.print_help()
        sys.exit(2)
    
    sample2Dict = SeqIO.index(options.i_fq_2, "fastq")
    # write paired ends     
    sample1OutHandle = open(options.o_fq_1, 'w')
    sample2OutHandle = open(options.o_fq_2, 'w')
    brokenOutHandle = open(options.o_broken_pairs, 'w')
    if options.trim_PE:
        for r in SeqIO.parse(options.i_fq_1, "fastq"):
            trimmedRecordId = re.sub("/1$", "/2", r.id)
            if trimmedRecordId in sample2Dict:
                SeqIO.write(r, sample1OutHandle, "fastq")
                SeqIO.write(sample2Dict[trimmedRecordId], sample2OutHandle, "fastq")
                del sample2Dict[trimmedRecordId]
            else:
                SeqIO.write(r, brokenOutHandle, "fastq")
    else:
        for r in SeqIO.parse(options.i_fq_1, "fastq"):
            if r.id in sample2Dict:
                SeqIO.write(r, sample1OutHandle, "fastq")
                SeqIO.write(sample2Dict[r.id], sample2OutHandle, "fastq")
                # delete matched pair [does not work on all version of biopython, tested only with 1.54. does not work with 1.63]
                del sample2Dict[r.id]
            else:
                # write to orphaned
                SeqIO.write(r, brokenOutHandle, "fastq")
    sample1OutHandle.close()                
    sample2OutHandle.close()
    # Handle broken pairs in 2.fastq
    for r in sample2Dict:
        SeqIO.write(sample2Dict[r], brokenOutHandle, "fastq")
    brokenOutHandle.close()


if __name__ == '__main__':
    main()

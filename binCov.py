#!/usr/bin/env python

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Calculates non-duplicate primary-aligned binned coverage of a chromosome from an input BAM file
"""

#Import libraries
import argparse
import sys
from collections import defaultdict, Counter, namedtuple
import pysam
import pybedtools

#Function to evaluate 

#Main function
def main():
    #Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ibam', type=pysam.AlignmentFile,
                        help='Input bam')
    parser.add_argument('chr', help='Contig to evaluate')
    parser.add_argument('cov_out', help='Output bed file of raw coverage')
    parser.add_argument('-n', '--norm_out', default = None,
    					help='Output bed file of normalized coverage')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size in bp (default: 1000)')
    parser.add_argument('-t', '--type', default='nucleotide',
                        help='Evaluate nucleotide or physical coverage (default: nucleotide)')
    parser.add_argument('-x', '--blacklist', default = None,
    	                help='BED file of regions to ignore')
    args = parser.parse_args()

    #Sanity check arguments
    if args.type != 'nucleotide' and args.type != 'physical':
    	print('ERROR: --type option must be either \'nucleotide\' or \'physical\'')
    	sys.exit(1)

    #Open outfiles
    fcovout = open(args.cov_out, 'w')
    if args.norm_out is not None:
	    fnormout = open(args.norm_out, 'w')

    #Run code here
    #TBD
    #TBD

    #Close outfiles
    fcovout.close()
    if args.norm_out is not None:
	    fnormout.close()


#Main block
if __name__ == '__main__':
    main()

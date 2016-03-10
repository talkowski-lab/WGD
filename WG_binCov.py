#!/usr/bin/env python

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Wrapper to run binCov.py on a list of specified contigs or on all contigs 
present in reference
"""

#Import libraries
import argparse
import pysam

#Main function
def main():
    #Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', type=pysam.AlignmentFile,
                        help='Input bam')
    parser.add_argument('ID', help='Sample ID')
    parser.add_argument('OUTDIR', help='Output directory for coverage files')
    parser.add_argument('-L', '--contigs', help='List of contigs to evaluate')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size in bp (default: 1000)')
    parser.add_argument('-m', '--mode', default='nucleotide',
                        choices = ['nucleotide', 'physical'],
                        help='Evaluate nucleotide or physical coverage '
                             '(default: nucleotide)')
    parser.add_argument('-x', '--blacklist', type=str,
                        help='BED file of regions to ignore')
    parser.add_argument('-v', '--overlap', nargs=1, type=float, default=0.05,
                           help='Maximum tolerated blacklist overlap before '
                                 'excluding bin')
    args = parser.parse_args()

    #
    coverage = binCov(args.bam, args.chr, args.binsize,
                      args.mode, args.overlap, args.blacklist)
    coverage.saveas(args.cov_out)
    call('sort -Vk1,1 -k2,2n -o ' + args.cov_out + ' ' + args.cov_out,
         shell=True)

    #Normalize coverage (if optioned) & write out
    if args.norm_out is not None:
        ncoverage = coverage.to_dataframe(names = 'chr start end cov'.split())
        medcov = ncoverage.loc[ncoverage['cov'] > 0, 'cov'].median()
        ncoverage['cov'] = ncoverage['cov'] / medcov
        ncoverage.to_csv(args.norm_out, sep='\t', index=False, header=False)
        call(' '.join(['sort -Vk1,1 -k2,2n -o', args.norm_out, 
                        args.norm_out]), shell=True)


#Main block
if __name__ == '__main__':
    main()

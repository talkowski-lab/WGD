#!/usr/bin/env python

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Calculates non-duplicate primary-aligned binned coverage
of a chromosome from an input BAM file
"""

#Import libraries
import argparse
import sys
from subprocess import call
import pysam
import pybedtools
import pandas as pd

#Define exception class for invalid coverage modes
class InvalidModeError(Exception):
    """Invalid coverage mode"""

#Function to return read or fragment intervals from pysam.AlignmentFile
def filter_mappings(bam, mode='nucleotide'):
    """
    Generates bed intervals from a bam for a specific chromosome corresponding
    either to read coordinates or fragment coordinates

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Input bam
    mode : str
        'physical' or 'nucleotide' (default: 'physical')

    Returns
    ------
    mappings : BedTool
        Read or fragment intervals (depending on mode)
    """

    #Sanity check mode
    if mode not in 'nucleotide physical'.split():
        raise InvalidModeError('Invalid mode: ' + mode + 
                               ' (options: nucleotide, physical)')

    #For nucleotide mode, return non-duplicate primary read mappings
    for read in bam:
        if not any([read.is_duplicate, read.is_unmapped,
                   read.is_secondary, read.is_supplementary]):
            if mode == 'nucleotide':
                yield '\t'.join([read.reference_name,
                                 str(read.reference_start),
                                 str(read.reference_end)]) + '\n'
            else:
                if read.is_read1 and read.is_proper_pair:
                    mincoord = min(str(read.reference_start),
                                   str(read.next_reference_start))
                    maxcoord = max(str(read.reference_start),
                                   str(read.next_reference_start))
                    yield '\t'.join([read.reference_name, 
                                     mincoord, maxcoord]) + '\n'


#Function to evaluate nucleotide or physical coverage
def binCov(bam, chr, binsize, mode='nucleotide', overlap=0.05, blacklist=None):
    """
    Generates non-duplicate, primary-aligned nucleotide or physical coverage 
    in regular bin sizes on a specified chromosome from a coordinate-sorted
    bamfile

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Input bam
    chr : string
        Chromosome to evaluate
    binsize : int
        Size of bins in bp
    mode : str
        Evaluate 'nucleotide' or 'physical' coverage
    overlap : float
        Maximum tolerated blacklist overlap before excluding bin
    blacklist : string
        Path to blacklist BED file

    Returns
    ------
    coverage : pybedtools.BedTool
        chr, start, end, coverage
    """
    
    #Create coverage bins and convert to BedTool
    maxchrpos = {d['SN']: d['LN'] for d in bam.header['SQ']}[chr]
    bin_starts = range(0, maxchrpos - binsize, binsize)
    bin_stops = range(binsize, maxchrpos, binsize)
    bins = []
    for i in range(0, len(bin_starts)-1):
        bins.append([chr, bin_starts[i], bin_stops[i]])
    bins = pybedtools.BedTool(bins)

    #Remove bins that have at least 5% overlap with blacklist by size
    blist = pybedtools.BedTool(blacklist[0])
    bins_filtered = bins.intersect(blist, v=True, f=overlap)

    #Filter bam
    mappings = filter_mappings(bam.fetch(chr), mode)
    bambed = pybedtools.BedTool(mappings)

    #Generate & return coverage
    coverage = bambed.coverage(bins_filtered, counts=True)
    return coverage


#Main function
def main():
    #Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', type=pysam.AlignmentFile,
                        help='Input bam')
    parser.add_argument('chr', help='Contig to evaluate')
    parser.add_argument('cov_out', help='Output bed file of raw coverage')
    parser.add_argument('-n', '--norm_out', nargs=1, type=str,
                        help='Output bed file of normalized coverage')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size in bp (default: 1000)')
    parser.add_argument('-m', '--mode', default='nucleotide',
                        choices = ['nucleotide', 'physical'],
                        help='Evaluate nucleotide or physical coverage '
                             '(default: nucleotide)')
    parser.add_argument('-x', '--blacklist', nargs=1, type=str,
                        help='BED file of regions to ignore')
    parser.add_argument('-v', '--overlap', nargs=1, type=float, default=0.05,
                           help='Maximum tolerated blacklist overlap before '
                                 'excluding bin')
    args = parser.parse_args()

    #Get coverage & write out
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
        call(' '.split(['sort -Vk1,1 -k2,2n -o', args.norm_out, 
                        args.norm_out]), shell=True)


#Main block
if __name__ == '__main__':
    main()

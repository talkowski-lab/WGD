#!/usr/bin/env python

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Calculates non-duplicate primary-aligned binned coverage of a chromosome from an input BAM file
"""

#Import libraries
import argparse
import sys
import pysam
import pybedtools

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
    mode : str, optional
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
    if mode == 'nucleotide':
        for read in bam:
            if not any([read.is_duplicate, read.is_unmapped,
                       read.is_secondary, read.is_supplementary]):
                yield '\t'.join([read.reference_name, str(read.reference_start),
                                 str(read.reference_end)]) + '\n'

#Function to evaluate nucleotide coverage
def nuc_binCov(bam, chr, binsize, blacklist=None):
	"""
    Generates non-duplicate, primary-aligned nucleotide coverage in regular bin
    sizes on a specified chromosome from a coordinate-sorted bamfile

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Input bam
    chr : string
        Chromosome to evaluate
    binsize : int
    	Size of bins in bp
    blacklist : string
    	Path to blacklist BED file

    Returns
    ------
    coverage : list
        chr, start, end, coverage
	"""

	#Define read filtering criteria
	#If True, read fails filtering
	def _filter(read):
		return (read.is_secondary or read.is_duplicate or
				read.is_supplementary or read.is_unmapped)


	#Subset bam for relevant reads and convert to BedTool
	subbam = pybedtools.BedTool().bam_to_bed(read for read in bam.fetch(str(chr)) if _filter(read) is False)
	
	#Instantiate coverage bins and convert to BedTool
	maxchrpos = {d['SN']: d['LN'] for d in bam.header['SQ']}[str(chr)]
	bin_starts = range(0, maxchrpos - binsize, binsize)
	bin_stops = range(binsize, maxchrpos, binsize)
	bins = []
	for i in range(0, len(bin_starts)-1):
		bins.append([chr, bin_starts[i], bin_stops[i]])
	bins = pybedtools.BedTool(bins)

	#Remove bins that have at least 5% overlap with blacklist by size
	bins_filtered = bins.intersect(blacklist, v=True, f=0.05)

	#Generate & return coverage
	coverage = subbam.coverage(bins_filtered, counts=True)
	return coverage

# #Function to evaluate physical coverage
# def phys_binCov(bam, chr, binsize, blacklist = None):
# 	"""
#     Generates non-duplicate, primary-aligned proper pair physical coverage
#     in regular bin sizes on a specified chromosome from a coordinate-sorted
#     bamfile

#     Parameters
#     ----------
#     bam : pysam.AlignmentFile
#         Input bam
#     chr : string
#         Chromosome to evaluate
#     binsize : int
#     	Size of bins in bp
#     blacklist : string
#     	Path to blacklist BED file

#     Returns
#     ------
#     coverage : list
#         chr, start, end, coverage
# 	"""

# 	#Define read filtering criteria
# 	#If True, read fails filtering
# 	def _filter(read):
# 		return (read.is_secondary or read.is_duplicate or
# 				read.is_supplementary or read.is_unmapped or
# 				read.mate_is_unmapped or (not read.is_proper_pair))

# 	#Subset bam for relevant reads and convert to BedTool
# 	nbam = bam.fetch(str(chr))
# 	subbam = pybedtools.BedTool(read for read in bam.fetch(str(chr)) if _filter(read) is False)
	
	# #Convert bam to bed of proper fragments
	# fragments =  subbam.bam_to_bed(bedpe=True)

	# #Instantiate coverage bins and convert to BedTool
	# maxchrpos = {d['SN']: d['LN'] for d in bam.header['SQ']}[str(chr)]
	# bin_starts = range(0, maxchrpos - binsize, binsize)
	# bin_stops = range(binsize, maxchrpos, binsize)
	# bins = []
	# for i in range(0, len(bin_starts)-1):
	# 	bins.append([chr, bin_starts[i], bin_stops[i]])
	# bins = pybedtools.BedTool(bins)

	# #Remove bins that have at least 5% overlap with blacklist by size
	# bins_filtered = bins.intersect(blacklist, v=True, f=0.05)

	# #Generate & return coverage
	# coverage = bins_filtered.coverage(subbam, counts=True, sorted=True)
	# return coverage

#Main function
def main():
    #Add arguments
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam', type=pysam.AlignmentFile,
                        help='Input bam')
    parser.add_argument('chr', help='Contig to evaluate')
    parser.add_argument('cov_out', help='Output bed file of raw coverage')
    parser.add_argument('-n', '--norm_out', nargs=1,
    					help='Output bed file of normalized coverage')
    parser.add_argument('-b', '--binsize', type=int, default=1000,
                        help='Bin size in bp (default: 1000)')
    parser.add_argument('-m', '--mode', default='nucleotide',
    					choices = ['nucleotide', 'physical'],
                        help='Evaluate nucleotide or physical coverage '
                             '(default: nucleotide)')
    parser.add_argument('-x', '--blacklist', nargs=1,
    	                help='BED file of regions to ignore')
    parser.add_argument('-v', '--overlap', nargs=1, type=float, default=0.05,
    	   				help='Maximum tolerated blacklist overlap before '
    	   				      'excluding bin')
    args = parser.parse_args()

    #Open outfiles
    fcovout = open(args.cov_out, 'w')
    if args.norm_out is not None:
	    fnormout = open(args.norm_out, 'w')

    #Get nucleotide coverage
    if args.mode == 'nucleotide':
    	coverage = nuc_binCov(args.bam, args.chr, args.binsize, args.blacklist)

    #Close outfiles
    fcovout.close()
    if args.norm_out is not None:
	    fnormout.close()


#Main block
if __name__ == '__main__':
    main()

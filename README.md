# WGD
**W**hole-**G**enome **D**osage: a suite of tools to evaluate dosage in whole-genome sequencing libraries

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.

## binCov.py
Iterates through a single chromosome of a bam file and calculates either nucleotide or physical coverage in regularly segmented bins.
```
usage: binCov.py [-h] [-n NORM_OUT] [-b BINSIZE] [-t {nucleotide,physical}]
                 [-x BLACKLIST] [-v OVERLAP]
                 ibam chr cov_out

Calculates non-duplicate primary-aligned binned coverage of a chromosome from
an input BAM file

positional arguments:
  ibam                  Input bam
  chr                   Contig to evaluate
  cov_out               Output bed file of raw coverage

optional arguments:
  -h, --help            show this help message and exit
  -n NORM_OUT, --norm_out NORM_OUT
                        Output bed file of normalized coverage
  -b BINSIZE, --binsize BINSIZE
                        Bin size in bp (default: 1000)
  -t {nucleotide,physical}, --type {nucleotide,physical}
                        Evaluate nucleotide or physical coverage (default:
                        nucleotide)
  -x BLACKLIST, --blacklist BLACKLIST
                        BED file of regions to ignore
  -v OVERLAP, --overlap OVERLAP
                        Maximum tolerated blacklist overlap before excluding
                        bin
```
**Usage Notes:**  
1. Input bam file must be coordinate-sorted and indexed.  
2. Normalized coverage is raw coverage per bin divided by median of all non-zero, non-blacklisted bins on the same contig
3. Bins will be ignored automatically if they share __any__ overlap with blacklisted regions (-x or --blacklist) 

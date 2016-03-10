# WGD
**W**hole-**G**enome **D**osage: a suite of tools to evaluate dosage in whole-genome sequencing libraries

**Contact:** Ryan Collins (rcollins@chgr.mgh.harvard.edu)

All code copyright (c) 2016 Ryan Collins and is distributed under terms of the MIT license.  
---  
### Example workflow  
#### Prerequisites  
The WGD pipeline requires the following:  
- Coordinate-sorted, indexed bams for all samples
- List of contigs to evaluate (only primary contigs recommended; e.g. 1...22, X, and Y for human)
- Bed-file of N-masked regions of the genome. These are available through the [UCSC Genome Browser](http://genome.ucsc.edu/ "UCSC Genome Browser")  
#### Step 1: Generate normalized coverage per chromosome on all libraries  
Normalized coverage is calculated by ```binCov.py``` on a per-chromosome basis. Parallelization of this process is strongly encouraged.  

There are two approaches to parallelization, depending on your available computational resources. Examples are given below using LSF as a scheduler, but could be easily configured to your scheduler/environment.  

**Parallel submission of all chromosomes from all samples**
```
#!/bin/bash  
while read sample; do
  while read contig; do
    bsub "binCov.py "
  done < list_of_contigs.txt
done < list_of_samples.txt
```

---
### binCov.py
Iterates through a single chromosome of a bam file and calculates either nucleotide or physical coverage in regularly segmented bins.
```
usage: binCov.py [-h] [-n NORM_OUT] [-b BINSIZE] [-m {nucleotide,physical}]
                 [-x BLACKLIST] [-v OVERLAP]
                 bam chr cov_out

Calculates non-duplicate primary-aligned binned coverage of a chromosome from
an input BAM file

positional arguments:
  bam                   Input bam
  chr                   Contig to evaluate
  cov_out               Output bed file of raw coverage

optional arguments:
  -h, --help            show this help message and exit
  -n NORM_OUT, --norm_out NORM_OUT
                        Output bed file of normalized coverage
  -b BINSIZE, --binsize BINSIZE
                        Bin size in bp (default: 1000)
  -m {nucleotide,physical}, --mode {nucleotide,physical}
                        Evaluate nucleotide or physical coverage (default:
                        nucleotide)
  -x BLACKLIST, --blacklist BLACKLIST
                        BED file of regions to ignore
  -v OVERLAP, --overlap OVERLAP
                        Maximum tolerated blacklist overlap before excluding
                        bin
```
**Usage Notes:**  
- Input bam file must be coordinate-sorted and indexed.  
- Only non-duplicate primary-aligned reads or proper pairs are considered for 'nucleotide' and 'physical' mode, respectively.  
- Normalized coverage is raw coverage per bin divided by median of all non-zero, non-blacklisted bins on the same contig.  
- Bins will be ignored automatically if they share at least ```-v``` percent overlap by size with blacklisted regions (```-x``` or ```--blacklist```).  
- Currently uses ```bedtools coverage``` syntax assuming ```bedtools``` version pre-2.24.0 (i.e. ```-a``` is features and ```-b``` is intervals for which to calculate coverage; this was reversed starting in ```bedtools v2.24.0```)

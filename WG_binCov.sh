#!/bin/bash

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

#Usage statement
usage(){
cat <<EOF

usage: WG_binCov.sh [-h] [-b BINSIZE] [-m MODE] 
                    [-L CONTIGS] [-x BLACKLIST] [-v OVERLAP] 
                    BAM ID OUTDIR

Wrapper for serialized execution of binCov.py across multiple chromosomes

Positional arguments:
  BAM     Input bam
  ID      Sample ID
  OUTDIR  Output directory

Optional arguments:
  -h  HELP         Show this help message and exit
  -b  BINSIZE      Bin size in bp (default: 1000)
  -m  MODE         Evaluate physical or nucleotide coverage (default: nucleotide)
  -L  CONTIGS      List of contigs to evaluate (default: all contigs in bam header)
  -x  BLACKLIST    BED file of regions to ignore
  -v  OVERLAP      Maximum tolerated blacklist overlap before excluding bin

EOF
}


#Read arguments
nperm=10000
while getopts ":n:h" opt; do
	case "$opt" in
		n)
			nperm=$OPTARG
			;;
		h)
			usage
			exit 0
			;;
	esac
done
shift $(( OPTIND - 1))
l1=$1 #first gene list
l2=$2 #second gene list
ref=$3 #reference gene list

#Check positional arguments
if [ -z ${l1} ] || [ -z ${l2} ] || [ -z ${ref} ]; then
	usage
	exit 0
fi


bam=$1
ID=$2
OUTDIR=$3

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

#Parse arguments
binsize=1000
mode=nucleotide
contigs=DEFAULT
blist=NONE
v=0.05
while getopts ":b:m:L:x:v:h" opt; do
	case "$opt" in
		b)
			binsize=${OPTARG}
			;;
		m)
			mode=${OPTARG}
			;;
		L)
			contigs=${OPTARG}
			;;
		x)
			blist=${OPTARG}
			;;
		v)
			v=${OPTARG}
			;;
		h)
			usage
			exit 0
			;;
	esac
done
shift $(( ${OPTIND} - 1))
bam=$1
ID=$2
OUTDIR=$3

#Check positional arguments
if [ -z ${bam} ] || [ -z ${ID} ] || [ -z ${OUTDIR} ]; then
	usage
	exit 0
fi

#Determine list of contigs to use (note: requires samtools)
if [ ${contigs} == "DEFAULT" ]; then
	contigs_list=mktemp
	samtools view -H ${bam} | fgrep -w "@SQ" | awk '{ print $2 }' | cut -d\: -f2 > ${contigs_list}
else
	contigs_list=${contigs}
fi

#Run binCov.py on all contigs
spath=$( readlink -f $0 )
echo ${spath}
#while read contig; do





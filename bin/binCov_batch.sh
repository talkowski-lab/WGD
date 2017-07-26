#!/bin/bash

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Batch script to iterate binCov over a supplied list of contigs

#Usage statement
usage(){
cat <<EOF

usage: binCov_batch.sh [-h] [-b BINSIZE] [-m MODE] [-n] 
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
  -n  NORMALIZED   Also generate normalized coverage values 
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
norm=0
while getopts ":b:m:nL:x:v:h" opt; do
	case "$opt" in
		b)
			binsize=${OPTARG}
			;;
		m)
			mode=${OPTARG}
			;;
    n)
      norm=1
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
	samtools view -H ${bam} | fgrep -w "@SQ" | \
  awk '{ print $2 }' | cut -d\: -f2 > \
  ${contigs_list}
else
	contigs_list=${contigs}
fi

#Run binCov.py on all contigs
spath=$( dirname $( readlink -f $0 ) )
echo ${spath}
while read contig; do
	if [ ${blist} != "NONE" ]; then
    if [ ${norm} == 1 ]; then
  		${spath}/binCov.py \
      -n ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.normCov.bed \
      -b ${binsize} \
      -m ${mode} \
      -x ${blist} \
      -v ${v} \
      ${bam} \
      ${contig} \
      ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.rawCov.bed
    else
      ${spath}/binCov.py \
      -b ${binsize} \
      -m ${mode} \
      -x ${blist} \
      -v ${v} \
      ${bam} \
      ${contig} \
      ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.rawCov.bed
    fi
	else
    if [ ${norm} == 1 ]; then
      ${spath}/binCov.py \
      -n ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.normCov.bed \
      -b ${binsize} \
      -m ${mode} \
      -v ${v} \
      ${bam} \
      ${contig} \
      ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.rawCov.bed
    else
      ${spath}/binCov.py \
      -b ${binsize} \
      -m ${mode} \
      -v ${v} \
      ${bam} \
      ${contig} \
      ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.rawCov.bed
    fi
	fi
done < ${contigs_list}

#!/bin/bash

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Wrapper for bedtools unionbedg to create a multi-sample coverage
# matrix from the output of binCov.py

#Usage statement
usage(){
cat <<EOF

usage: makeMatrix.sh [-h] [-o OUTFILE] SAMPLES

Helper tool to automate creation of sorted coverage matrices from
binCov.py output bed files

Positional arguments:
  SAMPLES     List of samples and coverage files (tab-delimmed)

Optional arguments:
  -h  HELP      Show this help message and exit
  -o  OUTFILE   Output file (default: stdout)

EOF
}

#Parse arguments
OUTFILE=/dev/stdout
while getopts ":o:h" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    o)
      OUTFILE=${OPTARG}
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
SAMPLES=$1

#Check for required input
if [ -z ${SAMPLES} ] || ! [ -e ${SAMPLES} ]; then
  usage
  exit 0
fi

#Use bedtools to create matrix
/usr/bin/env bedtools unionbedg \
  -header \
  -names $( while read ID cov; do echo "${ID}"; done < ${SAMPLES} | paste -s -d\  )\
  -i $( while read ID cov; do echo "${cov}"; done < ${SAMPLES} | paste -s -d\  )|\
  sed 's/chrom/chr/g' > ${OUTFILE}
#!/bin/bash

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

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
OUTFILE = /dev/stdout
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

#Iterate over samples, sort, and write to tmpfile
WRKDIR=$( mktemp -d )
while read ID cov; do
  sort -Vk1,1 -k2,2n ${cov} | cut -f4 > ${WRKDIR}/${ID}_covVals.txt
done < ${SAMPLES}

#

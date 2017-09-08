#!/bin/bash

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Wrapper script to run cn.MOPS workflow from a list of binCov input files

#Usage statement
usage(){
cat <<EOF

usage: cnMOPS_workflow.sh [-h] [-r REBIN] [-o OUTDIR] [-c CLEANUP] BINCOVS

Wrapper script to run cn.MOPS workflow from a list of binCov input files

Positional arguments:
  BINCOVS   list of sample IDs and full paths to binCov files, tab-delimmed

Optional arguments:
  -h  HELP      Show this help message and exit
  -r  REBIN     Bin compression ratio (default: 1)
  -o  OUTDIR    Output directory (default: pwd)
  -p  PREFIX    Name attached to all files (default: cnMOPS)
  -c  CLEANUP   Automatically delete all intermediate 
                  files (default: keep all intermediate files)

EOF
}

#Parse arguments
REBIN=1
OUTDIR=`pwd`
PREFIX="cnMOPS"
CLEANUP=0
while getopts ":r:o:p:ch" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    r)
      REBIN=${OPTARG}
      ;;
    o)
      OUTDIR=${OPTARG}
      ;;
    p)
      PREFIX=${OPTARG}
      ;;
    c)
      CLEANUP=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
BINCOVS=$1

#Get path to WGD bin
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

##DEV PARAMS ON BROAD CLUSTER
# REBIN=10
# OUTDIR=${TMPDIR}/cnMOPS_workflow_test
# PREFIX="cnMOPS_test"
# CLEANUP=0
# BINCOVS=${WRKDIR}/analysis/cnMOPS_benchmarking/matrix/XX_plus_20/XX_plus_20.makeMatrix_input.txt
# BIN=${WRKDIR}/bin/WGD/bin/

#Check for required input
if [ -z ${BINCOVS} ]; then
  usage
  exit 0
fi

#Attempts to create OUTDIR and directory tree
for DIR in ${OUTDIR} ${OUTDIR}/calls; do
  if ! [ -e ${DIR} ]; then
    mkdir ${DIR}
  fi
done

#Creates binCov matrix
${BIN}/makeMatrix.sh -z \
-o ${OUTDIR}/${PREFIX}.raw_matrix.bed \
${BINCOVS}

#Recompresses binCov matrix (if optioned)
if [ ${REBIN} -gt 1 ]; then
  ${BIN}/compressCov.sh -z -s \
  -o ${OUTDIR}/${PREFIX}.compressed_matrix.bed \
  ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
  ${REBIN}
else
  mv ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
  ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz
fi

#Run cn.MOPS
${BIN}/runcnMOPS.R \
-I ${PREFIX} \
${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz \
${OUTDIR}/calls/

#Format cn.MOPS calls
${WRKDIR}/bin/WGD/bin/cleancnMOPS.sh -z \
-o ${WRKDIR}/analysis/cnMOPS_benchmarking/calls/ \
<( cut -f1 ${BINCOVS} ) \
<( echo -e "${OUTDIR}/calls/${PREFIX}.cnMOPS.gff" )

#Clean up (if optioned)
if [ ${CLEANUP} -eq 1 ]; then
  for file in ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
              ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz; do
    rm ${file}
  done
fi


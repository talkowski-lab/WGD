#!/bin/bash

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Wrapper script to run cn.MOPS workflow from a list of binCov input files

#Usage statement
usage(){
cat <<EOF

usage: cnMOPS_workflow.sh [-h] [-r REBIN] [-o OUTDIR] [-c CLEANUP] [-x EXCLUDE] [-S SUBTRACT] BINCOVS

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
  -x  BLACKLIST Intervals to be hard-filtered from binCov matrix
  -S  SUBTRACT  Intervals to subtract from any overlapping calls 
                  (e.g. N-masked reference gaps) 

EOF
}

#Parse arguments
REBIN=1
OUTDIR=`pwd`
PREFIX="cnMOPS"
CLEANUP=0
BLACKLIST=0
SUBTRACT=0
while getopts ":r:o:p:cx:S:h" opt; do
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
    x)
      BLACKLIST=${OPTARG}
      ;;
    S)
      SUBTRACT=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
BINCOVS=$1

#Get path to WGD bin
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

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
echo -e "STATUS | $( date +"%T (%m-%d-%y)" ) | GENERATING COVERAGE MATRIX..."
if [ ${BLACKLIST} != "0" ]; then
  ${BIN}/makeMatrix.sh -z \
  -r ${BLACKLIST} \
  -o ${OUTDIR}/${PREFIX}.raw_matrix.bed \
  ${BINCOVS}
else
  ${BIN}/makeMatrix.sh -z \
  -o ${OUTDIR}/${PREFIX}.raw_matrix.bed \
  ${BINCOVS}
fi

#Recompresses binCov matrix (if optioned)
if [ ${REBIN} -gt 1 ]; then
  echo -e "STATUS | $( date +"%T (%m-%d-%y)" ) | REBINNING COVERAGE MATRIX..."
  ${BIN}/compressCov.sh -z -s \
  -o ${OUTDIR}/${PREFIX}.compressed_matrix.bed \
  ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
  ${REBIN}
else
  mv ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
  ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz
fi

#Blacklist binCov matrix (if optioned)
if [ ${BLACKLIST} != "0" ]; then
  bedtools intersect -header -v -wa \
  -a ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz \
  -b ${BLACKLIST} > \
  ${OUTDIR}/${PREFIX}.compressed_matrix.bed2
  mv ${OUTDIR}/${PREFIX}.compressed_matrix.bed2 \
  ${OUTDIR}/${PREFIX}.compressed_matrix.bed
  bgzip -f ${OUTDIR}/${PREFIX}.compressed_matrix.bed
  tabix -f ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz
fi

#Run cn.MOPS
echo -e "STATUS | $( date +"%T (%m-%d-%y)" ) | RUNNING cn.MOPS..."
${BIN}/runcnMOPS.R \
-I ${PREFIX} \
${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz \
${OUTDIR}/calls/

#Write lists of samplesa and input GFFs
cut -f1 ${BINCOVS} > ${OUTDIR}/samples.list
echo -e "${OUTDIR}/calls/${PREFIX}.cnMOPS.gff" > ${OUTDIR}/GFFs.list

#Format cn.MOPS calls
# echo -e "STATUS | $( date +"%T (%m-%d-%y)" ) | FORMATTING cn.MOPS CALLS..."
# if [ ${SUBTRACT} != "0" ]; then
  # ${BIN}/cleancnMOPS.sh -z \
  # -o ${OUTDIR}/calls/ \
  # -S ${SUBTRACT} \
  # ${OUTDIR}/samples.list \
  # ${OUTDIR}/GFFs.list
# else
  # ${BIN}/cleancnMOPS.sh -z \
  # -o ${OUTDIR}/calls/ \
  # ${OUTDIR}/samples.list \
  # ${OUTDIR}/GFFs.list
# fi  

#Clean up (if optioned)
if [ ${CLEANUP} -eq 1 ]; then
  for file in ${OUTDIR}/${PREFIX}.raw_matrix.bed.gz \
              ${OUTDIR}/${PREFIX}.compressed_matrix.bed.gz \
              ${OUTDIR}/samples.list \
              ${OUTDIR}/GFFs.list; do
    rm ${file}
  done
fi


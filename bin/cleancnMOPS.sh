#!/bin/bash

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Simple script to split cn.MOPS calls by sample and copy state and merge
# across multiple runs (e.g. at different resolutions)

#Usage statement
usage(){
cat <<EOF

usage: cleancnMOPS.sh [-h] [-z] [-o OUTDIR] SAMPLES GFFS

Helper tool to split cn.MOPS calls by sample and copy state and merge
across multiple runs (e.g. at different resolutions)

Positional arguments:
  SAMPLES   list of sample IDs
  GFFS      list of full paths to all cn.MOPS outputs to be considered

Optional arguments:
  -h  HELP        Show this help message and exit
  -z  GZIP        Gzip output file
  -o  OUTDIR      Output directory

EOF
}

#Parse arguments
OUTDIR=`pwd`
GZ=0
while getopts ":o:zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    o)
      OUTDIR=${OPTARG}
      ;;
    z)
      GZ=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
SAMPLES=$1
GFFS=$2

#Check for required input
if [ -z ${SAMPLES} ] || [ -z ${GFFS} ]; then
  usage
  exit 0
fi

#Check for gzip optioned only if output file specified
if [ ${OUTDIR} == "/dev/stdout" ] && [ ${GZ} == 1 ]; then
  echo -e "\nOUTDIR required for zip-compressed output"
  usage
  exit 0
fi

#Scrub ".gz" from output filename if provided by user
if [ ${GZ} == 1 ] && [ ${OUTDIR: -3} == ".gz" ]; then
  OUTDIR=$( echo "${OUTDIR}" | sed 's/\.gz//g' )
fi

#Unzip input file if gzipped
GZI=0
if [ $( file ${INPUT} | fgrep "gzip compressed" | wc -l ) -gt 0 ]; then
  GZI=1
  TMPI=`mktemp`; mv ${TMPI} ${TMPI}.gz; TMPI=${TMPI}.gz
  cp ${INPUT} ${TMPI}
  gunzip ${TMPI}
  INPUT=$( echo "${TMPI}" | sed 's/\.gz/\t/g' | cut -f1 )
fi

#Gather number of samples present in file
NSAMP=$( head -n1 ${INPUT} | awk '{ print NF-3 }' )
COLS=$( seq 4 $((${NSAMP}+3)) | paste -s -d, )

#Calculate new bin size
OBIN=$( fgrep -v "#" ${INPUT} | head -n1 | awk '{ print $3-$2 }' )
NBIN=$((${RATIO}*${OBIN}))

#Print header from original file to OUTDIR
head -n1 ${INPUT} > ${OUTDIR}

#Iterate over contigs present in input file
while read CONTIG; do
  #Get min and max coordinate present in input file for contig
  MIN=$( awk -v CONTIG=${CONTIG} '{ if ($1==CONTIG) print $2 }' ${INPUT} | head -n1 )
  MAX=$( awk -v CONTIG=${CONTIG} '{ if ($1==CONTIG) print $3 }' ${INPUT} | tail -n1 )
  #Perform bedtools map function (operation dependent on norm/non-norm)
  bedtools map -c ${COLS} -o ${BEDOP} \
  -a <( paste <( seq ${MIN} ${NBIN} ${MAX} ) \
              <( seq $((${MIN}+${NBIN})) ${NBIN} $((${MAX}+${NBIN})) ) | \
        awk -v CONTIG=${CONTIG} -v OFS="\t" '{ print CONTIG, $1, $2 }' ) \
 -b <( awk -v CONTIG=${CONTIG} -v OFS="\t" '{ if ($1==CONTIG) print $0 }' \
  ${INPUT} ) | awk '{ if ($4!=".") print $0 }'
done < <( fgrep -v "#" ${INPUT} | cut -f1 | sort -Vk1,1 | uniq ) >> ${OUTDIR}

#Gzip output if optioned
if [ ${GZ} == 1 ]; then
  gzip --force ${OUTDIR}
fi

#Clean up
if [ GZI == 1 ]; then
  rm -rf ${TMPI}
  rm -rf ${INPUT}
fi

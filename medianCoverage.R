#!/usr/bin/env Rscript

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Script to calculate median bin coverage per sample for all samples in a bincov
# matrix

# Note: loads entire coverage matrix into memory. This may pose a problem for
# large matrices or on small-memory machines. Other workarounds exist from the
# command line, but most are slower.

# Load library
library(optparse)

# Get command-line arguments
args <- parse_args(OptionParser(usage="%prog covMatrix.bed OUTFILE", 
	                            option_list=list()),
                   positional_arguments=TRUE)

# Checks for appropriate positional arguments
if(length(args$args) != 2) 
  {cat("Incorrect number of required positional arguments\n\n")
   stop()}

# Read matrix
cov <- read.table(args$args[1], header=T)

#Downsample to 1M random rows if nrows > 1M (for computational efficiency)
if(nrow(cov)>1000000){
	cov <- cov[sample(1:nrow(cov), 1000000),]
}

# Get medians with and without zero-cov bins
zerobins <- which(as.integer(apply(cov, 1, median)) == 0)
withzeros <- apply(cov[,-c(1:3)], 2, median)
withoutzeros <- apply(cov[-zerobins,-c(1:3)], 2, median)

# write out results
res <- data.frame("#ID"=names(cov[,-c(1:3)]),
	              "Med_withZeros"=withzeros,
	              "Med_withoutZeros"=withoutzeros)
write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)

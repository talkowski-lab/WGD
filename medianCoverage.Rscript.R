#!/usr/bin/env Rscript

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Script to calculate median bin coverage per sample for all samples in a bincov matrix

# Load library
library(optparse)

# Get command-line arguments
args <- parse_args(OptionParser(usage="%prog covMatrix.bed OUTFILE", option_list=list()),positional_arguments=TRUE)

# checks for appropriate positional arguments
if(length(args$args) != 2) 
  {cat("Incorrect number of required positional arguments\n\n")
   stop()}

# writes args & opts to vars
path.to.matrix <- args$args[1]
OUTFILE <- args$args[2]

# read matrix
cov <- read.table(path.to.matrix,header=T)

# get medians with and without zero-cov bins
zerobins <- which(apply(cov,1,median)==0)
withzeros <- apply(cov[,-c(1:3)],2,median)
withoutzeros <- apply(cov[-zerobins,-c(1:3)],2,median)

# write out results
res <- data.frame("#ID"=names(cov[,-c(1:3)]),
	              "Med_withZeros"=withzeros,
	              "Med_withoutZeros"=withoutzeros)
write.table(res,OUTFILE,sep="\t",col.names=T,row.names=F,quote=F)
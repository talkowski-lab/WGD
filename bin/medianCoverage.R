#!/usr/bin/env Rscript

# Copyright (c) 2016 Ryan Collins <rcollins@chgr.mgh.harvard.edu>
# Distributed under terms of the MIT license.

# Script to calculate median bin coverage per sample for all samples in a bincov
# matrix

# Note: loads entire coverage matrix into memory. This may pose a problem for
# large matrices or on small-memory machines. Other workarounds exist from the
# command line, but most are slower. One suggested alternative is to split
# the input coverage matrices by chromosome prior to computing medians.

# Load library
require(optparse)

# Define options
option_list <- list(
  make_option(c("-b","--binwise"), action="store_true", default=FALSE,
              help="compute medians of all samples per bin [default: median of all bins per sample]"),
  make_option(c("-s","--stdev"), action="store_true", default=FALSE,
              help="compute standard deviation of all bins per sample [default: FALSE]"))

# Get command-line arguments and options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2)
{cat("Incorrect number of required positional arguments\n\n")
  stop()}

# Read matrix
cov <- read.table(args$args[1], header=T)

# Function to compute medians per sample
covPerSample <- function(cov,downsample=1000000,sd=F){
  # Downsample to 1M random rows if nrows > 1M (for computational efficiency)
  if(nrow(cov)>1000000){
    cov <- cov[sample(1:nrow(cov), downsample),]
  }
  # Get medians with and without zero-cov bins
  zerobins <- which(as.integer(apply(cov, 1, median)) == 0)
  withzeros <- as.numeric(apply(cov[,-c(1:3)], 2, median))
  withoutzeros <- as.numeric(apply(cov[-zerobins,-c(1:3)], 2, median))
  #Get SDs with and without zero-cov bins (if optioned)
  if(sd==T){
    withzeros.sd <- as.numeric(apply(cov[,-c(1:3)], 2, sd))
    withoutzeros.sd <- as.numeric(apply(cov[-zerobins,-c(1:3)], 2, sd))
  }
  # compile results df to return
  if(sd==T){
    res <- data.frame("ID"=names(cov[,-c(1:3)]),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros,
                      "SD_withZeros"=withzeros.sd,
                      "SD_withoutZeros"=withoutzeros.sd)
  }else{
    res <- data.frame("ID"=names(cov[,-c(1:3)]),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros)
  }
  # Return output df
  return(res)
}

# Function to compute medians per bin
covPerBin <- function(cov,downsample=500,sd=F){
  # Downsample to 500 random samples if nsamples > 500 (for computational efficiency)
  if(ncol(cov)>503){
    cov <- cov[,sample(1:ncol(cov), downsample)]
  }
  # Get medians with and without zero-cov samples
  meds <- t(apply(cov[,-c(1:3)], 1, function(vals){
    withzeros <- median(vals)
    if(any(vals>0)){
      withoutzeros <- median(vals[which(vals>0)])
    }else{
      withoutzeros <- NA
    }
    return(c(withzeros,withoutzeros))
  }))
  # Get standard deviations (if optioned)
  sds <- t(apply(cov[,-c(1:3)], 1, function(vals){
    withzeros <- sd(vals)
    if(any(vals>0)){
      withoutzeros <- sd(vals[which(vals>0)])
    }else{
      withoutzeros <- NA
    }
    return(c(withzeros,withoutzeros))
  }))
  # compile results df to return
  if(sd==T){
    res <- data.frame("#chr"=cov[,1],"start"=cov[,2],"end"=cov[,3],
                      "Med_withZeros"=meds[,1],
                      "Med_withoutZeros"=meds[,2],
                      "SD_withZeros"=sds[,1],
                      "SD_withoutZeros"=sds[,2])
  }else{
    res <- data.frame("#chr"=cov[,1],"start"=cov[,2],"end"=cov[,3],
                      "Med_withZeros"=meds[,1],
                      "Med_withoutZeros"=meds[,2])
  }
  # Return output df
  return(res)
}

# Compute appropriate medians & write out
if(opts$binwise==TRUE){
  res <- covPerBin(cov,sd=opts$stdev)
  names(res)[1] <- "#chr"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}else{
  res <- covPerSample(cov,sd=opts$stdev)
  names(res)[1] <- "#ID"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}

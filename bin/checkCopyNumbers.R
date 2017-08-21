#!/usr/bin/env Rscript

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Script to analyze a matrix of median binCov values per contig per sample
# for approximate copy numbers (with visualization)

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)

###################################################
#####Helper function to load median coverage matrix
###################################################
readMatrix <- function(INFILE){
  dat <- read.table(INFILE,comment.char="",header=T)
  colnames(dat)[1] <- "ID"
  dat[,-1] <- t(apply(dat[,-1],1,as.numeric))
  return(dat)
}

#############################################################
#####Helper function to normalize contigs for a single sample
#############################################################
normalizeContigsPerSample <- function(vals,exclude){
  #Convert vals to numeric
  vals <- as.numeric(vals)

  #Iterate over values
  newVals <- sapply(1:length(vals),function(i){
    #Compute mean of values excluding current value & any other specified
    excl.mean <- mean(vals[-c(i,exclude)])

    #Normalize current value versus appropriate mean
    return(vals[i]/excl.mean)
  })

  #Return normalized values
  return(newVals)
}

#########################################################################
#####Helper function to normalize contigs for an entire matrix of samples
#########################################################################
normalizeContigsPerMatrix <- function(dat,exclude,ploidy=2){
  #Iterate over samples & normalize
  dat[,-1] <- ploidy*t(apply(dat[,-1],1,normalizeContigsPerSample,exclude=23:24))

  #Return transformed data
  return(dat)
}

#############################################################################
#####Helper function to recenter distributions per chromosome once normalized
#############################################################################
recenterContigDists <- function(dat,exclude){
  #Gather summary statistics
}


# ADD: FUNCTION TO ASSIGN SEXES (CLUSTERING)

###############################################################
#####Helper function to plot distribution of samples per contig
###############################################################
boxplotsPerContig <- function(dat,contigLabels=paste("chr",c(1:22,"X","Y"),sep="")){
  #Load required library
  require(vioplot)

  #Prepare plot area
  plot(x=c(0,ncol(dat)-1),y=c(0,max(4,max(dat[,-1],na.rm=T))),
       type="n",yaxt="n",xaxt="n",yaxs="i",xaxs="i",ylab="",xlab="")

  #Iterate over contigs and plot violins & jitters
  sapply(1:(ncol(dat)-1),function(i){
    #Jitter
    points(x=jitter(rep(i-0.5,times=nrow(dat)),amount=0.3),y=dat[,i+1],
           pch=19,col="#80a1d6",cex=0.5)

  })
}

##########################
#####Rscript functionality
##########################
require(optparse)
#List of command-line options
option_list <- list(
  make_option(c("-O", "--OUTDIR"),type="character",default=NULL,
              help="output directory [default: pwd]",
              metavar="character"),
  make_option(c("-p", "--noplot"),action="store_true",default=FALSE,
              help="disable copy number visualization [default: FALSE]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] median_coverage_matrix",
                                option_list=option_list),
                   positional_arguments=TRUE)
INFILE <- args$args[1]
OUTDIR <- args$options$OUTDIR
noplot <- args$options$noplot
if(is.null(OUTDIR)){
  OUTDIR <- ""
}

#Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Must supply an input median coverage matrix\n")
}

#Loads data
dat <- readMatrix(INFILE)

#Transforms data to predicted copy numbers
dat <- normalizeContigsPerMatrix(dat,exclude=23:24,ploidy=2)










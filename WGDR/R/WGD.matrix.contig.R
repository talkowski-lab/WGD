# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.plotcontig
################
# Plots smoothed dosage for a specified contig from 1-2 groups of samples
# in a WGD matrix. Plot can be underlaid with annotation tracks if specified.
# Note: designed for direct integration into WGD.matrix.plot, but can also
# be called independently
################

WGD.matrix.contig <- function(mat,                       #matrix object from which to plot. Must be read with WGD.readmatrix
                              assignments,               #vector of sample assignments. Can be generated with WGD.matrix.fingerprint
                              contig,                    #contig to plot
                              OUTDIR,                    #output directory
                              grouplabs=c("A","B"),      #labels for each group
                              groupcols=c("blue","red"), #colors for each group
                              filename=NULL,             #add custom filename for output
                              fraction=0.25              #plot top N% of each groups
){
  #Validate faction parameter
  if(is.numeric(fraction)==F | fraction<=0 | fraction>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'fraction' must be float ~ (0,1]",sep=""))
  }

  #Validate assignments vector
  if(!(is.vector(assignments)) | length(assignments)!=(ncol(mat$log2)-3)){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'assignments' be a numeric vector with one entry for each sample in 'mat'",sep=""))
  }

  #Ensure output directory exists
  if(file.exists(OUTDIR)==F){
    cat(paste("WGDR::WARNING [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "] ",match.call()[[1]]," specified output directory (",
              OUTDIR,") does not exist. Attempting to create...\n",
              sep=""))
    dir.create(OUTDIR,recursive=T)
  }

  #Disable scientific notation
  options(scipen=1000)

  #Determine number of groups
  groups <- length(unique(assignments))

  #Determine top fraction of each group
  sapply(groups,function(g){
    grp <- WGD.matrix.postprocess(mat$mat[,c(1:3,which(assignments==g)+3)])

  })









}

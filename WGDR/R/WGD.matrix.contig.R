# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.contig
################
# Plots smoothed dosage for a specified contig from one or more groups of samples
# in WGD matrices. Plots can represent all samples or distribution of samples.
# Plot can be underlaid with annotation tracks if specified. Note: designed
# for direct integration into WGD.matrix.plot, but can also be called independently
################

WGD.matrix.contig <- function(mats,                      #list of matrixes from which to plot. Must be read with WGD.readmatrix
                              contig,                    #contig to plot
                              OUTDIR,                    #output directory
                              coords=NULL,               #specify limits on coordinates as vector
                              covlims=c(-0.25,0.25),     #specify limits on coverage y-axis
                              mode=c("d","a","e"),       #plot mode: d=distribution; a=all; e=extremes
                              smoothing=0.25,            #smoothing factor ~ (0,1]
                              fraction=0.25,             #fraction of samples to plot (note: no effect whe mode=d)
                              grouplabs=c("A","B"),      #labels for each group
                              groupcols=c("blue","red"), #colors for each group
                              filename=NULL              #add custom filename for output
){
  #Validate faction parameter
  if(is.numeric(fraction)==F | fraction<=0 | fraction>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'fraction' must be float ~ (0,1]",sep=""))
  }

  #Validate smoothing parameter
  if(is.numeric(smoothing)==F | smoothing<=0 | smoothing>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'smoothing' must be float ~ (0,1]",sep=""))
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

  #Determine plotting coordinates
  if(!is.null(coords))
  maxCoord <- max(unlist(lapply(mats,function(mat){
    return(max(mat$mat$end))
  })))
  coords <- c(0,maxCoord)

  #Prepare plot









}

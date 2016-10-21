# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.tranche
################
# Subsets WGD matrix to retain middle X% of bins
################
# Returns an eight-item list:
#  $mat : matrix of original values
#  $res : matrix of residuals
#  $log2 : matrix of log2-fold changes
#  $stat : matrix of per-bin distribution statistics
#  $rstat : matrix of per-bin residual distribution statistics
#  $lstat : matrix of per-bin log2-fold change distribution statistics
#  $sstat.res : matrix of per-sample residual distribution statistics
#  $sstat.log2 : matrix of per-sample residual distribution statistics
################

WGD.matrix.tranche <- function(mat,            #matrix object from which to read. Must be read with WGD.readmatrix
                               tranche=0.99,   #fraction of bins to keep
                               quiet=F         #option to disable verbose output
){
  #Sanity check tranche parameter
  if(is.numeric(tranche)==F | tranche<=0 | tranche>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'tranche' must be float ~ (0,1]",sep=""))
  }

  #Disable scientific notation
  options(scipen=1000)

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: applying tranche partition (tranche=",tranche,")\n",sep=""))
  }

  #Determines min and max cutoffs
  Q1.min <- quantile(mat$stat$Q1,probs=(1-tranche)/2)
  Q3.max <- quantile(mat$stat$Q3,probs=1-((1-tranche)/2))

  #Filters bins
  pass <- which(mat$stat$Q1>Q1.min
                & mat$stat$Q3<Q3.max)
  nmat <- WGD.matrix.postprocess(mat$mat[pass,])

  #Returns filtered matrix
  return(nmat)
}

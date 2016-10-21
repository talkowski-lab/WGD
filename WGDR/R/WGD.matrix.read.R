# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.read
################
# Helper function for standardized import of WGD normalized coverage
# matrix and automated summary statistic collection.
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

WGD.matrix.read <- function(path,        #full path to matrix file
                           allosomes=F,  #option to auto-exclude non-numeric contigs
                           norm=F,       #option to normalize coverage matrix; only necessary if using raw binCov matrix
                           quiet=F       #option to disable verbose output
){
  #Prohibit scientific notation & auto string factorization
  options(scipen=1000,
          stringsAsFactors=F)

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: processing data from ",path,"\n",sep=""))
  }

  #Read matrix
  mat <- read.table(path,header=T)

  #Process matrix
  return(WGD.matrix.postprocess(mat,
                                norm=norm,
                                allosomes=allosomes))
}

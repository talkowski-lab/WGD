# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.read
################
# Helper function for standardized import of WGD normalized coverage
# matrix and automated summary statistic collection.
################
# Returns a four-item list:
#  $mat : matrix of original values
#  $res : matrix of residuals
#  $stat : matrix of per-bin distribution statistics
#  $rstat : matrix of per-bin residual distribution statistics
#  $sstat : matrix of per-sample residual distribution statistics
################

WGD.matrix.read <- function(path,        #full path to matrix file
                           allosomes=F,  #option to auto-exclude non-numeric contigs
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
  return(WGD.matrix.postprocess(mat))
}

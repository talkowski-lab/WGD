# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.read
################
# Helper function for standardized import of WGD normalized coverage
# matrix and automated summary statistic collection.
################
# Returns an five-item list:
#  $mat : matrix of original values
#  $res : matrix of residuals
#  $stat : matrix of per-bin distribution statistics
#  $rstat : matrix of per-bin residual distribution statistics
#  $sstat.res : matrix of per-sample residual distribution statistics
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
              "]: importing data from ",path,"\n",sep=""))
  }

  #Read matrix & relabel bed column headers
  mat <- read.table(path,header=T,comment.char="")
  names(mat)[1:3] <- c("chr","start","end")

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: imported ",ncol(mat)-3," samples\n",sep=""))
  }

  #Process matrix
  return(WGD.matrix.postprocess(mat,
                                norm=norm,
                                allosomes=allosomes,
                                quiet=quiet))
}

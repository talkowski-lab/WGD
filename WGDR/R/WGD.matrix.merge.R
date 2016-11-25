# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.merge
################
# Helper function for standardized merger of multiple WGD normalized
# coverage matrices and automated summary statistic collection.
################
# Returns an five-item list:
#  $mat : matrix of original values
#  $res : matrix of residuals
#  $stat : matrix of per-bin distribution statistics
#  $rstat : matrix of per-bin residual distribution statistics
#  $sstat.res : matrix of per-sample residual distribution statistics
################

WGD.matrix.merge <- function(matrices,     #list of matrices imported with WGD.matrix.read
                             innerjoin=T,  #option to only retains bins that appear in all matrices
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
              "]: merging data from ",length(matrices)," matrices\n",sep=""))
  }

  #Merge first two matrices
  if(innerjoin==T){
    mat <- merge(matrices[[1]]$mat,
                 matrices[[2]]$mat,
                 sort=F)
  }else{
    mat <- merge(matrices[[1]]$mat,
                 matrices[[2]]$mat,
                 all=T,sort=F)
  }

  #Iteratively merge third:final matrices
  if(length(matrices)>2){
    for(i in 3:length(matrices)){
      if(innerjoin==T){
        mat <- merge(mat,
                     matrices[[i]]$mat)
      }else{
        mat <- merge(mat,
                     matrices[[i]]$mat,
                     all=T)
      }
    }
  }

  #Process new merged matrix & return
  return(WGD.matrix.postprocess(mat,quiet=quiet))

}

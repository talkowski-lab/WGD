# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.slice
################
# Helper function to extract the top Nth-centile of most
# or least biased libraries from a WGD matrix
################
# Returns an five-item list:
#  $mat : matrix of original values
#  $res : matrix of residuals
#  $stat : matrix of per-bin distribution statistics
#  $rstat : matrix of per-bin residual distribution statistics
#  $sstat.res : matrix of per-sample residual distribution statistics
################

WGD.matrix.slice <- function(mat,          #WGD matrix object. Must be read with WGD.matrix.read
                             slice=0.1,    #returns N% most or least variable libraries (based on sd)
                             top=T,        #top=T returns most biased, top=F returns least biased
                             quiet=F       #option to disable verbose output
){
  #Prohibit scientific notation & auto string factorization
  options(scipen=1000,
          stringsAsFactors=F)

  #Validate slice parameter
  if(is.numeric(slice)==F | slice<=0 | slice>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'slice' must be float ~ (0,1]",sep=""))
  }

  #Prints status
  if(top==T){
    direction <- c("top","most","≥")
  }else{
    direction <- c("bottom","least","<")
  }
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: slicing ",direction[1]," ",round(100*slice,2),"% ",
              direction[2]," biased samples\n",sep=""))
  }

  #Determine samples to slice
  if(top==T){
    slice=1-slice
    sdthresh <- quantile(mat$sstat.res$sd,slice)
    keep <- which(mat$sstat.res$sd>=sdthresh)+3
  }else{
    sdthresh <- quantile(mat$sstat.res$sd,slice)
    keep <- which(mat$sstat.res$sd<sdthresh)+3
  }

  #Prints outcome
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: retaining ",length(keep),"/",ncol(mat$mat)-3," samples with binwise sd ",
              direction[3]," ",round(sdthresh,4),"\n",sep=""))
  }

  #Reprocess sliced matrix
  return(WGD.matrix.postprocess(mat$mat[,c(1:3,keep)],
                                quiet=quiet))
}

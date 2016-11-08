# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.fingerprint.read
################
# Reads a WGD fingerprint from file
################
# Returns a three-item list:
#  $fp.full : a BED-style data frame of all bins used in fingerprinting with mean and sd per assignment group
#  $fp.final : a BED-style data frame of the final set of bins used for training the GMM classifier
#  $model : GMM classifier fit to fp.final
################

WGD.fingerprint.read <- function(FPDIR,     #input directory containing all fingerprint files
                                 tag="WGD", #tag (string) to be matched at the front of each file name
                                 quiet=F    #option to disable verbose output
){
  #Require mclust
  require(mclust)

  #Assigns paths
  fp.full.path <- paste(FPDIR,"/",tag,".","fingerprint_full.bed",sep="")
  fp.final.path <- paste(FPDIR,"/",tag,".","fingerprint_final.bed",sep="")
  training.data.path <- paste(FPDIR,"/",tag,".","fingerprint_training_data.txt",sep="")

  #Checks if full fingerprint file exists
  if(!(file.exists(fp.full.path))){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] unable to find full fingerprint file (",
               fp.full.path,")\n",sep=""))
  }

  #Checks if final fingerprint file exists
  if(!(file.exists(fp.final.path))){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] unable to find final fingerprint file (",
               fp.final.path,")\n",sep=""))
  }

  #Checks if fingerprint training data exists
  if(!(file.exists(training.data.path))){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] unable to find fingerprint training data (",
               training.data.path,")\n",sep=""))
  }

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: reading saved fingerprint from file\n",
              "WGDR::STATUS [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "]: please be patient; this may take a while...\n",
               sep=""))
  }

  #Read fingerprint files
  fp.full <- read.table(fp.full.path,header=T)
  fp.final <- read.table(fp.final.path,header=T)

  #Read training data & parse
  training.data <- read.table(training.data.path,header=T)
  assignments <- as.factor(training.data[,ncol(training.data)])
  train.matrix <- training.data[,-ncol(training.data)]

  #Fit GMM
  GMM <- MclustDA(data=train.matrix,class=assignments)

  #Return fingerprint
  return(list("fp.full"=fp.full,
              "fp.final"=fp.final,
              "model"=GMM))
}

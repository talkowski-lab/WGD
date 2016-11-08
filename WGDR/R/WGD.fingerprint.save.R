# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.fingerprint.save
################
# Saves a WGD fingerprint to file for future use
################

WGD.fingerprint.save <- function(fingerprint,  #WGD fingerprint object to save
                                 OUTDIR,       #output directory -- will be created if does not exist
                                 tag="WGD",    #tag (string) to be added to the front of each file name
                                 quiet=F       #option to disable verbose output
){
  #Create output directory if does not exist
  if(!(dir.exists(OUTDIR))){
    #Prints status
    if(quiet==F){
      cat(paste("WGDR::STATUS [",
                strsplit(as.character(Sys.time()),split=" ")[[1]][2],
                "]: creating output directory: ",OUTDIR,"\n",sep=""))
    }
    dir.create(OUTDIR)
  }

  #Use scientific notation for decimals longer than six digits
  options(scipen=6)

  #Write full fingerprint to file
  write.table(fingerprint$fp.full,
              paste(OUTDIR,"/",tag,".","fingerprint_full.bed",sep=""),
              col.names=T,row.names=F,sep="\t",quote=F)

  #Write final fingerprint to file
  write.table(fingerprint$fp.final,
              paste(OUTDIR,"/",tag,".","fingerprint_final.bed",sep=""),
              col.names=T,row.names=F,sep="\t",quote=F)

  #Write training data to file
  training.dat <- cbind(fingerprint$model$data,
                        as.character(fingerprint$model$class))
  colnames(training.dat)[ncol(training.dat)] <- "class"
  write.table(training.dat,
              paste(OUTDIR,"/",tag,".","fingerprint_training_data.txt",sep=""),
              col.names=T,row.names=T,sep="\t",quote=F)
}

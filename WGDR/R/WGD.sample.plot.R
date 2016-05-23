# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.plotsample
################
# Plots genome-wide dosage for a single sample in a WGD matrix.
# Writes as pdf to OUTDIR
# Note: uses many undocumented helper plot functions, see
# package contents for function details.
################

WGD.sample.plot <- function(mat,            #matrix object from which to plot. Must be read with WGD.readmatrix
                            ID,             #ID of sample to plot. Must match column header in mat
                            OUTDIR,         #output directory for plot
                            filename=NULL,  #add custom filename for output
                            sampling=25     #sampling & smoothing rate ([1,100], 1 = least, 100 = most)
){
  #Validate sampling parameter
  if(is.numeric(sampling)==F | sampling<1 | sampling>100){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'sampling' must be integer ~ [1,100]",sep=""))
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

  #Gather identity and order of contigs
  contigs <- unique(mat$mat$chr)

  #Prepare layout
  if(!is.null(filename)){
    outname <- paste(OUTDIR,filename,sep="")
  }else{
    outname <- paste(OUTDIR,ID,".WGD_report.pdf",sep="")
  }
  pdf(outname,
      height=10,width=8,
      title=paste("WGD_REPORT_",ID,sep=""))
  par(oma=c(0.25,0.25,0.25,0.25))
  layout(matrix(c(1,1,1,1,
                  2,3,4,5,
                  6,7,8,9,
                  10,11,12,13,
                  14,15,16,17,
                  18,19,20,21,
                  22,23,25,25,
                  24,24,25,25),
                byrow=T,ncol=4))

  #Plot genome-wide strip at top
  WGD.sample.gstrip(mat,ID,sampling=1)

  #Plot first 22 contigs
  if(length(contigs)>=22){
    for(contig in contigs[1:22]){
      WGD.sample.contig(mat,ID,contig,sampling)
    }
  }else{
    for(contig in contigs){
      WGD.sample.contig(mat,ID,contig,sampling)
    }
    for(i in c((22-length(contigs)):22)){
      par(bty="n")
      plot(0,0,type="n",
           xaxt="n",yaxt="n",xlab="",ylab="")
    }
    par(bty="o")
  }

  #Summary plots
  WGD.sample.summaryplots(mat,ID)

  #Close output
  dev.off()
}


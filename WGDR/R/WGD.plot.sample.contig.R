# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.plot.sample.contig
################
# Plots dosage for a single sample in a WGD matrix.
# Writes as jpeg to OUTDIR
################

WGD.plot.sample.contig <- function(mat,                #matrix object from which to plot. Must be read with WGD.readmatrix
                                   ID,                 #ID of sample to plot. Must match column header in mat
                                   OUTDIR,             #output directory for plot
                                   filename=NULL,      #add custom filename for output (must be .jpg or .jpeg)
                                   roll=11,            #bins used in rolling mean
                                   sampling=5,         #sampling rate for bins to plot
                                   res=150,            #jpeg resolution passed to jpeg() call
                                   ylims=NULL,         #limits on y-axis; defaults to displaying 95% of plot data
                                   color="cyan",       #color to plot sample
                                   bg.track=NA,        #annotation track to plot on right Y-axis
                                   bg.roll=1           #bins used in rolling mean of bg.track
){
  #Validate roll parameter
  if(is.numeric(roll)==F | roll<1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'roll' must be positive, whole integer\n",sep=""))
  }

  #Validate bg.roll parameter
  if(is.numeric(bg.roll)==F | bg.roll<1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'bg.roll' must be positive, whole integer\n",sep=""))
  }

  #Require zoo for rolling means
  require(zoo)

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

  #Slice sample of interest from mat
  pd <- mat$res[,c(1:3,which(colnames(mat$res)==ID))]

  #Apply rolling mean to normalized coverage residuals
  if(roll>1){
    pd[,4] <- rollmean(pd[,4],roll,na.pad=T)
  }

  #Downsample normalized coverage residuals
  pd <- pd[seq(1,nrow(pd),sampling),]

  #Instantiate y-axis limits
  # if(is.null(ylims)){
    ymax <- 2*round(10*max(abs(c(quantile(pd[,4],0.025,na.rm=T),quantile(pd[,4],0.975,na.rm=T)))),0)/10
    ylims <- c(-ymax,ymax)
  # }

  plot(pd[,2],pd[,4],ylim=ylims,col=color,cex=0.2)


  #Prepare layout
  if(!is.null(filename)){
    outname <- paste(OUTDIR,filename,sep="")
  }else{
    outname <- paste(OUTDIR,ID,".WGD_report.jpg",sep="")
  }
  jpeg(outname,
       height=300*11*0.8,width=300*8.5*0.8,quality=100,res=300,
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


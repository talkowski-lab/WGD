# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.sample.gstrip
################
# Plots genome-wide summary strip of dosage for a single sample in a WGD matrix.
# Note: undocumented; designed for direct integration into WGD.plotsample
# See WGD.plotsample documentation for relevant info
################

WGD.sample.gstrip <- function(mat,        #matrix object from which to plot. Must be read with WGD.readmatrix
                             ID,          #ID of sample to plot. Must match column header in mat
                             sampling=1  #sampling & smoothing rate ([1,100], 1 = least, 100 = most)
){
  #Disable scientific notation
  options(scipen=1000)

  #Libraries
  require(grDevices)

  #Get plotting values
  vals <- mat$mat[seq(1,nrow(mat$mat),by=sampling),
                  which(names(mat$mat)==ID)]
  binsize <- median(mat$mat[,3]-mat$mat[,2])

  #Get chromosome breaks
  chrbreaks <- as.data.frame(t(sapply(unique(mat$mat$chr),
                                      function(contig){
                                        minBin <- min(which(mat$mat$chr==contig))
                                        maxBin <- max(which(mat$mat$chr==contig))
                                        return(c(contig,minBin,maxBin))
                                      })))
  names(chrbreaks) <- c("chr","min","max")

  #Set color palette
  colGen <- colorRampPalette(c("red","white","blue"))(100)
  cvals <- round((vals*100)-50,0)
  cvals[which(cvals<1)] <- 1
  cvals[which(cvals>100)] <- 100

  #Plot parameters
  par(mar=c(0.6,4.1,1.6,0.1))

  #Plot
  plot(vals,
       pch=15,col=colGen[cvals],
       xaxt="n",yaxt="n",xlab="",yaxs="i",xaxs="i",
       ylim=c(0,2),
       main=paste("WGD Report: Sample ",ID,sep=""),
       ylab="Norm. Cov. Deviance",
       cex.lab=0.9,cex=0.8,
       panel.first=c(rect(xleft=par("usr")[1],
                          xright=par("usr")[2],
                          ybottom=par("usr")[3],
                          ytop=par("usr")[4],
                          border=NA,col="gray97"),
                     rect(xleft=par("usr")[1],
                          xright=par("usr")[2],
                          ybottom=0.25,
                          ytop=1.75,
                          border=NA,col="gray90"),
                     rect(xleft=par("usr")[1],
                          xright=par("usr")[2],
                          ybottom=0.5,
                          ytop=1.5,
                          border=NA,col="gray80"),
                     rect(xleft=par("usr")[1],
                          xright=par("usr")[2],
                          ybottom=0.75,
                          ytop=1.25,
                          border=NA,col="gray60")))
  abline(h=c(1,0.75,1.25,0.5,1.5,0.25,1.75),
         lty=c(1,rep(2,6)),
         lwd=c(2,rep(1,6)),
         col=c("gray10",rep("gray30",6)))
  abline(v=chrbreaks[,3])
  axis(2,at=seq(0,2,by=0.25),
       labels=c("-100%","","-50%","","0%",
                "","+50%","","+100%"),
       las=2,cex.axis=0.8)
  text(x=(chrbreaks[,3]+chrbreaks[,2])/2,
       y=par("usr")[3],pos=3,
       labels=chrbreaks[,1],
       font=2,cex=0.7)
  mtext(text=paste(prettyNum(binsize,big.mark=",")," bp Bins",sep=""),
        adj=1,cex=0.6,font=3)
}

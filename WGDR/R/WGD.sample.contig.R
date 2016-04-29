# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.sample.plotcontig
################
# Plots smoothed dosage for a specified contig from a single sample
# in a WGD matrix. Plot is backlaid with summary stats of cohort for
# interpretability.
# Note: undocumented; designed for direct integration into WGD.plotsample
# See WGD.plotsample documentation for relevant info
################

WGD.sample.contig <- function(mat,        #matrix object from which to plot. Must be read with WGD.readmatrix
                              ID,         #ID of sample to plot. Must match column header in mat
                              contig,     #contig to plot
                              sampling=25 #sampling & smoothing rate ([1,100], 1 = least, 100 = most)
){
  #Disable scientific notation
  options(scipen=1000)

  #Libraries
  require(grDevices)

  #Get plotting values
  vals <- lowess(mat$mat[which(mat$mat$chr==contig),
                         which(names(mat$mat)==ID)],
                 f=0.1*(sampling/100))

  #Get IQR values
  lowest <- lowess(mat$stat$min[which(mat$mat$chr==contig)],
                   f=0.1*(sampling/100))
  lower <- lowess(mat$stat$Q1[which(mat$mat$chr==contig)],
                  f=0.1*(sampling/100))
  med <- lowess(mat$stat$med[which(mat$mat$chr==contig)],
                  f=0.1*(sampling/100))
  upper <- lowess(mat$stat$Q3[which(mat$mat$chr==contig)],
                  f=0.1*(sampling/100))
  highest <- lowess(mat$stat$max[which(mat$mat$chr==contig)],
                   f=0.1*(sampling/100))

  #Set color palette
  colGen <- colorRampPalette(c("red","white","blue"))(100)
  cvals <- (2*round((vals$y*100)-100,0))+50
  cvals[which(cvals<1)] <- 1
  cvals[which(cvals>100)] <- 100

  #Plot parameters
  par(mar=c(1.1,4.1,2.1,0.1))

  #Plot
  plot(x=range(vals$x),y=c(0.5,1.5),
       type="n",
       main=contig,
       ylab="Norm. Cov. Deviance",
       xaxt="n",yaxt="n",xlab="",xaxs="i",yaxs="i",
       panel.first=c(rect(xleft=par("usr")[1],
                          xright=par("usr")[2],
                          ybottom=par("usr")[3],
                          ytop=par("usr")[4],
                          border=NA,col="gray95")))
  polygon(x=c(lower$x,rev(upper$x)),
          y=c(med$y-(2*(med$y-lower$y)),
              rev(med$y+(2*(upper$y-med$y)))),
          border="navajowhite3",col="moccasin")
  polygon(x=c(lower$x,rev(upper$x)),
          y=c(lower$y,rev(upper$y)),
          border="navajowhite3",col="navajowhite2")
  points(med,type="l",col="white")
  points(vals,type="l",lwd=3)
  abline(h=c(0.75,1,1.25),
         lty=c(2,1,2),lwd=c(1,2,1),
         col=c("gray30","black","gray30"))
  segments(x0=vals$x[1:(length(vals$x)-1)],
           x1=vals$x[2:length(vals$x)],
           y0=vals$y[1:(length(vals$y)-1)],
           y1=vals$y[2:length(vals$y)],
           col=colGen[cvals],lwd=1)
  points(vals$x[seq(1,length(vals$x),by=sampling)],
         vals$y[seq(1,length(vals$x),by=sampling)],
         pch=19,cex=0.8,
         col=colGen[cvals[seq(1,length(vals$x),by=sampling)]])
  points(vals$x[seq(1,length(vals$x),by=sampling)],
         vals$y[seq(1,length(vals$x),by=sampling)],
         pch=21,cex=0.8)
  axis(2,at=seq(0.5,1.5,by=0.25),
       labels=c("-50%","-25%","0%","+25%","+50%"),
       las=2,cex.axis=0.9)
}

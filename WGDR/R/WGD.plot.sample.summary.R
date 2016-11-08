# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.plot.sample.summary
################
# Plots a boxplot of bins per contig and a distribution of all
# normalized coverage genome-wide
# Note: undocumented; designed for direct integration into WGD.plotsample
# See WGD.plotsample documentation for relevant info
################

WGD.plot.sample.summary <- function(mat,  #matrix object from which to plot. Must be read with WGD.readmatrix
                                    ID    #ID of sample to plot. Must match column header in mat
){
  #Disable scientific notation
  options(scipen=1000)

  #Libraries
  require(grDevices)

  #Set color palette
  colGen <- colorRampPalette(c("red","white","blue"))(100)

  #Gather identity and order of contigs
  contigs <- unique(mat$mat$chr)

  #Get sample values
  vals <- mat$mat[,which(names(mat$mat)==ID)]

  #Boxplot per contig
  par(mar=c(3.1,4.1,0.6,0.1))
  boxplot(vals ~ as.factor(mat$mat$chr),
          col=c("powderblue","gold"),outline=F,ylim=c(0,2),lty=1,
          xaxs="i",yaxs="i",yaxt="n",xaxt="n",
          ylab="Norm. Cov. Dev.")
  abline(h=c(0.5,1.5),
         lty=2,
         col=c("red","blue"))
  axis(2,at=seq(0,2,by=0.25),
       labels=c("-100%","","-50%","","0%",
                "","+50%","","+100%"),
       las=2,cex.axis=0.8)
  text(x=seq(1:length(contigs)),
       y=par("usr")[3],
       cex=0.8,pos=3,
       labels=contigs,
       font=2)

  #Histogram of all bins in genome
  par(mar=c(4.1,4.1,2.1,1.1))
  hdat <- hist(vals,
               breaks=seq(round(min(vals)-1,0),round(max(vals)+1,0),by=0.025),
               plot=F)
  histColVals <- 2*(round(50*hdat$mids)-50)+50
  histColVals[which(histColVals<1)] <- 1
  histColVals[which(histColVals>100)] <- 100
  plot(x=c(0,2),y=c(0,1.05*max(hdat$counts)),
       type="n",
       main="Distribution of Normalized Coverage Deviance",
       xlab="Norm. Cov. Deviance",xaxs="i",xaxt="n",
       ylab="Bins (Count)",yaxs="i",yaxt="n")
  abline(v=seq(0.25,1.75,by=0.25),
         lty=2,col=c(rep("red",3),"white",
                     rep("blue",3)))
  rect(xleft=hdat$breaks[1:(length(hdat$breaks)-1)],
       xright=hdat$breaks[2:length(hdat$breaks)],
       ybottom=0,
       ytop=hdat$counts,
       col=colGen[histColVals])
  abline(v=1,lwd=2)
  #Annotations
  axis(1,at=seq(0,2,by=0.25),
       labels=c("-100%","","-50%","","0%",
                "","+50%","","+100%"),
       cex.axis=0.9)
  axis(2,at=seq(0,par("usr")[4],
                by=round(par("usr")[4]/6,-2)),
       labels=prettyNum(seq(0,par("usr")[4],
                            by=round(par("usr")[4]/6,-2)),
                        big.mark=","),
       las=2,cex=0.8)
  rect(xleft=par("usr")[1],
       xright=0.45,
       ybottom=0.45*par("usr")[4],
       ytop=par("usr")[4],
       col="white")
  text(x=par("usr")[1]+0.025,
       y=0.975*par("usr")[4],
       adj=c(0,1),cex=0.8,font=2,
       labels=paste("Skewness: ",
                    round(mat$sstat[which(rownames(mat$sstat)==ID),9],2),"\n",
                    "Kurtosis: ",
                    round(mat$sstat[which(rownames(mat$sstat)==ID),10],2),"\n",
                    "  Mean: ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),4],2),"%\n",
                    "  Median: ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),3],2),"%\n",
                    "  Mode: ",
                    100*round(mean(hdat$mids[which(hdat$counts==max(hdat$counts))])-1,2),"%\n",
                    "StDev: ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),7],2),"%\n",
                    "MAD: ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),8],2),"%\n",
                    "IQR: ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),2],2),"% - ",
                    100*round(mat$sstat[which(rownames(mat$sstat)==ID),5],2),"%",
                    sep=""))
  points(x=c(rep(0.03,3),
             mean(hdat$mids[which(hdat$counts==max(hdat$counts))]),
             1+mat$sstat[which(rownames(mat$sstat)==ID),3],
             1+mat$sstat[which(rownames(mat$sstat)==ID),4]),
         y=c(par("usr")[4]*c(0.71,0.77,0.83,
                             rep(0.99,3))),
         pch=25,col="black",bg=c("aquamarine",
                                 "yellow",
                                 "forestgreen"))
  #Normal Distrib
  nfit <- dnorm(hdat$mids,mean=1,sd=0.2)
  lines(x=hdat$mids,
        y=0.95*nfit*par("usr")[4]/max(nfit),
        lwd=2)
}

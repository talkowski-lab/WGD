# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.fingerprint.create
################
# Isolates top N% of most discordant bins from a pair of WGDR matrices
# and fits a GMM to those bins from two training matrices
################
# Returns a two-item list:
#  $fp.full : a BED-style data frame of all bins used in fingerprinting with mean and sd per assignment group
#  $fp.final : a BED-style data frame of the final set of bins used for training the GMM classifier
#  $model : GMM classifier fit to fp.final
################

WGD.fingerprint.create <- function(matA,matB,   #matrix objects from which to read. Must be read with WGD.readmatrix
                            alpha=1E-8,         #Bonferroni-corrected false discovery rate of selected bins will be controlled to this probability
                            minDiff=0.6,        #minimum difference between means for two-sample t-test
                            maxBins=100,        #restricts maximum number of bins for training
                            QCplot=F,           #option to plot fingerprint QC plots
                            OUTDIR=NA,          #output directory for fingerprint QC plots
                            quiet=F             #option to disable verbose output
){
  #Sanity check alpha parameter
  if(is.numeric(alpha)==F | alpha<=0 | alpha>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'alpha' must be float ~ (0,1]",sep=""))
  }

  #Sanity check minDiff parameter
  if(is.numeric(minDiff)==F | minDiff<=0 | minDiff>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'minDiff' must be float ~ (0,1]",sep=""))
  }

  #Allow scientific notation for decimals with more than 6 digits
  options(scipen=6)

  #Load required packages
  require(mclust)

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: fingerprinting coverage matrices (alpha=",alpha,"; minDiff=",minDiff,")\n",sep=""))
  }

  #Restrict to bins that exist in both matA and matB
  bins.in.both <- merge(matA$mat,matB$mat,sort=F)[,1:3]
  if(nrow(bins.in.both)!=nrow(matA$mat)){
    matA$mat <- merge(matA$mat,bins.in.both,sort=F)
    matA <- WGD.matrix.postprocess(matA$mat,allosomes=T,quiet=T)
  }
  if(nrow(bins.in.both)!=nrow(matB$mat)){
    matB$mat <- merge(matB$mat,bins.in.both,sort=F)
    matB <- WGD.matrix.postprocess(matB$mat,allosomes=T,quiet=T)
  }

  #Isolates bins to test based on minimum absolute difference in means
  dAB <- matA$rstat$mean.rank-matB$rstat$mean.rank
  bins.to.test <- which(abs(dAB)>=minDiff)

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: ",prettyNum(length(bins.to.test),big.mark=","),
              "/",prettyNum(nrow(bins.in.both),big.mark=","),
              " bins will be tested as candidates\n",sep=""))
  }

  #Two-sample t-test between groups for all bins to select for significantly different bins
  bins.t <- sapply(bins.to.test,function(i){
    if(matA$stat$mean.rank[i] >= matB$stat$mean.rank[i]){
      return(t.test(x=matA$res.ranks[i,-c(1:3)],
                    y=matB$res.ranks[i,-c(1:3)],
                    mu=minDiff,alternative="greater")$p.value)
    }else{
      return(t.test(x=matB$res.ranks[i,-c(1:3)],
                    y=matA$res.ranks[i,-c(1:3)],
                    mu=minDiff,alternative="greater")$p.value)
    }
  })
  bins.t <- p.adjust(bins.t,method="bonferroni")

  #Identifies fingerprinting bins
  bins.to.keep <- bins.to.test[which(bins.t<=alpha)]
  fingerprint <- bins.in.both[bins.to.keep,]

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: kept ",prettyNum(length(bins.to.keep),big.mark=","),
              "/",prettyNum(length(bins.to.test),big.mark=","),
              " candidate bins\n",sep=""))
  }

  #Merges matA and matB to determine training distributions
  all <- WGD.matrix.merge(list(matA,matB),quiet=T)

  #Add priors to fingerprint bins
  fingerprint$dAB <- dAB[bins.to.keep]
  fingerprint$p <- bins.t[which(bins.t<=alpha)]
  fingerprint$mean <- all$rstat$mean[bins.to.keep]
  fingerprint$sd <- all$rstat$sd[bins.to.keep]
  fingerprint$rankMean <- all$rstat$mean.rank[bins.to.keep]
  fingerprint$rankSD <- apply(all$res.ranks[bins.to.keep,-c(1:3)],1,sd)
  fingerprint$A.mean <- matA$rstat$mean[bins.to.keep]
  fingerprint$A.sd <- matA$rstat$sd[bins.to.keep]
  fingerprint$A.rankMean <- matA$rstat$mean.rank[bins.to.keep]
  fingerprint$A.rankSD <- apply(matA$res.ranks[bins.to.keep,-c(1:3)],1,sd)
  fingerprint$B.mean <- matB$rstat$mean[bins.to.keep]
  fingerprint$B.sd <- matB$rstat$sd[bins.to.keep]
  fingerprint$B.rankMean <- matB$rstat$mean.rank[bins.to.keep]
  fingerprint$B.rankSD <- apply(matB$res.ranks[bins.to.keep,-c(1:3)],1,sd)

  #Generates fingerprint QC plot (if optioned)
  if(QCplot==T){

    #Sets default filename
    filename=paste("WGD.fingerprintQC.full.",Sys.Date(),".jpg",sep="")

    #Prints status
    if(quiet==F){
      cat(paste("WGDR::STATUS [",
                strsplit(as.character(Sys.time()),split=" ")[[1]][2],
                "]: plotting fingerprinting QC to ",
                OUTDIR,"\n",sep=""))
    }

    #Prepares plots
    jpeg(paste(OUTDIR,"/",filename,sep=""),height=500,width=1000,quality=150)
    par(mar=c(0.3,3.6,3.1,0.6),bty="o",mfrow=c(2,1))

    #First plot - ranks
    plot(x=range(matA$mat$start),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="",
         main=paste("WGD Bin Fingerprint Selection\n",
                    "(",prettyNum(nrow(fingerprint),big.mark=","),"/",
                    prettyNum(nrow(matA$mat),big.mark=","),
                    " bins selected; alpha=",alpha,", Min. Delta=±",minDiff,")",sep=""))
    abline(h=seq(0,1,0.1),col="gray90")
    abline(h=0.5)
    segments(x0=fingerprint$start,x1=fingerprint$end,
             y0=fingerprint$A.rankMean,y1=fingerprint$A.rankMean,
             col="blue",lwd=2)
    rect(xleft=fingerprint$start,xright=fingerprint$end,
         ybottom=fingerprint$A.rankMean-fingerprint$A.rankSD,
         ytop=fingerprint$A.rankMean+fingerprint$A.rankSD,
         col=adjustcolor("blue",alpha=0.25),border=NA)
    segments(x0=fingerprint$start,x1=fingerprint$end,
             y0=fingerprint$B.rankMean,y1=fingerprint$B.rankMean,
             col="red",lwd=2)
    rect(xleft=fingerprint$start,xright=fingerprint$end,
         ybottom=fingerprint$B.rankMean-fingerprint$B.rankSD,
         ytop=fingerprint$B.rankMean+fingerprint$B.rankSD,
         col=adjustcolor("red",alpha=0.25),border=NA)
    axis(2,at=seq(0,1,0.1),labels=seq(0,100,10),las=2,cex.axis=0.65)
    mtext(text="Ranked Coverage Percentiles",2,line=1.8)

    #Second plot - normalized coverage residuals
    par(mar=c(3.1,3.6,0.3,0.6))
    ylim <- ceiling(10*quantile(abs(c(fingerprint$A.mean,fingerprint$B.mean)),0.99))/10
    plot(x=range(matA$mat$start),y=c(-ylim,ylim),
         type="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="",main="")
    abline(h=seq(-ylim,ylim,0.1),col="gray90")
    abline(h=0)
    segments(x0=fingerprint$start,x1=fingerprint$end,
             y0=fingerprint$A.mean,y1=fingerprint$A.mean,
             col="blue",lwd=2)
    rect(xleft=fingerprint$start,xright=fingerprint$end,
         ybottom=fingerprint$A.mean-fingerprint$A.sd,
         ytop=fingerprint$A.mean+fingerprint$A.sd,
         col=adjustcolor("blue",alpha=0.25),border=NA)
    segments(x0=fingerprint$start,x1=fingerprint$end,
             y0=fingerprint$B.mean,y1=fingerprint$B.mean,
             col="red",lwd=2)
    rect(xleft=fingerprint$start,xright=fingerprint$end,
         ybottom=fingerprint$B.mean-fingerprint$B.sd,
         ytop=fingerprint$B.mean+fingerprint$B.sd,
         col=adjustcolor("red",alpha=0.25),border=NA)
    axis(2,at=seq(-ylim,ylim,0.1),las=2,cex.axis=0.65,
         labels=paste(round(seq(-ylim,ylim,0.1)*100,0),"%",sep=""))
    mtext(text="Normalized Coverage Deviance",2,line=2.3)
    axis(1,at=seq(min(matA$mat$start),max(matA$mat$end),
                  by=diff(range(matA$mat$start))/12),cex.axis=0.75,
         labels=paste(round(seq(min(matA$mat$start),max(matA$mat$end),
                                by=diff(range(matA$mat$start))/12)/1000000,0),"Mb",sep=""))
    mtext(text="Genomic Coordinate",1,line=1.8)
    dev.off()
  }

  #Further filter fingerprint bins if at least three-fold more bins selected than will be used as features
  if(4*maxBins<length(bins.to.keep)){
    #remove lowest 25% of p-values and dABs, and remove 25% of bins with largest sd
    maxsds <- apply(cbind(fingerprint$A.sd,fingerprint$B.sd),1,max)
    fingerprint.subset <- fingerprint[which(fingerprint$p<=quantile(fingerprint$p,0.75) &
                                              abs(fingerprint$dAB)>=quantile(abs(fingerprint$dAB),0.25) &
                                              maxsds<=quantile(maxsds,0.75)),]
  }

  #Subset final training bins evenly distributed across the chromosome
  if(nrow(fingerprint.subset)>maxBins){
    step <- round(nrow(bins.in.both)/maxBins,0)
    midpoints <- bins.in.both$start[seq(1,nrow(bins.in.both),step)]
    spaced.bins <- unlist(sapply(midpoints,function(val){
      which(abs(fingerprint.subset$start-val)==min(abs(fingerprint.subset$start-val)))
    }))
    repeat{
      spaced.bins <- sort(c(unique(spaced.bins),spaced.bins[duplicated(spaced.bins)]+1))
      if(length(unique(spaced.bins))>=maxBins){
        break
      }
    }
    fingerprint.subset <- fingerprint.subset[spaced.bins,]

    #Prints status
    if(quiet==F){
      cat(paste("WGDR::STATUS [",
                strsplit(as.character(Sys.time()),split=" ")[[1]][2],
                "]: filtered final fingerprint to ",nrow(fingerprint.subset),
                " most informative and evenly spaced bins\n",sep=""))
    }

    #Generates final fingerprint QC plot (if optioned)
    if(QCplot==T){

      #Sets default filename
      filename=paste("WGD.fingerprintQC.final.",Sys.Date(),".jpg",sep="")

      #Prepares plots
      jpeg(paste(OUTDIR,"/",filename,sep=""),height=500,width=1000,quality=150)
      par(mar=c(0.3,3.6,3.1,0.6),bty="o",mfrow=c(2,1))

      #First plot - ranks
      plot(x=range(matA$mat$start),y=c(0,1),type="n",
           xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="",
           main=paste("WGD Bin fingerprint.subset Selection (Final)\n",
                      "(",prettyNum(nrow(fingerprint.subset),big.mark=","),"/",
                      prettyNum(nrow(matA$mat),big.mark=","),
                      " bins applied)",sep=""))
      abline(h=seq(0,1,0.1),col="gray90")
      abline(h=0.5)
      # segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
      #          y0=fingerprint.subset$A.rankMean,y1=fingerprint.subset$B.rankMean,
      #          col="gray50")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
               y0=fingerprint.subset$A.rankMean-fingerprint.subset$A.rankSD,
               y1=fingerprint.subset$A.rankMean+fingerprint.subset$A.rankSD,
               col="blue")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
               y0=fingerprint.subset$B.rankMean-fingerprint.subset$B.rankSD,
               y1=fingerprint.subset$B.rankMean+fingerprint.subset$B.rankSD,
               col="red")
      points(x=fingerprint.subset$start,y=fingerprint.subset$A.rankMean,
             pch=23,col="blue",bg="white",cex=0.7)
      points(x=fingerprint.subset$start,y=fingerprint.subset$B.rankMean,
             pch=23,col="red",bg="white",cex=0.7)
      axis(2,at=seq(0,1,0.1),labels=seq(0,100,10),las=2,cex.axis=0.65)
      mtext(text="Ranked Coverage Percentiles",2,line=1.8)

      #Second plot - normalized coverage residuals
      par(mar=c(3.1,3.6,0.3,0.6))
      ylim <- ceiling(10*quantile(abs(c(fingerprint.subset$A.mean,fingerprint.subset$B.mean)),0.99))/10
      plot(x=range(matA$mat$start),y=c(-ylim,ylim),
           type="n",xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="",main="")
      abline(h=seq(-ylim,ylim,0.1),col="gray90")
      abline(h=0)
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
               y0=fingerprint.subset$A.mean,y1=fingerprint.subset$B.mean,
               col="gray50")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
               y0=fingerprint.subset$A.mean-fingerprint.subset$A.sd,
               y1=fingerprint.subset$A.mean+fingerprint.subset$A.sd,
               col="blue")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$end,
               y0=fingerprint.subset$A.mean-fingerprint.subset$A.sd,
               y1=fingerprint.subset$A.mean-fingerprint.subset$A.sd,
               col="blue")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$end,
               y0=fingerprint.subset$A.mean+fingerprint.subset$A.sd,
               y1=fingerprint.subset$A.mean+fingerprint.subset$A.sd,
               col="blue")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$start,
               y0=fingerprint.subset$B.mean-fingerprint.subset$B.sd,
               y1=fingerprint.subset$B.mean+fingerprint.subset$B.sd,
               col="red")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$end,
               y0=fingerprint.subset$B.mean-fingerprint.subset$B.sd,
               y1=fingerprint.subset$B.mean-fingerprint.subset$B.sd,
               col="red")
      segments(x0=fingerprint.subset$start,x1=fingerprint.subset$end,
               y0=fingerprint.subset$B.mean+fingerprint.subset$B.sd,
               y1=fingerprint.subset$B.mean+fingerprint.subset$B.sd,
               col="red")
      points(x=fingerprint.subset$start,y=fingerprint.subset$A.mean,
             pch=23,col="blue",bg="white",cex=0.7)
      points(x=fingerprint.subset$start,y=fingerprint.subset$B.mean,
             pch=23,col="red",bg="white",cex=0.7)
      axis(2,at=seq(-ylim,ylim,0.1),las=2,cex.axis=0.65,
           labels=paste(round(seq(-ylim,ylim,0.1)*100,0),"%",sep=""))
      mtext(text="Normalized Coverage Deviance",2,line=2.3)
      axis(1,at=seq(min(matA$mat$start),max(matA$mat$end),
                    by=diff(range(matA$mat$start))/12),cex.axis=0.75,
           labels=paste(round(seq(min(matA$mat$start),max(matA$mat$end),
                                  by=diff(range(matA$mat$start))/12)/1000000,0),"Mb",sep=""))
      mtext(text="Genomic Coordinate",1,line=1.8)
      dev.off()
    }
  }

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: fitting GMM to training data\n",sep=""))
  }

  #Create transposed, merged matrix of all training data for GMM
  train.matrix <- t(merge(all$res,fingerprint.subset[,1:3],sort=F)[,-c(1:3)])
  colnames(train.matrix) <- apply(fingerprint.subset[,1:3],1,function(coords){
    return(paste(coords,collapse="_"))
  })
  train.matrix <- as.data.frame(apply(train.matrix,2,function(vals){
    if(!(any(is.na(vals)))){
      vals <- scale(as.numeric(vals),center=T)
      return(vals)
    }
  }))
  rownames(train.matrix) <- colnames(all$mat[,-c(1:3)])

  #Fits model to training data
  assignments <- as.factor(c(rep("A",ncol(matA$mat)-3),
                             rep("B",ncol(matB$mat)-3)))
  GMM <- MclustDA(data=train.matrix,class=assignments)
#
#   #Prints NxN clustering matrix if QC plots are optioned
#   if(QCplot==T){
#
#     #Sets default filename
#     filename=paste("WGD.GMM_QC.",Sys.Date(),".jpg",sep="")
#
#     #Plots clusters
#     jpeg(paste(OUTDIR,"/",filename,sep=""),
#          height=50*ncol(train.matrix),
#          width=50*ncol(train.matrix))
#     clPairs(train.matrix,assignments)
#     dev.off()
#   }

  #Returns fingerprint & model
  return(list("fp.full"=fingerprint,
              "fp.final"=fingerprint.subset,
              "model"=GMM))
}

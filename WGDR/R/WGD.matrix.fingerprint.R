# Copyright (c) 2016 Ryan L. Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# WGD R Companion Library (WGDR) Function

################
# WGD.matrix.fingerprint
################
# Isolates top N% of most variable bins from WGDR matrix
################
# Returns a fingerprint list object:
#  $fingerprint : BED-style data frame of all bins used in fingerprinting with median per assignment group
#  $assignments : vector of assignments for all samples in input matrix
################

WGD.matrix.fingerprint <- function(mat,            #matrix object from which to read. Must be read with WGD.readmatrix
                                   fraction=0.10,  #fraction of bins to use for finterprinting
                                   k=2,            #number of clusters to use during PCA
                                   heatmap=F,      #option to plot a heatmap of resulting fingerprints
                                   OUTDIR=NA,      #output directory for heatmap plot
                                   filename=NULL,  #add custom filename for output
                                   quiet=F         #option to disable verbose output
){
  #Sanity check fraction parameter
  if(is.numeric(fraction)==F | fraction<=0 | fraction>1){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'fraction' must be float ~ (0,1]",sep=""))
  }

  #Sanity check k parameter
  if(is.numeric(k)==F | k<1 | k>2){
    stop(paste("WGDR::ERROR [",
               strsplit(as.character(Sys.time()),split=" ")[[1]][2],
               "] ",match.call()[[1]]," parameter 'k' must be integer ~ [1,2]",sep=""))
  }

  #Disable scientific notation
  options(scipen=1000)

  #Require pcaMethods and RColorBrewewr libraries
  require(pcaMethods)
  require(RColorBrewer)

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: fingerprinting coverage matrix (fraction=",fraction,")\n",sep=""))
  }

  #Isolates fingerprint bins
  fingerprint <- mat$stat[which(mat$lstat$range95pct>=quantile(mat$lstat$range95pct,1-fraction)),1:3]

  #Prints status
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: isolated ",prettyNum(nrow(fingerprint),big.mark=","),
              " most variable bins for fingerprinting\n",sep=""))
  }
  if(quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: assigning samples versus log2-fold change fingerprints\n",sep=""))
  }

  #Fingerprints samples
  lstat.finger <- mat$lstat[which(mat$lstat$range95pct>=quantile(mat$lstat$range95pct,1-fraction)),]
  log2.finger <- mat$log2[which(mat$lstat$range95pct>=quantile(mat$lstat$range95pct,1-fraction)),]
  log2.finger[,-c(1:3)] <- apply(log2.finger[,-c(1:3)],2,function(vals){
    vals[which(!(is.finite(vals)))] <- NA
    return(vals)
  })

  #Runs imputation for infinite & NA values, then performs PCA & kmeans
  pc <- pca(log2.finger[,-c(1:3)],nPcs=4,method="ppca")
  imputed <- completeObs(pc)
  PCA <- prcomp(imputed,center=T,scale=T)
  if(k>1){
    assignments <- kmeans(PCA$rotation[,1:2],2)$cluster
  }else{
    assignments <- rep(1,times=ncol(imputed))
  }

  #Calculate average log2-fold change per group at each fingerprinting
  group.medians <- t(sapply(which(mat$lstat$range95pct>=quantile(mat$lstat$range95pct,1-fraction)),
                            function(row){
                              if(k>1){
                                g1 <- as.numeric(mat$log2[row,which(assignments==1)+3])
                                g1[which(!(is.finite(g1)))] <- NA
                                g1 <- median(g1,na.rm=T)
                                g2 <- as.numeric(mat$log2[row,which(assignments==2)+3])
                                g2[which(!(is.finite(g2)))] <- NA
                                g2 <- median(g2,na.rm=T)
                                return(c(g1,g2))
                              }else{
                                g <- as.numeric(mat$log2[row,-c(1:3)])
                                g[which(!(is.finite(g)))] <- NA
                                return(median(g,na.rm=T))
                              }
                            }))
  fingerprint <- cbind(fingerprint,group.medians)
  names(fingerprint)[4:5] <- c("g1","g2")

  #Orders samples by cluster & average absolute value fingerprint
  assignments.mod <- (2*assignments)-3
  abssums <- apply(log2.finger[,-c(1:3)],2,function(vals){
    vals[which(!(is.finite(vals)))] <- NA
    return(median(abs(vals),na.rm=T))
  })
  abssums.mod <- abssums*assignments.mod
  assignments.order <- order(abssums.mod)
  finger.order <- rev(order(fingerprint[,4]-fingerprint[,5]))

  #Plots heatmap if optioned
  if(heatmap==T){
    if(is.na(OUTDIR)){
      stop(paste("WGDR::ERROR [",
                 strsplit(as.character(Sys.time()),split=" ")[[1]][2],
                 "] ",match.call()[[1]]," parameter 'OUTDIR' must be set for heatmap plotting",sep=""))
    }
    if(quiet==F){
      cat(paste("WGDR::STATUS [",
                strsplit(as.character(Sys.time()),split=" ")[[1]][2],
                "]: plotting heatmap of sample fingerprints by assignment\n",sep=""))
    }
    if(!(dir.exists(OUTDIR))){
      dir.create(OUTDIR)
    }
    #Set diverging color palette
    pal <- colorRampPalette(c("red","white","blue"))
    #Restrict log2 fingerprints on -2x to 2x (-1 to 1)
    log2.finger.plot <- t(log2.finger[finger.order,assignments.order+3])
    log2.finger.plot[which(!(is.finite(log2.finger.plot)))] <- NA
    log2.finger.plot[which(log2.finger.plot < -1)] <- -1
    log2.finger.plot[which(log2.finger.plot>1)] <- 1
    #Prepare layout
    if(!is.null(filename)){
      outname <- paste(OUTDIR,filename,sep="")
    }else{
      outname <- paste(OUTDIR,"WGD_assignment_heatmap.jpg",sep="")
    }
    jpeg(outname,quality=100,res=300,
         height=300*8.5*0.8,width=300*8.5*0.8,
         title=paste("WGD_PCA",sep=""))
    par(mar=c(2.1,2.1,0.6,0.6),bty="n")
    plot(x=c(1,ncol(log2.finger.plot)),
         y=c(1,nrow(log2.finger.plot)),
         xaxs="i",yaxs="i",ylab="",xlab="",xaxt="n",yaxt="n")
    for(i in 1:nrow(log2.finger.plot)){
      colvals <- ceiling(128*(log2.finger.plot[i,]+1))
      colvals[which(colvals==0)] <- 1
      rect(xleft=0:(ncol(log2.finger.plot)-1),
           xright=1:ncol(log2.finger.plot),
           ybottom=i-1,ytop=i,border=NA,
           col=pal(256)[colvals])
    }
    mtext(1,line=0,text="Genomic Bins")
    mtext(2,line=0,text="Samples")
    dev.off()
  }

  #Returns fingerprint bins and assignments
  return(list("fingerprint"=fingerprint,
              "assignments"=assignments))
}

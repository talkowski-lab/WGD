#!/usr/bin/env Rscript

# Copyright (c) 2017 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Script to analyze a matrix of median binCov values per contig per sample
# for approximate copy numbers (with visualization)

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)

###################################################
#####Helper function to load median coverage matrix
###################################################
readMatrix <- function(INFILE){
  dat <- read.table(INFILE,comment.char="",header=T)
  colnames(dat)[1] <- "ID"
  dat[,-1] <- t(apply(dat[,-1],1,as.numeric))
  return(dat)
}

#############################################################
#####Helper function to normalize contigs for a single sample
#############################################################
normalizeContigsPerSample <- function(vals,exclude=NA,ploidy){
  #Convert vals to numeric
  vals <- as.numeric(vals)

  #Iterate over values
  newVals <- sapply(1:length(vals),function(i){
    #Compute mean of values excluding current value & any other specified
    excl.mean <- mean(vals[-c(i,exclude)],na.rm=T)

    #Normalize current value versus appropriate mean
    return(vals[i]/excl.mean)
  })

  #Scale to expected ploidy
  newVals <- ploidy*newVals

  #Return normalized values
  return(newVals)
}

#########################################################################
#####Helper function to normalize contigs for an entire matrix of samples
#########################################################################
normalizeContigsPerMatrix <- function(dat,exclude=NA,scale.exclude=NA,
                                      genome.ploidy=2,contig.ploidy){
  #Iterate over samples & normalize
  suppressWarnings(if(is.na(exclude)){
    dat[,-1] <- t(apply(dat[,-1],1,normalizeContigsPerSample,genome.ploidy))
  }else{
    dat[,-1] <- t(apply(dat[,-1],1,normalizeContigsPerSample,exclude=exclude-1,genome.ploidy))
  })

  #Scale mad to mad of first 12 chromosomes
  mad.others <- mad(unlist(dat[,2:13]),na.rm=T)

  #Iterate over contigs (minus scale.exclude) and scale
  scaledVals <- sapply(setdiff(2:ncol(dat),scale.exclude),function(i){
    #Calculate & apply adjustments
    median.adjust <- median(dat[,i],na.rm=T)-contig.ploidy[i-1]
    newvals <- dat[,i]-median.adjust
    return(newvals)
  })
  suppressWarnings(if(is.na(scale.exclude)){
    dat[,-1] <- scaledVals
  }else{
    dat[,-c(1,scale.exclude)] <- scaledVals
  })

  #Round up values that were normalized below zero
  dat[,-1] <- apply(dat[,-1],2,function(vals){
    vals[which(vals<0 & !is.na(vals))] <- 0
    return(vals)
  })

  #Return transformed data
  return(dat)
}

#######################################################################
#####Helper function to test each contig per sample for evidence of CNA
#######################################################################
testCNs <- function(dat,FDR=T){
  dat.out <- dat

  #Iterate over contigs
  left.p <- apply((apply(dat[,-1],2,scale,center=T,scale=T)),2,pnorm)
  right.p <- apply((apply(dat[,-1],2,scale,center=T,scale=T)),2,pnorm,lower.tail=F)

  #Choose minimum p-value between left and right tails
  dat.out[,-1] <- t(sapply(1:nrow(left.p),function(row){
    pvals <- sapply(1:ncol(left.p),function(col){
      return(min(left.p[row,col],right.p[row,col],na.rm=T))
    })
    return(pvals)
  }))

  #FDR correct all p-values (if optioned)
  if(FDR==T){
    dat.out[,-1] <- apply(dat.out[,-1],2,p.adjust,method="fdr")
  }

  #Return
  return(dat.out)
}

###########################################
#####Helper function to assign & plot sexes
###########################################
assignSex <- function(dat,sexChr=24:25,
                      sexAssign.df, #four-column df: cn(X), cn(Y), label, color
                      mosaicThresh="Bonferroni",
                      plot=T,axLim=3){
  #Exclude incomplete entries
  sample.exclude <- unlist(sapply(1:nrow(dat),function(i){
    if(any(is.na(dat[i,sexChr]))){
      return(i)
    }
  }))
  dat.mod <- dat[-sample.exclude,c(1,sexChr)]
  if(nrow(dat)!=nrow(dat.mod)){
    warning(paste(nrow(dat)-nrow(dat.mod),
                  " samples missing sex chromosome coverage information and were excluded from sex assignments.",
                  sep=""))
  }

  #Round X & Y CN-types to nearest whole integer
  dat.mod.rounded <- dat.mod
  dat.mod.rounded[,-1] <- apply(dat.mod[,-1],2,round,digits=0)

  #Create output data frame with sex assignments
  sexes <- as.data.frame(t(unlist(sapply(unique(dat[,1]),function(ID){
    #Iterate over all IDs
    if(ID %in% dat.mod[,1]){
      #Get rounded CNs
      CN.X <- round(dat.mod.rounded[which(dat.mod.rounded$ID==ID),2],0)
      CN.Y <- round(dat.mod.rounded[which(dat.mod.rounded$ID==ID),3],0)

      #Assign to existing sex copy-type, or "OTHER"
      if(length(which(sexAssign.df$CN.X==CN.X & sexAssign.df$CN.Y==CN.Y))==1){
        return(as.vector(c(ID,CN.X,CN.Y,
                           sexAssign.df[which(sexAssign.df$CN.X==CN.X & sexAssign.df$CN.Y==CN.Y),3])))
      }else{
        return(as.vector(c(ID,CN.X,CN.Y,"OTHER")))
      }
    }else{
      return(as.vector(c(ID,
               dat[which(dat$ID==ID),sexChr[1]],
               dat[which(dat$ID==ID),sexChr[2]],
               NA)))
    }
  }))))
  colnames(sexes) <- c("ID","chrX.CN","chrY.CN","Assignment")
  rownames(sexes) <- 1:nrow(sexes)

  #Gather sd of X and Y assignments from males
  sd.X <- sd(dat.mod[which(sexes$Assignment[-sample.exclude]=="MALE"),2],na.rm=T)
  sd.Y <- sd(dat.mod[which(sexes$Assignment[-sample.exclude]=="MALE"),3],na.rm=T)

  #Run mosaic check per sample
  pMosaic <- as.data.frame(t(unlist(apply(sexes,1,function(vals){
    if(any(is.na(vals))){
      pMos.X <- NA
      pMos.Y <- NA
    }else{
      #Get raw CNs
      CN.X <- dat.mod[which(dat.mod.rounded$ID==vals[1]),2]
      CN.Y <- dat.mod[which(dat.mod.rounded$ID==vals[1]),3]

      #Calculate p-value for chrX mosaicism
      if(CN.X>as.numeric(vals[2])){
        pMos.X <- pnorm(CN.X,mean=as.numeric(vals[2]),
                        sd=sd.X*as.numeric(vals[2]),lower.tail=F)
      }else{
        pMos.X <- pnorm(CN.X,mean=as.numeric(vals[2]),
                        sd=sd.X*as.numeric(vals[2]),lower.tail=T)
      }

      #Calculate p-value for chrY mosaicism
      if(CN.Y>as.numeric(vals[3])){
        pMos.Y <- pnorm(CN.Y,mean=as.numeric(vals[3]),
                        sd=sd.Y*max(as.numeric(vals[3]),1),lower.tail=F)
      }else{
        pMos.Y <- pnorm(CN.Y,mean=as.numeric(vals[3]),
                        sd=sd.Y*max(as.numeric(vals[3]),1),lower.tail=T)
      }
    }

    #Return p-values
    return(c(pMos.X,pMos.Y))
  }))))
  colnames(pMosaic) <- c("pMos.X","pMos.Y")
  pMosaic$qMos.X <- p.adjust(pMosaic$pMos.X,method="fdr")
  pMosaic$qMos.Y <- p.adjust(pMosaic$pMos.Y,method="fdr")

  #Add mosaic check to sexes output
  sexes <- cbind(sexes,pMosaic)

  #####Plot sex assignments
  if(plot==T){
    #Prepare plot area
    par(mar=c(3.5,3.5,0.5,0.5))
    plot(x=c(0,axLim),y=c(0,axLim),type="n",
         xlab="",ylab="",xaxt="n",yaxt="n")

    #Add grid lines
    abline(h=0:axLim,v=0:axLim,lty=3,col="gray50")

    #Add sex karyotype labels behind each centroid
    sexCNs.df <- data.frame("X"=rep(0:axLim,axLim+1),
                            "Y"=as.vector(unlist(sapply(0:axLim,rep,times=axLim+1))))
    sexCNs.df$karyo <- apply(sexCNs.df,1,function(vals){
      return(paste(paste(rep("X",vals[1]),collapse=""),
                   paste(rep("Y",vals[2]),collapse=""),
                   sep=""))
    })
    text(x=sexCNs.df[,1],y=sexCNs.df[,2],
         labels=sexCNs.df$karyo,
         font=2,col="gray95",cex=0.8)

    #Assign colors for sex plotting
    if(mosaicThresh=="FDR"){
      colVect <- apply(sexes[-sample.exclude,c(4,7:8)],1,function(vals){
        if(min(as.numeric(vals[2:3]))<0.05){
          return("#53edd0")
        }else if(vals[1] %in% sexAssign.df$label){
          return(sexAssign.df[which(sexAssign.df$label==vals[1]),]$color)
        }else{
          return("#8F1336")
        }
      })
    }else{
      colVect <- apply(sexes[-sample.exclude,c(4:6)],1,function(vals){
        if(min(as.numeric(vals[2:3]))<0.05/nrow(dat.mod)){
          return("#53edd0")
        }else if(vals[1] %in% sexAssign.df$label){
          return(sexAssign.df[which(sexAssign.df$label==vals[1]),]$color)
        }else{
          return("#8F1336")
        }
      })
    }

    #Plot points
    points(dat.mod[,-1],pch=19,col=colVect,cex=0.5)

    #Add x-axis
    axis(1,at=0:axLim)
    mtext(1,line=2.2,text="chrX Copy Number")

    #Add y-axis
    axis(2,at=0:axLim,las=2)
    mtext(2,line=2.2,text="chrY Copy Number")

    #Add legend
    legendLabs <- apply(sexAssign.df[,1:3],1,function(vals){
      paste(vals[3]," (",
            paste(rep("X",times=vals[1]),collapse=""),
            paste(rep("Y",times=vals[2]),collapse=""),
            ")",sep="")
    })
    legend("topright",bg="white",
           legend=c(legendLabs,"MOSAIC","OTHER"),
           pch=19,col=c(sexAssign.df$color,"#53edd0","#8F1336"),
           cex=0.75,pt.cex=1)
  }

  #Return sex assignments
  return(sexes)
}

###############################################################
#####Helper function to plot distribution of samples per contig
###############################################################
boxplotsPerContig <- function(dat,exclude,genome.ploidy=2,contig.ploidy,
                              contigLabels=paste("chr",c(1:22,"X","Y"),sep=""),
                              colorSignif=T,connect=F){
  #Load library
  require(beeswarm)

  #Create plot color dataframe
  if(colorSignif==T){
    #Calculate p-values
    pvals <- testCNs(dat,FDR=T)

    #Iterate over pvalue matrix & compute plot color
    colMat <- sapply(2:ncol(pvals),function(col){
      sapply(1:nrow(pvals),function(row){
        if(is.na(pvals[row,col])){
          return(NA)
        }else{
          if(pvals[row,col]<0.05){
            if(dat[row,col]>contig.ploidy[col-1]){
              return("blue")
            }else{
              return("red")
            }
          }else{
            return("#838393")
          }
        }
      })
    })
  }else{
    colMat <- dat[,-1]
    colMat <- "#838393"
  }

  #Get max y-value
  ymax <- max(4,max(dat[,-1],na.rm=T))

  #Prepare plot area
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(x=c(0,ncol(dat)-1),y=c(0,ymax),
       type="n",yaxt="n",xaxt="n",xaxs="i",ylab="",xlab="")

  #Add shading
  rect(xleft=seq(0,ncol(dat)+1,2),xright=seq(1,ncol(dat)+2,2),
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border=NA,col="gray95")

  #Add gridlines
  abline(h=seq(0,ymax,0.5),lty=2,col="gray85")
  abline(h=0:ymax,col="gray80")
  abline(h=c(0,genome.ploidy))

  #Iterate over contigs and plot all data
  if(connect==T){
    sapply(1:nrow(dat),function(i){
      vals <- dat[i,-1]
      cols <- colMat[i,]
      sapply(1:length(vals),function(j){
        if(cols[j]!="#838393" & !is.na(cols[j])){
          if(j==1){
            segments(x0=j-0.5,x1=j+0.5,
                     y0=as.numeric(vals[j]),
                     y1=as.numeric(vals[j+1]),
                     col="gray80",lwd=0.4)
          }else if(j==length(vals)){
            segments(x0=j-1.5,x1=j-0.5,
                     y0=as.numeric(vals[j-1]),
                     y1=as.numeric(vals[j]),
                     col="gray80",lwd=0.4)
          }else{
            segments(x0=j-1.5,x1=j-0.5,
                     y0=as.numeric(vals[j-1]),
                     y1=as.numeric(vals[j]),
                     col="gray80",lwd=0.4)
            segments(x0=j-0.5,x1=j+0.5,
                     y0=as.numeric(vals[j]),
                     y1=as.numeric(vals[j+1]),
                     col="gray80",lwd=0.4)
          }
        }
      })
    })
  }
  sapply(1:(ncol(dat)-1),function(i){
    if(connect==F){
      #Jitter
      points(x=jitter(rep(i-0.5,times=nrow(dat)),amount=0.3),y=dat[,i+1],
             pch=19,col=colMat[,i],cex=0.25)
      # #Swarm
      # beeswarm(dat[,i+1],add=T,at=i-0.5,method="swarm",
      #          corral="wrap",corralWidth=0.6,
      #          pch=19,pwcol=colMat[,i],cex=0.2)
    }else{
      points(x=rep(i-0.5,times=nrow(dat)),y=dat[,i+1],
             pch=19,col=colMat[,i],cex=0.25)
    }
  })

  #Add boxplots
  if(!is.na(exclude)){
    boxplot(dat[,-c(1,exclude)],at=(1:(ncol(dat)-1))[-c(exclude-1)]-0.5,
            add=T,outline=F,col=NA,lwd=0.75,lty=1,staplewex=0,
            yaxt="n",xaxt="n",ylab="",xlab="")
  }else{
    boxplot(dat[,-1],at=(1:(ncol(dat)-1))-0.5,
            add=T,outline=F,col=NA,lwd=0.75,lty=1,staplewex=0,
            yaxt="n",xaxt="n",ylab="",xlab="")
  }

  #Add x-axis labels
  axis(1,at=(1:length(contigLabels))-0.5,tick=F,line=-0.8,las=2,labels=contigLabels)
  mtext(1,text="Chromosome",line=2.2)

  #Add y-axis labels
  axis(2,at=0:ymax,las=2)
  mtext(2,text="Estimated Copy Number",line=2.2)

  #Add text in top-left indicating number of samples
  nSamp <- nrow(dat)
  nSamp.hasNA <- length(unique(unlist(sapply(1:nrow(dat[,-1]),function(i){
    if(any(is.na(dat[i,-1]))){
      return(i)
    }else{
      return(NA)
    }
  }))))
  text(x=par("usr")[1],y=0.975*ymax,
       labels=paste("N=",prettyNum(nSamp,big.mark=",")," samples (",
                    prettyNum(nSamp.hasNA,big.mark=",")," incomplete)",
                    sep=""),pos=4)

  #Add border cleanup
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA,border="black",lwd=2)
}

##########################
#####Rscript functionality
##########################
require(optparse)
#List of command-line options
option_list <- list(
  make_option(c("-O", "--OUTDIR"),type="character",default=NULL,
              help="output directory [default: pwd]",
              metavar="character"),
  make_option(c("-p", "--noplot"),action="store_true",default=FALSE,
              help="disable copy number visualization [default: FALSE]"),
  make_option(c("-z", "--gzip"),action="store_false",default=TRUE,
              help="gzip output files [default: TRUE]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] median_coverage_matrix",
                                option_list=option_list),
                   positional_arguments=TRUE)
INFILE <- args$args[1]
OUTDIR <- args$options$OUTDIR
noplot <- args$options$noplot
gzip <- args$options$gzip
if(is.null(OUTDIR)){
  OUTDIR <- ""
}

#Create OUTDIR if it doesn't already exist
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

#Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Must supply an input median coverage matrix\n")
}

#Loads data
dat <- readMatrix(INFILE)

#Transforms data to predicted copy numbers
dat.norm <- normalizeContigsPerMatrix(dat,exclude=24:25,scale.exclude=24:25,
                                      genome.ploidy=2,contig.ploidy=rep(2,24))

#Set sex assignment table
sexAssign.df <- data.frame("CN.X"=c(1,2,1,3,2,1),
                           "CN.Y"=c(1,0,0,0,1,2),
                           "label"=c("MALE","FEMALE","TURNER",
                                     "TRIPLE X","KLINEFELTER","JACOBS"),
                           "color"=c("#00BFF4","#fd8eff","#e02006",
                                     "#7B2AB3","#FF6A09","#29840f"))

#Assign sexes
sexes <- assignSex(dat.norm,plot=F,sexAssign.df=sexAssign.df)

#Split data by chrX copy number (≥2 or <2) & evenly distribute NAs among M & F
sex.males <- as.vector(which(round(as.numeric(sexes$chrX.CN,0))<2 & !is.na(sexes$chrX.CN)))
sex.females <- as.vector(which(round(as.numeric(sexes$chrX.CN,0))>=2 & !is.na(sexes$chrX.CN)))
sex.NAs <- as.vector(which(is.na(sexes$chrX.CN)))
sex.NAs.first <- sex.NAs[1:floor((length(sex.NAs)/2))]
sex.NAs.second <- sex.NAs[(floor((length(sex.NAs)/2))+1):length(sex.NAs)]
dat.males <- dat[c(sex.males,sex.NAs.first),]
dat.females <- dat[c(sex.females,sex.NAs.second),]

#Normalize CN - males & females separately
males.norm <- normalizeContigsPerMatrix(dat.males,exclude=24:25,scale.exclude=NA,
                                        genome.ploidy=2,contig.ploidy=c(rep(2,22),1,1))
females.norm <- normalizeContigsPerMatrix(dat.females,exclude=24:25,scale.exclude=NA,
                                          genome.ploidy=2,contig.ploidy=c(rep(2,22),2,0))

#Plot CN per contig - Males
png(paste(OUTDIR,"/estimated_CN_per_contig.chrX_lessThan_2copies.with_contours.png",sep=""),
    height=1250,width=2500,res=300)
boxplotsPerContig(males.norm,exclude=NA,contig.ploidy=c(rep(2,22),1,1),connect=T)
dev.off()
png(paste(OUTDIR,"/estimated_CN_per_contig.chrX_lessThan_2copies.no_contours.png",sep=""),
    height=1250,width=2500,res=300)
boxplotsPerContig(males.norm,exclude=NA,contig.ploidy=c(rep(2,22),1,1),connect=F)
dev.off()

#Plot CN per contig - Females
png(paste(OUTDIR,"/estimated_CN_per_contig.chrX_atLeast_2copies.with_contours.png",sep=""),
    height=1250,width=2500,res=300)
boxplotsPerContig(females.norm,exclude=NA,contig.ploidy=c(rep(2,22),2,0),connect=T)
dev.off()
png(paste(OUTDIR,"/estimated_CN_per_contig.chrX_atLeast_2copies.no_contours.png",sep=""),
    height=1250,width=2500,res=300)
boxplotsPerContig(females.norm,exclude=NA,contig.ploidy=c(rep(2,22),2,0),connect=F)
dev.off()

#Recombine independently normalized data by sex & reassign sexes
merged.norm <- rbind(males.norm,females.norm)

#Plots sex assignment dotplot
png(paste(OUTDIR,"/sex_assignments.png",sep=""),
    height=1500,width=1500,res=300)
sexes <- assignSex(merged.norm,sexAssign.df=sexAssign.df)
dev.off()

#Write table of sexes
sexes <- sexes[match(dat$ID,sexes$ID),]
colnames(sexes)[1] <- "#ID"
write.table(sexes,paste(OUTDIR,"/sample_sex_assignments.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
if(gzip==T){
  system(paste("gzip -f ",OUTDIR,"/sample_sex_assignments.txt",sep=""),intern=F,wait=F)
}

#Generate p-values, q-values, and rounded CNs for males/females
males.p <- testCNs(males.norm,FDR=F)
males.q <- testCNs(males.norm,FDR=T)
males.CN <- males.norm
males.CN[,-1] <- apply(males.CN[,-1],2,round,digits=2)
females.p <- testCNs(females.norm,FDR=F)
females.q <- testCNs(females.norm,FDR=T)
females.CN <- females.norm
females.CN[,-1] <- apply(females.CN[,-1],2,round,digits=2)

#Merge male/female p-values and rounded CNs
merged.p <- rbind(males.p,females.p)
merged.p <- merged.p[match(dat$ID,merged.p$ID),]
colnames(merged.p) <- c("#ID",paste("chr",c(1:22,"X","Y"),"_pValue",sep=""))
merged.q <- rbind(males.q,females.q)
merged.q <- merged.q[match(dat$ID,merged.q$ID),]
colnames(merged.q) <- c("#ID",paste("chr",c(1:22,"X","Y"),"_qValue",sep=""))
merged.CN <- rbind(males.CN,females.CN)
merged.CN <- merged.CN[match(dat$ID,merged.CN$ID),]
colnames(merged.CN) <- c("#ID",paste("chr",c(1:22,"X","Y"),"_CopyNumber",sep=""))

#Write merged p-values
write.table(merged.p,paste(OUTDIR,"/CNA_pValues.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
if(gzip==T){
  system(paste("gzip -f ",OUTDIR,"/CNA_pValues.txt",sep=""),intern=F,wait=F)
}

#Write merged q-values
write.table(merged.q,paste(OUTDIR,"/CNA_qValues.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
if(gzip==T){
  system(paste("gzip -f ",OUTDIR,"/CNA_qValues.txt",sep=""),intern=F,wait=F)
}

#Write merged copy number estimates
write.table(merged.CN,paste(OUTDIR,"/estimated_copy_numbers.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
if(gzip==T){
  system(paste("gzip -f ",OUTDIR,"/estimated_copy_numbers.txt",sep=""),intern=F,wait=F)
}

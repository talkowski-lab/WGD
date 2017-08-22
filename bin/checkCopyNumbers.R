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
normalizeContigsPerSample <- function(vals,exclude){
  #Convert vals to numeric
  vals <- as.numeric(vals)

  #Iterate over values
  newVals <- sapply(1:length(vals),function(i){
    #Compute mean of values excluding current value & any other specified
    excl.mean <- mean(vals[-c(i,exclude)],na.rm=T)

    #Normalize current value versus appropriate mean
    return(vals[i]/excl.mean)
  })

  #Return normalized values
  return(newVals)
}

#########################################################################
#####Helper function to normalize contigs for an entire matrix of samples
#########################################################################
normalizeContigsPerMatrix <- function(dat,exclude,ploidy=2){
  #Iterate over samples & normalize
  dat[,-1] <- ploidy*t(apply(dat[,-1],1,normalizeContigsPerSample,exclude=exclude-1))

  #Scale sd to sd of first 12 chromosomes
  sd.others <- sd(unlist(dat[,2:13]),na.rm=T)

  #Iterate over contigs (minus excluded) and scale
  scaledVals <- sapply(setdiff(2:ncol(dat),exclude),function(i){
    #Calculate & apply adjustments
    mean.adjust <- mean(dat[,i],na.rm=T)-ploidy
    newvals <- dat[,i]-mean.adjust
    sd.adjust <- sd.others/sd(dat[,i],na.rm=T)
    newvals <- 2+((newvals-2)*sd.adjust)
    return(newvals)
  })
  dat[,-c(1,exclude)] <- scaledVals

  #Return transformed data
  return(dat)
}

###########################################
#####Helper function to assign & plot sexes
###########################################
assignSex <- function(dat,sexChr=24:25,
                      sexCNs=list(c(1,1),c(2,0)),
                      sexLabs=c("MALE","FEMALE"),
                      sexColors=c("#2777f7","#f027f7"),
                      plot=T,axLim=3){
  #Create data frame of expected sex CNs
  sexCNs.df <- as.data.frame(t(matrix(unlist(sexCNs),ncol=length(sexChr))))
  colnames(sexCNs.df) <- colnames(dat[,sexChr])

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

  #Perform k-means clustering
  set.seed(20)
  clusters <- kmeans(x=dat.mod[,-1],centers=sexCNs.df)

  #Create output data frame with sex assignments
  sexes <- sapply(1:nrow(dat),function(i){
    #Get sex assignment
    if(as.character(i) %in% names(clusters$cluster)){
      cIdx <- clusters$cluster[which(names(clusters$cluster)==as.character(i))]
      sex <- sexLabs[cIdx]
    }else{
      sex <- NA
    }

    #Return info
    return(sex)
  })
  names(sexes) <- dat[,1]

  #####Plot sex assignments
  #Prepare plot area
  par(mar=c(3.5,3.5,0.5,0.5))
  plot(x=c(0,axLim),y=c(0,axLim),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")

  #Add grid lines
  abline(h=0:axLim,v=0:axLim,lty=3,col="gray50")

  #Get color vector
  colVect <- as.vector(unlist(sapply(sexes,function(sex){
    if(is.na(sex)){
      col=NA
    }else{
      col <- sexColors[which(sexLabs==sex)]
    }
    return(col)
  })))

  #Plot points
  points(dat[,sexChr],pch=19,col=colVect,cex=0.5)

  #Add x-axis
  axis(1,at=0:axLim)
  mtext(1,line=2.2,text="chrX Copy Number")

  #Add y-axis
  axis(2,at=0:axLim,las=2)
  mtext(2,line=2.2,text="chrY Copy Number")

  #Add legend
  legend("topright",bg="white",legend=sexLabs,
         pch=c(19,19),col=sexColors,cex=2)

  #Return sex assignments
  return(sexes)
}

# ########################################
# #####Helper function to calculate CN fit
# ########################################
# fitCN <- function(dat,exclude){
#   #Exclude incomplete entries
#   sample.exclude <- unlist(sapply(1:nrow(dat),function(i){
#     if(all(is.na(dat[i,-1]))){
#       return(i)
#     }
#   }))
#   dat.mod <- dat[-sample.exclude,-c(1,exclude)]
#   if(nrow(dat)!=nrow(dat.mod)){
#     warning(paste(nrow(dat)-nrow(dat.mod),
#                   " samples had no coverage information and were excluded from CN modeling.",
#                   sep=""))
#   }
#
#   #Compute PCs
#   PCs <- prcomp(dat.mod,retx=T)
#
#   #Bind first 10 PCs to dat
#   dat.mod <- cbind(dat.mod,PCs$x[,1:10])
#
#   #Iterate over remaining contigs and fit linear model
#   sapply(1:(ncol(dat.mod)-10),function(i){
#     #Fit linear model
#     fit <- lm(as.formula(paste(colnames(dat.mod)[i],"~",
#                         paste(colnames(dat.mod)[-i],
#                               collapse="+"),sep = "")),
#        dat=dat.mod)
#     #Calculate fits on full data
#     fit.vals <- predict.lm(fit,dat.mod)
#   })
#
# }

###############################################################
#####Helper function to plot distribution of samples per contig
###############################################################
boxplotsPerContig <- function(dat,ploidy=2,exclude,
                              contigLabels=paste("chr",c(1:22,"X","Y"),sep="")){
  #Load library
  require(beeswarm)

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
  abline(h=c(0,ploidy))

  #Iterate over contigs and plot violins & jitters
  sapply(1:(ncol(dat)-1),function(i){
    #Jitter
    points(x=jitter(rep(i-0.5,times=nrow(dat)),amount=0.3),y=dat[,i+1],
           pch=19,col="#3e68ad",cex=0.25)
    # #Swarm
    # beeswarm(dat[,i+1],add=T,at=i-0.5,method="swarm",
    #          corral="wrap",corralWidth=0.6,
    #          pch=19,col="#3e68ad",cex=0.2)

  })

  #Add boxplots
  boxplot(dat[,-c(1,exclude)],at=(1:(ncol(dat)-1))[-c(exclude-1)]-0.5,
          add=T,outline=F,col=NA,lwd=0.75,lty=1,staplewex=0,
          yaxt="n",xaxt="n",ylab="",xlab="")

  #Add x-axis labels
  axis(1,at=(1:length(contigLabels))-0.5,tick=F,line=-0.8,las=2,labels=contigLabels)
  mtext(1,text="Chromosome",line=2.2)

  #Add y-axis labels
  axis(2,at=0:ymax,las=2)
  mtext(2,text="Estimated Copy Number",line=2.2)
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
              help="disable copy number visualization [default: FALSE]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] median_coverage_matrix",
                                option_list=option_list),
                   positional_arguments=TRUE)
INFILE <- args$args[1]
OUTDIR <- args$options$OUTDIR
noplot <- args$options$noplot
if(is.null(OUTDIR)){
  OUTDIR <- ""
}

#Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Must supply an input median coverage matrix\n")
}

#Loads data
dat <- readMatrix(INFILE)

#Transforms data to predicted copy numbers
dat <- normalizeContigsPerMatrix(dat,exclude=24:25,ploidy=2)

#Plots boxplots per contig
png(paste(OUTDIR,"/estimated_CN_per_contig.png",sep=""),
    height=1250,width=2500,res=300)
boxplotsPerContig(dat,exclude=24:25)
dev.off()

#Plots sex assignment dotplot
png(paste(OUTDIR,"/sex_assignments.png",sep=""),
    height=1500,width=1500,res=300)
sexes <- assignSex(dat)
dev.off()







#====================================================
#extract passed parameters
#----------------------------------------------------
args <- commandArgs(TRUE)
file <- args[1]
xMin <- as.numeric(args[2])
xMax <- as.numeric(args[3])
geneEnd <- args[4]
combinedStrands <- as.numeric(args[5])
#====================================================


#====================================================
#get data
#----------------------------------------------------
data <- read.table(file,header=TRUE,sep=",")
data <- data[data$bin > xMin & data$bin < xMax,]
#====================================================


#====================================================
#common plot properties
#----------------------------------------------------
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlab <- paste('Distance from Annotated', geneEnd, '(bp)')
xlim <- c(xMin,xMax)
#====================================================


#====================================================
#fracUV plot
#----------------------------------------------------
if(combinedStrands == 1){
    jpgFile <- paste(file,".fracUV.jpg",sep="")
    ylab <- 'Normalized Fraction UV'
    ylim <- c(0,1)
    bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
    plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    lines(c(0,0),ylim)
    lines(data$bin,data$fracUV50,col=1)
    lines(data$bin,data$fracUV5,col=1,lty=2) #dashed lines
    lines(data$bin,data$fracUV95,col=1,lty=2)
}
#====================================================


#====================================================
#counts plots
#----------------------------------------------------
ylab <- 'Relative Bin Density'
#BLUE = control
#GREEN = UV
green <- rgb(0,0.5,0)
#----------------------------------------------------
#----------------------------------------------------
#medians with 5th and 95th percentile lines
#----------------------------------------------------
nMax <- max(data$nRelDens95,na.rm=TRUE)
uMax <- max(data$uRelDens95,na.rm=TRUE)
yMax <- max(c(nMax,uMax))
nMin <- min(data$nRelDens5,na.rm=TRUE)
uMin <- min(data$uRelDens5,na.rm=TRUE)
yMin <- min(c(0,nMin,uMin))
ylim <- c(yMin,yMax)
jpgFile <- paste(file,".Counts.withRange.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),ylim)
lines(xlim,c(0,0))
strand <- 1
toPlot <- data[data$strand == strand,]
lines(toPlot$bin,toPlot$nRelDens50,col="blue")
lines(toPlot$bin,toPlot$nRelDens5,col="blue",lty=2) #dashed lines
lines(toPlot$bin,toPlot$nRelDens95,col="blue",lty=2)
lines(toPlot$bin,toPlot$uRelDens50,col=green)
lines(toPlot$bin,toPlot$uRelDens5,col=green,lty=2) #dashed lines
lines(toPlot$bin,toPlot$uRelDens95,col=green,lty=2)
if(combinedStrands == 0){
    strand <- -1
    toPlot <- data[data$strand == strand,]
    lines(toPlot$bin,toPlot$nRelDens50,col="blue")
    lines(toPlot$bin,toPlot$nRelDens5,col="blue",lty=2) #dashed lines
    lines(toPlot$bin,toPlot$nRelDens95,col="blue",lty=2)
    lines(toPlot$bin,toPlot$uRelDens50,col=green)
    lines(toPlot$bin,toPlot$uRelDens5,col=green,lty=2) #dashed lines
    lines(toPlot$bin,toPlot$uRelDens95,col=green,lty=2)
}
#----------------------------------------------------
#----------------------------------------------------
#medians without 5th and 95th percentile lines
#----------------------------------------------------
nMax <- max(data$nRelDens50,na.rm=TRUE)
uMax <- max(data$uRelDens50,na.rm=TRUE)
yMax <- max(c(nMax,uMax))
nMin <- min(data$nRelDens50,na.rm=TRUE)
uMin <- min(data$uRelDens50,na.rm=TRUE)
yMin <- min(c(0,nMin,uMin))
ylim <- c(yMin,yMax)
jpgFile <- paste(file,".Counts.noRange.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),ylim)
lines(xlim,c(0,0))
strand <- 1
toPlot <- data[data$strand == strand,]
lines(toPlot$bin,toPlot$nRelDens50,col="blue")
lines(toPlot$bin,toPlot$uRelDens50,col=green)
if(combinedStrands == 0){
    strand <- -1
    toPlot <- data[data$strand == strand,]
    lines(toPlot$bin,toPlot$nRelDens50,col="blue")
    lines(toPlot$bin,toPlot$uRelDens50,col=green)
}
#----------------------------------------------------
#----------------------------------------------------
#means
#----------------------------------------------------
nMax <- max(data$nRelDensMean,na.rm=TRUE)
uMax <- max(data$uRelDensMean,na.rm=TRUE)
yMax <- max(c(nMax,uMax))
nMin <- min(data$nRelDensMean,na.rm=TRUE)
uMin <- min(data$uRelDensMean,na.rm=TRUE)
yMin <- min(c(0,nMin,uMin))
ylim <- c(yMin,yMax)
jpgFile <- paste(file,".Counts.means.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),ylim)
lines(xlim,c(0,0))
strand <- 1
toPlot <- data[data$strand == strand,]
lines(toPlot$bin,toPlot$nRelDensMean,col="blue")
lines(toPlot$bin,toPlot$uRelDensMean,col=green)
if(combinedStrands == 0){
    strand <- -1
    toPlot <- data[data$strand == strand,]
    lines(toPlot$bin,toPlot$nRelDensMean,col="blue")
    lines(toPlot$bin,toPlot$uRelDensMean,col=green)
}
#====================================================


#====================================================
#nGenes plot
#----------------------------------------------------
if(combinedStrands == 1){
    jpgFile <- paste(file,".nGenes.jpg",sep="")
    ylab <- '# of genes contributing to bin'
    yMax <- max(data$nGenes,na.rm=TRUE)
    ylim <- c(0,yMax)
    bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
    plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
    lines(c(0,0),ylim)
    lines(data$bin,data$nGenes,col=1)
}
#====================================================



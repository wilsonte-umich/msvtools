#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]

#get data
data <- read.table(file,header=TRUE,sep=",")

#common plot properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
#pch <- 20
pch <- "."
cex <- 0.3

#size vs. fracSpliced plot
jpgFile <- paste(file,".size.jpg",sep="")
xlab <-'Gene Size (bp)'
ylab <- 'Fraction Spliced'
xmin <- min(data$geneSize,na.rm=TRUE)
xmax <- max(data$geneSize,na.rm=TRUE)
xlim <- c(xmin,xmax)
ylim <- c(0,2)
log <- 'x'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(xmin,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
points(data$geneSize,data$splicedFrac,pch=pch,col=1)


#exon Density vs. intron Density
jpgFile <- paste(file,".density.jpg",sep="")
xlab <-'Exon Density'
ylab <- 'Intron Density'
eMin <- min(data$eDensity,na.rm=TRUE)
eMax <- max(data$eDensity,na.rm=TRUE)
iMin <- min(data$iDensity,na.rm=TRUE)
iMax <- max(data$iDensity,na.rm=TRUE)
xmin <- min(eMin,iMin)
xmax <- max(eMax,iMax)
lim <- c(xmin,xmax)
xlim <- lim
ylim <- lim
log <- 'xy'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(eMin,iMin,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(xmin,xmax),c(xmin,xmax))
points(data$eDensity,data$iDensity,pch=pch,col=1)

#fracSpliced histogram
jpgFile <- paste(file,".fracSpliced.jpg",sep="")
xlab <-'Percent Unspliced'
ylab <- 'Frequency'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
histData <- round(data[data$splicedFrac < 2,'splicedFrac'] * 100)
hist(histData,prob=TRUE,ylab=ylab,xlab=xlab,breaks=100)








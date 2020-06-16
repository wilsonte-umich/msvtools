#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
xMin <- as.numeric(args[2])
xMax <- as.numeric(args[3])
geneEnd <- args[4]

#get data
data <- read.table(file,header=TRUE,sep=",")
data <- data[data$bin > xMin & data$bin < xMax,]

#common plot properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlab <- paste('Distance from Annotated', geneEnd, '(bp)')
xlim <- c(xMin,xMax)

#fracUV plot
jpgFile <- paste(file,".fracUV.jpg",sep="")
ylab <- 'Normalized Fraction UV'
ylim <- c(0,1)
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),c(0,1))
lines(data$bin,data$fracUV50,col=1)
lines(data$bin,data$fracUV5,col=1,lty=2) #dashed lines
lines(data$bin,data$fracUV95,col=1,lty=2)

#counts plot
jpgFile <- paste(file,".Counts.jpg",sep="")
ylab <- 'Relative Bin Density'
nMax <- max(data$nRelDens50,na.rm=TRUE)
uMax <- max(data$uRelDens50,na.rm=TRUE)
yMax <- max(c(nMax,uMax))
ylim <- c(0,yMax)
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),c(0,yMax))
lines(data$bin,data$nRelDens50,col="blue")
lines(data$bin,data$uRelDens50,col=rgb(0,0.5,0))

#nGenes plot
jpgFile <- paste(file,".nGenes.jpg",sep="")
ylab <- '# of genes contributing to bin'
yMax <- max(data$nGenes,na.rm=TRUE)
ylim <- c(0,yMax)
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(0,0),c(0,yMax))
lines(data$bin,data$nGenes,col=1)



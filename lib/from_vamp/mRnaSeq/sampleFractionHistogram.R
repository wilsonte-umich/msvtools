#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
ratioField <- args[2]
inputPath <- args[3]

#get data
data <- read.table(file,header=TRUE,sep=",")

#set derived image properties
jpgFileRoot <- paste(inputPath,'/',ratioField,".sampleRatioHist",sep="")
main <- paste(ratioField,'sample ratio histogram')
maxY <- max(data$Y,na.rm=TRUE)
minY <- 0
ylim <- c(minY,maxY)
#maxX <- max(data$X,na.rm=TRUE)
#minX <- min(data$X,na.rm=TRUE)
maxX <- 3
minX <- 0
xlim <- c(minX,maxX)

#set fixed image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlab <- 'Sample Normalized Density Ratio'
ylab <- 'Frequency'

#linear plot
jpgFile <- paste(jpgFileRoot,".linear.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$X,data$Y,pch=pch,cex=cex,col='red')

#log plot
log <- 'x'
xlim <- c(0.1,maxX)

jpgFile <- paste(jpgFileRoot,".log.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$X,data$Y,pch=pch,cex=cex,col='red')


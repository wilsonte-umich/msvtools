#extract passed parameters
args <- commandArgs(TRUE)
file1 <- args[1]
file2 <- args[2]
sample1 <- args[3]
sample2 <- args[4]
inputPath <- args[5]

#get data
data1 <- read.table(file1,header=TRUE,sep=",")
data2 <- read.table(file2,header=TRUE,sep=",")

#set derived image properties
jpgFile <- paste(inputPath,'/',sample1,'_',sample2,".normDensHist.jpg",sep="")
main <- paste(sample1,sample2,'gene normalized density histogram')
maxY <- max(data1$Y,data2$Y,na.rm=TRUE)
minY <- 0
ylim <- c(minY,maxY)
maxX <- max(data1$X,data2$X,na.rm=TRUE)
minX <- 1/1000
xlim <- c(minX,maxX)

#set fixed image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
log <- 'x'
xlab <- 'normalized density'
ylab <- 'fraction of genes'

#plot
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data1$X,data1$Y,pch=pch,cex=cex,col='red')
lines(data2$X,data2$Y,pch=pch,cex=cex,col='blue')




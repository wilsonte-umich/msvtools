#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
sample1 <- args[2]
sample2 <- args[3]

#set derived image properties
jpgFile <- paste(file,".corr.jpg",sep="")
xlab <- paste(sample1,'normalized density')
ylab <- paste(sample2,'normalized density')
main <- paste(sample1,sample2,'gene correlation plot')

#set fixed image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
maxAxis <- 1000
minAxis <- 1 / maxAxis
xlim <- c(minAxis,maxAxis)
ylim <- c(minAxis,maxAxis)
log <- 'xy'

#get data
data <- read.table(file,header=TRUE,sep=",")

#plot
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(xlim,ylim)
points(data[[sample1]],data[[sample2]],pch=pch,cex=cex)


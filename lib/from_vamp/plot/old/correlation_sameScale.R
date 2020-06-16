#generic script for making correlation plots
#where both axes have the same scale and range

#extract passed parameters
args <- commandArgs(TRUE)
dataFile <- args[1]
main <- args[2]
sample1 <- args[3]
sample2 <- args[4]
log <- args[5]
minSampleValue <- as.numeric(args[6])
if(!(log=="xy")){log=""}

#get data
data <- read.table(dataFile,header=TRUE,sep=",")

#set derived image properties
jpgFile <- paste(dataFile,".corr.jpg",sep="")
xlab <- sample1
ylab <- sample2
minValue <- min(data[[sample1]], data[[sample2]], na.rm=TRUE)
if(log=="xy"){minValue <- max(minSampleValue,minValue)}
maxValue <- max(data[[sample1]], data[[sample2]], na.rm=TRUE)
xlim <- c(minValue,maxValue)
ylim <- c(minValue,maxValue)

#set fixed image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3

#plot
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(xlim,ylim)
points(data[[sample1]],data[[sample2]],pch=pch,cex=cex,col="blue")


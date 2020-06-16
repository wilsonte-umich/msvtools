#generic script for making histogram plots
#assumes input file already has columns val and freq

#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
main <- args[2]
xlab <- args[3]
log <- args[4]
if(!(log=="x")){log=""}

#get data
data <- read.table(file,header=TRUE,sep=",")

#set derived image properties
jpgFile <- paste(file,".hist.jpg",sep="")
maxY <- max(data$freq,na.rm=TRUE)
minY <- 0
ylim <- c(minY,maxY)
maxX <- max(data$val,na.rm=TRUE)
minX <- min(data$val,na.rm=TRUE)
xlim <- c(minX,maxX)

#set fixed image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
ylab <- 'Frequency'

#linear plot
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$val,data$freq,pch=pch,cex=cex,col='blue')



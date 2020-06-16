#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
maxGeneDistance <- as.numeric(args[2])

#get data
data <- read.table(file,header=TRUE,sep=",")
data <- data[data$bin > -maxGeneDistance & data$bin < maxGeneDistance,]

#===================================================
#fracUV plot
#===================================================
jpgFile <- paste(file,".Counts.jpg",sep="")
xlab <- 'Distance past polyA (bp)'
ylab <- 'Relative Bin Density'
main <- 'Hit Counts'
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlim <- c(-maxGeneDistance,maxGeneDistance)
#nMax <- max(data$nCount95,na.rm=TRUE)
#uMax <- max(data$uCount95,na.rm=TRUE)
#ylim <- c(0,max(c(nMax,uMax)))
nMax <- max(data$nCount50,na.rm=TRUE)
uMax <- max(data$uCount50,na.rm=TRUE)
ylim <- c(0,max(c(nMax,uMax)))


bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$nCount50,col=1)
#lines(data$bin,data$nCount5,col=1,lty=2) #dashed lines
#lines(data$bin,data$nCount95,col=1,lty=2)
lines(data$bin,data$uCount50,col=2)
#lines(data$bin,data$uCount5,col=2,lty=2) #dashed lines
#lines(data$bin,data$uCount95,col=2,lty=2)



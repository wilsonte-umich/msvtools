#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
maxGeneDistance <- as.numeric(args[2])
startPadding <- as.numeric(args[3])

#get data
data <- read.table(file,header=TRUE,sep=",")
data <- data[data$bin > -startPadding & data$bin < maxGeneDistance,]

#===================================================
#fracUV plot
#===================================================
jpgFile <- paste(file,".Counts.jpg",sep="")
xlab <- 'Distance from Annontated Start (bp)'
ylab <- 'Relative Bin Density'
main <- 'Hit Counts'
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlim <- c(-startPadding,maxGeneDistance)
nMax <- max(data$nRelDens50,na.rm=TRUE)
uMax <- max(data$uRelDens50,na.rm=TRUE)
ylim <- c(0,max(c(nMax,uMax)))

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$nRelDens50,col=1)
lines(data$bin,data$uRelDens50,col=2)




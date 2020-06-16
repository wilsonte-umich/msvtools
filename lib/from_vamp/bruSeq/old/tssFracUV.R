#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
maxGeneDistance <- as.numeric(args[2])

#get data
data <- read.table(file,header=TRUE,sep=",")

#===================================================
#fracUV plot
#===================================================
jpgFile <- paste(file,".fracUV.jpg",sep="")
xlab <- 'Distance from TSS (bp)'
ylab <- 'Normalized Fraction UV'
main <- 'Fraction of UV Gene Hits'
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlim <- c(-maxGeneDistance,maxGeneDistance)
ylim <- c(0,1)

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$forward50,col=1)
lines(data$bin,data$reverse50,col=2)
lines(data$bin,data$forward5,col=1,lty=2) #dashed lines
lines(data$bin,data$reverse5,col=2,lty=2)
lines(data$bin,data$forward95,col=1,lty=2)
lines(data$bin,data$reverse95,col=2,lty=2)

#===================================================
#nGenes plot
#===================================================
#set derived image properties
jpgFile <- paste(file,".nGenes.jpg",sep="")
ylab <- 'Number of Genes'
main <- 'Number of Genes'
maxForward = max(data$forwardNGenes,na.rm=TRUE) 
maxReverse = max(data$reverseNGenes,na.rm=TRUE) 
ylim <- c(0, max(c(maxForward,maxReverse),na.rm=TRUE))

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$forwardNGenes,col=1)
lines(data$bin,data$reverseNGenes,col=2)



#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
xMin <- as.numeric(args[2])
xMax <- as.numeric(args[3])
geneEnd <- args[4]

#get data
data <- read.table(file,header=TRUE,sep=",")
data <- data[data$bin > xMin & data$bin < xMax,]

#===================================================
#fracUV plot
#===================================================
jpgFile <- paste(file,".Counts.jpg",sep="")
xlab <- paste('Distance from Annontated', geneEnd, '(bp)')
ylab <- 'Relative Bin Density'
main <- 'Hit Counts'
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlim <- c(xMin,xMax)
nMax <- max(data$nRelDens50,na.rm=TRUE)
uMax <- max(data$uRelDens50,na.rm=TRUE)
yMax <- max(c(nMax,uMax))
ylim <- c(0,yMax)

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(c(0,0),c(0,yMax))
lines(data$bin,data$nRelDens50,col="blue")
lines(data$bin,data$uRelDens50,col="red")




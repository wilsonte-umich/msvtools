#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
sample1 <- args[2]
sample2 <- args[3]

#get data
data <- read.table(file,header=TRUE,sep=",")

#set derived image properties
jpgFile <- paste(file,".rank.jpg",sep="")
xlab <- 'rank'
ylab <- paste(sample2,'/', sample1,'ratio')
main <- paste(sample1,sample2,'gene rank ratio plot')
maxX <- max(data$X,na.rm=TRUE)
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
maxY <- 100
minY <- 1 / maxY
ylim <- c(minY,maxY)
log <- 'y'

#plot
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(xlim,c(1,1))
points(data$X,data$Y,pch=pch,cex=cex)




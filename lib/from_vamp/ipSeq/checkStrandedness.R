#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]
sample <- args[2]

#get data
data <- read.table(file,header=TRUE,sep=",")

#set derived image properties
jpgFile <- paste(file,".jpg",sep="")
xlab <- 'genome bin position (bp)'
ylab <- 'bin hit count'
main <- sample
maxTop <- max(data$top, na.rm=TRUE)
maxBottom <- max(abs(data$bottom), na.rm=TRUE)
maxY <- max(c(maxTop,maxBottom), na.rm=TRUE)
ylim <- c(-maxY,maxY)
minX <- min(data$bin, na.rm=TRUE)
maxX <- max(data$bin, na.rm=TRUE)
xlim <- c(minX,maxX)

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
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$top,col='blue')
lines(data$bin,data$bottom,col='red')





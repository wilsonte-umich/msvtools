#extract passed parameters
args <- commandArgs(TRUE)
file <- args[1]

#parameters
minValue <- 10

#get data
data <- read.table(file,header=TRUE,sep=",")

#common plot properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
minCount <- min(data$nCount_sense,data$nCount_antisense,data$uCount_sense,data$uCount_antisense,na.rm=TRUE)
minCount <- max(minValue,minCount)
maxCount <- max(data$nCount_sense,data$nCount_antisense,data$uCount_sense,data$uCount_antisense,na.rm=TRUE)
lim <- c(minCount,maxCount)
xlim <- lim
ylim <- lim
log <- 'xy'

#uvSenseAntisense plot
jpgFile <- paste(file,".uvSenseAntisense.jpg",sep="")
xlab <- 'UV Sense Normalized Count'
ylab <- 'UV Antisense Normalized Count'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(minCount,minCount),c(maxCount,maxCount))
points(data$uCount_sense,data$uCount_antisense,pch=pch)

#sumSum plot
data$uSum <- max(data$uCount_sense - data$uCount_antisense,minCount,na.rm=TRUE)
data$nSum <- max(data$nCount_sense - data$nCount_antisense,minCount,na.rm=TRUE)
jpgFile <- paste(file,".uvNonUV.jpg",sep="")
xlab <- 'UV Sense - Antisense'
ylab <- 'nonUV Sense - Antisense'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
lines(c(minCount,minCount),c(maxCount,maxCount))
points(data$uSum,data$nSum,pch=pch)

#uFracSense plot
jpgFile <- paste(file,".uFracSense.jpg",sep="")
xlab <- 'nonUV Normalized Gene Density'
ylab <- 'UV Fraction Sense'
minND <- min(data$normD,na.rm=TRUE)
maxND <- max(data$normD,na.rm=TRUE)
xlim <- c(minND,maxND)
ylim <- c(0,1)
log <- 'x'
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
points(data$normD,data$uFracSense,pch=pch)




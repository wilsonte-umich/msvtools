#extract passed parameters
args <- commandArgs(TRUE)
coverageFile <- args[1]
fragsFile <- args[2]
jpgFile <- args[3]
root <- args[4]
chrom <- args[5]
start <- args[6]
end <- args[7]

coverageData <- read.table(coverageFile,header=TRUE,sep=",")
fragsData <- read.table(fragsFile,header=TRUE,sep=",")

width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlab <- paste('Distance from Center (bp)')
ylab <- 'Normalized Coverage'
main <- paste(root,": ",chrom,":",start,"-",end,sep="")
xMin <- min(coverageData$pos)
xMax <- max(coverageData$pos)
#yMin <- min(fragsData$y) - 0.1
yMin <- -0.3
yMax <- 1
xlim <- c(xMin,xMax)
ylim <- c(yMin,yMax)
col <- c("black","red","green3","blue","cyan3","magenta","darkgoldenrod1","darkorchid","gray","black","red","green3","blue","cyan3","magenta","darkgoldenrod1","darkorchid","gray")
#lty <- c(1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3)
lty <- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)


bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(c(xMin,xMax),c(-0.1,-0.1),col="gray")
lines(c(xMin,xMax),c(-0.2,-0.2),col="gray")
lines(c(xMin,xMax),c(-0.3,-0.3),col="gray")
lines(c(0,0),c(yMin,yMax))
lines(c(xMin,xMax),c(0,0))
lines(c(xMin,xMax),c(0.5,0.5))
ncol <- ncol(coverageData)
nsamples <- ncol - 1
for (i in 2:ncol){
    j <- i - 1
    lines(coverageData[,c(1,i)], col=col[j], lty=lty[j])
}

points(fragsData$forwardProm, fragsData$y, pch=">", col="black")
points(fragsData$forwardUnq, fragsData$y, pch=">", col="blue")
points(fragsData$reverseProm, fragsData$y, pch="<", col="black")
points(fragsData$reverseUnq, fragsData$y, pch="<", col="red")
points(fragsData$normalRooted, fragsData$y, pch="|", col="black")

sampleNames <- colnames(coverageData[1,2:ncol])
#cols = 1:nsamples
legend(xMin * 1.075, yMax + 0.05, sampleNames, fill=col, cex = 0.75)





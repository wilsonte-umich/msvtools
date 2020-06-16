
#extract passed histogram file name
args <- commandArgs(TRUE)
csvFile <- args[1]

#set curve-fitting parameters
varEst <- 0.1 #estimate sd start value as 10% of mean fragment size
guassian <- COUNT ~ a*(exp(-((SIZE-m)^2)/(2*(s^2)))) #a = amplitude, i.e. max pair count, m = mean fragment size, s = std dev

#recover the passed data
data <- read.table(csvFile, header=TRUE, sep=',')

#estimate curve fit start values from the tip of the peak
maxCount <- max(data$COUNT)
a <- maxCount
m <- median(data[data$COUNT==a, 'SIZE'])
s <- round(m * varEst)

#run the curve fit; iterate twice so that 2nd time uses more optimal data range
for (i in 1:2)  {
    peak <- data[data$SIZE >= (m - 2*s) & data$SIZE <= (m + 2*s),]
    curveFit <- nls(guassian, data=peak, start=list(a=a,m=m,s=s))
    coef <- coef(curveFit)
    a <- round(coef[['a']])
    m <- round(coef[['m']])
    s <- round(coef[['s']])
}

#return the results CalculateStatistics.pl
meanNormal <- m
stdDevNormal <- s
minNormal <- m - 3*s #main peak defined as mean +/- 3 std dev
maxNormal <- m + 3*s
amplitude <- a
cat(paste(meanNormal,stdDevNormal,minNormal,maxNormal,amplitude,"\n"))

#provide a plot of the data
jpgFile <- paste(csvFile,".jpg",sep="")
xlab <- 'Convergent Pair Size'
ylab <- 'Pair Count'
main <- 'Convergent Pairs Histogram Plot'
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3
xlim <- c(0,meanNormal * 3)
ylim <- c(0,maxCount)
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(c(minNormal,minNormal),c(0,maxCount))
lines(c(meanNormal,meanNormal),c(0,maxCount))
lines(c(maxNormal,maxNormal),c(0,maxCount))
points(data$SIZE,data$COUNT,pch=pch,cex=cex,col='blue')





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
xlim <- c(0,maxGeneDistance)
ylim <- c(0,1)

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$fracUV50,col=1)
lines(data$bin,data$fracUV5,col=1,lty=2) #dashed lines
lines(data$bin,data$fracUV95,col=1,lty=2)


#===================================================
#nGenes plot
#===================================================
#set derived image properties
jpgFile <- paste(file,".nGenes.jpg",sep="")
ylab <- 'Number of Genes'
main <- 'Number of Genes'
ylim <- c(0, max(data$nGenes,na.rm=TRUE))

bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(data$bin,data$nGenes,col=1)




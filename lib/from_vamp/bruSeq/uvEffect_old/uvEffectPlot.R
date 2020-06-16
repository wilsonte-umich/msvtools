#functions to derive file names
getCSV <- function(mod='') { getFile('csv', mod) }
getJPG <- function(mod='') { getFile('jpg', mod) }
getFile <- function(ext, mod='') {
	if (mod != '') { paste(fileRoot, '_', mod, '.', ext, sep='') } 
    else { paste(fileRoot, '.', ext, sep='') }
}

#extract passed parameters
args <- commandArgs(TRUE)
minX <- as.numeric(args[1])
maxX <- as.numeric(args[2])
binSize <- args[3]
outputPath <- args[4]
nPrevArgs <- 4
nFiles <- length(args) - nPrevArgs

for (stat in c('MEAN', 'MEDIAN')){

#read data for all inputs
for (i in 1:nFiles){
	fileRoot <- args[i + nPrevArgs]
	dataFile <- getCSV('stats')
	data <- read.table(dataFile,header=TRUE,sep=',')
	data <- data[data$DISTANCE >= minX & data$DISTANCE <= maxX,]
	if (i == 1) {
        nRows <- length(data$DISTANCE)  	    
        toPlot <- data.frame(row.names=c(1:nRows))   
	}
	toPlot[[i]] <- data[[stat]]
}
maxY <- ceiling(max(toPlot, na.rm=TRUE))
toPlot$DISTANCE <- data$DISTANCE  

#plot the data
fileRoot <- paste(outputPath,'/UVEffect_',binSize,sep="")
jpgFile <- getJPG(stat)
bitmap(file=jpgFile,type='jpeg',width=1000,height=1000,units='px') #don't use jpep, server doesn't recognize it
plot(0,0,xlim=c(minX,maxX),ylim=c(0,maxY),main='UV effect on brU gene distribution',xlab='distance from gene start',ylab='relative bin density',pch='')
for (i in 1:nFiles){ lines(toPlot$DISTANCE, toPlot[[i]], col=i)}
graphics.off()

}







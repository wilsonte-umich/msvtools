#functions to derive file names
getCSV <- function(mod='') { getFile('csv', mod) }
getFile <- function(ext, mod='') {
	if (mod != '') { paste(fileRoot, '_', mod, '.', ext, sep='') } 
    else { paste(fileRoot, '.', ext, sep='') }
}

#extract passed parameters
args <- commandArgs(TRUE)
fileRoot <- args[1]

#calculate derived parameters
dataFile <- getCSV()
statsFile <- getCSV('stats')

#perform data calculations
data <- read.table(dataFile,header=TRUE,sep=',')
nRows <- length(data$DISTANCE)
nCols <- length(data)
values <- data[,2:nCols]
output <- data.frame(row.names=c(1:nRows))
output$DISTANCE <- data$DISTANCE
output$MEAN <- apply(values, 1, mean, na.rm=TRUE)
output$MEDIAN <- apply(values, 1, median, na.rm=TRUE)
output$STDEV <- apply(values, 1, sd, na.rm=TRUE)

#output stats table
write.table(output,file=statsFile, row.names=FALSE, sep=',')








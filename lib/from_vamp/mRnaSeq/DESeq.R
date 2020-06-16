#takes data from vamp and runs DESeq as instructed

#extract passed parameters
args <- commandArgs(TRUE)
dataFile <- args[1]
conds <- args[2]
testName <- args[3]
refName <- args[4]
nSamples <- as.numeric(args[5])
FDR <- as.numeric(args[6])
inputPath <- args[7]
fileRoot <- args[8]

#get data
countsTable <- read.csv(dataFile, header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable) <- countsTable$name2
countsTable <- countsTable[,-1]
conds <- strsplit(conds, ',')
conds <- unlist(conds)
conds <- factor(conds)

#run DESeq
library(DESeq)
method = "normal"
if(nSamples == 2){method = "blind"} #no replicates present
cds <- newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds, method=method)
results <- nbinomTest(cds, refName, testName) #thus fold-change = testName/refName, induced genes > 1
meanBaseMean <- mean(results$baseMean,na.rm=TRUE)
results$normMeanCount <- results$baseMean / meanBaseMean

#set common image properties
width <- 7
height <- 7
units <- 'in' #w and h in inches
pointsize <- 12
res <- 300 #dpi
pch <- 20
cex <- 0.3

#set MA image properties
log="x"
xLimit <- 1000
xlim <- c(1/xLimit,xLimit)
yLimit <- 10
ylim <- c(-yLimit,yLimit)
xlab <- 'Normalized Average Adjusted Count (A)'
ylab <- 'log2 Fold Difference (M)'
main <- paste(testName, 'vs.', refName)

#make MA plot
jpgFile <- paste(inputPath, '/plots/MA/', fileRoot, ".MA.jpg",sep="")
bitmap(file=jpgFile,type='jpeg',width=width,height=height,unit=units,pointsize=pointsize,res=res)
plot(1,1,pch="",log=log,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)
lines(xlim,c(0,0),pch=pch,cex=cex)
points(results$normMeanCount,results$log2FoldChange,pch=pch,cex=cex,col=ifelse(results$padj<FDR,"red","black"))

#print results tables
results <- results[,c('id','normMeanCount','baseMeanA','baseMeanB','foldChange','log2FoldChange','padj')]

results$normMeanCount <- round(results$normMeanCount, digits=4)
results$baseMeanA <- round(results$baseMeanA, digits=0)
results$baseMeanB <- round(results$baseMeanB, digits=0)
results$foldChange <- round(results$foldChange, digits=4)
results$log2FoldChange <- round(results$log2FoldChange, digits=3)
colnames(results)[colnames(results)=='baseMeanA'] <- refName
colnames(results)[colnames(results)=='baseMeanB'] <- testName
FDRheader <- paste('FDR_', FDR, sep='')
results[[FDRheader]] <- ifelse(results$padj<FDR,ifelse(results$log2FoldChange>0,1,-1),0)
allFile <- paste(inputPath, '/tables/all/', fileRoot, ".all.csv",sep="")
write.table(results[order(results$log2FoldChange),],file=allFile,quote=FALSE,sep = ",",na="",row.names=FALSE,col.names=TRUE)
outliers <- results[results$padj<FDR,]
outlierFile <- paste(inputPath, '/tables/outliers/', fileRoot, ".outliers.csv",sep="")
write.table(outliers[order(outliers$log2FoldChange),],file=outlierFile,quote=FALSE,sep = ",",na="",row.names=FALSE,col.names=TRUE)


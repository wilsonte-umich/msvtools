
# get passed arguments
LIB_DIR     <- Sys.getenv('LIB_DIR')
R_LIB_DIR   <- Sys.getenv('R_LIB_DIR')
MODEL_NAME  <- Sys.getenv('MODEL_NAME')
DATAFILE    <- Sys.getenv('DATAFILE')
MEAN_1      <- as.numeric(Sys.getenv('MEAN_1'))
STDEV_1     <- as.numeric(Sys.getenv('STDEV_1'))
MAX_CN      <- as.numeric(Sys.getenv('MAX_CN'))
EPROB_FILE  <- Sys.getenv('EPROB_FILE')
TRAIN_RDATA <- Sys.getenv('TRAIN_RDATA')
HMM_FILE    <- Sys.getenv('HMM_FILE')
BIN_FILE    <- Sys.getenv('BIN_FILE')
PLOT_DIR    <- Sys.getenv('PLOT_DIR')
BIN_SIZE    <- as.numeric(Sys.getenv('BIN_SIZE'))

# set plot shared properties
source(paste(LIB_DIR, "plot_common.R", sep="/"))
width <- 5

# calculate/set parameters
CNs <- 1:MAX_CN # will only model these copy number states
nsd <- 3.5      # number of stdev used for plotting

# load dependencies
library('mixtools', lib.loc=R_LIB_DIR)

# load data
message('loading data')
d <- read.table(DATAFILE, header=FALSE, sep="\t")
inCol <- c('CHROM', 'START', 'END', 'READ_COUNT')
colnames(d) <- inCol
nunflt <- nrow(d)

# plot a wide histogram from zero to one above MAX_CN
# helpful for validating that MAX_CN, MEAN_1 and STDEV_1 are well-set
message('plotting wide read count histogram')
JPG_FILE <- paste(PLOT_DIR, "/", MODEL_NAME, ".train.HISTOGRAM_WIDE.jpg", sep="")
jpeg(JPG_FILE, 
       width=width, height=height, unit=units, pointsize=pointsize, res=resol)
dd <- d[d$READ_COUNT<=MEAN_1*(MAX_CN+1) + 3.5*STDEV_1*(MAX_CN+1),]
hist(dd$READ_COUNT, breaks=200)
abline(v=MEAN_1,col="blue")
abline(v=MEAN_1 - nsd*STDEV_1,col="red")
abline(v=MEAN_1 + nsd*STDEV_1,col="red")
abline(v=MEAN_1*MAX_CN + nsd*STDEV_1*MAX_CN,col="red")
graphics.off()

# subset the bins to those with usable copy numbers
# all other bins dropped _forever_ as
#   presumptive CN==0 (can never have a CNV)
#   unmappable bin (always useless)
#   excessively high copy number (could never reliably detect a CNV)
#   highly repetitive region with excessive mapping (always useless)
message('discarding unusable bins')
d   <- d[d$READ_COUNT>=MEAN_1        - nsd*STDEV_1 &
       d$READ_COUNT<=MEAN_1*MAX_CN + nsd*STDEV_1*MAX_CN,]
RCs    <- round(min(d$READ_COUNT),0):round(max(d$READ_COUNT),0)
RCagg  <- aggregate(d$READ_COUNT, list(RC=round(d$READ_COUNT,0)), length)
chroms <- sort(unique(d$CHROM)) # the input chromosomes
nBin <- nrow(d)
message(paste("    ", round(nBin / nunflt * 100, 2), "% of bins kept as usable"))

# downsample further to speed up mixtools
message('solving mixed copy number model')
ns   <- min(nBin, 1e5)
dd   <- d[sample(nrow(d), ns),]

# run mixtools to probabilistically correlate read counts to copy numbers
mdl <- normalmixEM(
    dd$READ_COUNT,
    mu          = CNs * MEAN_1,
    sigma       = CNs * STDEV_1,
    mean.constr = paste(CNs, 'a', sep=""), # constrain mean read count to be quantal
    sd.constr   = NULL                     # but not sd (doesn't converge if constrained...)
)
CNs_mdl <- data.frame(CN=CNs, frac=mdl$lambda, mean=mdl$mu, sd=mdl$sigma)
CNs_mdl

# make a plot depicting the final model
message('plotting final model')
JPG_FILE <- paste(PLOT_DIR, "/", MODEL_NAME, ".train.MODEL.jpg", sep="")
jpeg(JPG_FILE, 
       width=width, height=height, unit=units, pointsize=pointsize, res=resol)
hist(d$READ_COUNT, breaks=200, freq=FALSE) # histogram of all usable bins, not the downsample
dns <- matrix(NA, length(RCs), MAX_CN)
message(paste("model parameters, including bins from CN 1 to", MAX_CN))
for (CN in CNs) {
    y <- mdl$lambda[CN] * dnorm(RCs, mdl$mu[CN], mdl$sigma[CN])
    lines(RCs, y, lwd=2, col=CN)
    dns[,CN] <- y
}
lines(RCs, apply(dns, 1, sum), lty=2, lwd=2, col="red")
graphics.off()

# calculate all possible emission probabilities
message("calculating emission probabilities for unique read count/CN combinations")
ep         <- list() # hold calculated probabilities for lookup
for(CN in CNs){ 
    message(paste("  ", CN))
    ep[[CN]] <- list()
    for(RC in RCs){
        ep[[CN]][[RC]] <- suppressWarnings(
            max(1e-99, # ensure non-zero probabilities in all states
                diff( pnorm(RC+c(-1/2, 1/2), mdl$mu[CN], mdl$sigma[CN]) ),
                na.rm=TRUE)
        )  
    }
}

# associate emission probabilities with bins
message("assigning emission probabilities to individual bins")
setEP <- function(CN, RC) ep[[CN]][[round(RC,0)]]
for(CN in CNs){
    lbl <- paste('CN', CN, sep="")
    message(paste("  ", lbl))
    d[[lbl]] <- mapply(setEP, CN, d$READ_COUNT)
}

# run the HMM on each chromosome
message("solving HMM")
d2 <- data.frame()
d$CN <- 'CN_NA'
CN_lbls <- paste('CN', CNs, sep="")
for (chrom in chroms){
    message(paste("  ", chrom))
    dd <- d[d$CHROM==chrom,]
    write(do.call(paste, c(dd[,CN_lbls], sep=",")), file=EPROB_FILE, ncolumns=1)
    pipe=paste("cat ", EPROB_FILE, " | ", LIB_DIR, "/segment -z 0.1 -p 0.995", sep="")   # run Viterbi
    cni <- read.table(pipe(pipe), header=FALSE, sep="\t")[,1] + 1        # get output state calls
    dd$CN <- CNs[cni]                                                    # simplify to copy numbers                                       
    d2 <- rbind(d2, dd)                                                  # assemble the updated table
}
unlink(EPROB_FILE)
CNs_hmm <- aggregate(d2$CN, list(CN=d2$CN), length)
CNs_hmm

# plot the modeled CN to read count correlation
message("plotting boxplot of HMM output")
JPG_FILE <- paste(PLOT_DIR, "/", MODEL_NAME, ".train.HMM_BOXPLOT.jpg", sep="")
jpeg(JPG_FILE, 
       width=width, height=height, unit=units, pointsize=pointsize, res=resol)
boxplot(READ_COUNT ~ CN, d2)
graphics.off()

# save training data for use by refine
message("saving RData for 'refine'")
save(ep, RCagg, CNs_mdl, CNs_hmm, file=TRAIN_RDATA)
#load(TRAIN_RDATA)

# determine runs of contiguous bins (faster way to do this?)
message("finding segments of contiguous copy number states")
segN  <- 0
d2$segN <-
    sapply( # break on any copy number change OR bin discontinuity
        c(TRUE, d2[2:nrow(d2),'END'] -  d2[1:(nrow(d2)-1),'END'] != BIN_SIZE) |
        c(TRUE, d2[2:nrow(d2),'CN']  != d2[1:(nrow(d2)-1),'CN']),
    function(inc){
        if(inc) segN <<- segN + 1
        segN
    })

# collapse the segments and determine aggregate statistics
message("aggregating segment parameters")
seg <- data.frame(segN=sort(unique(d2$segN)))
by <- list(segN=d2$segN)
seg$CHROM   <- aggregate(d2$CHROM, by, function(x) x[1] )[[2]]
seg$START   <- aggregate(d2$START, by, min)[[2]] - 1
seg$END     <- aggregate(d2$END,   by, max)[[2]]
seg$NAME    <- "."
seg$CN      <- aggregate(d2$CN, by, function(x) x[1] )[[2]]
seg$STRAND  <- "+"
seg$N_BIN   <- aggregate(d2$READ_COUNT,   by, length)[[2]]
seg$RC_MEAN <- aggregate(d2$READ_COUNT,   by, mean)[[2]]
seg$RC_SD   <- aggregate(d2$READ_COUNT,   by, sd)[[2]]

# determine the likelihood of each CN state for each segment
message("calculating segment CN log likelihoods")
for(CN in CNs){
    lbl <- paste('CN', CN, sep="")
    message(paste("  ", lbl))
    Ls <- list()
    seg[[lbl]] <- aggregate(d2$READ_COUNT, by, function(RCs){
        sd <- round(ifelse(length(RCs)==1, mdl$sigma[CN], max(mdl$sigma[CN], sd(RCs))), 0)
        sum(sapply(RCs, function(RC) {
            RC <- round(RC,0)
            key <- paste(RC, sd, sep=",")
            if(!(key %in% names(Ls))){
                Ls[[key]] <- log(suppressWarnings(
                    max(0,
                        diff( pnorm(RC+c(-1/2, 1/2), mdl$mu[CN], sd) ),
                        na.rm=TRUE)
                ))
            }
            Ls[[key]]
        } ))
    })[[2]]
}
Ls <- list()

# convert to relative likelihoods
message("converting to relative likelihoods, i.e. weights")
maxLL <- sapply(1:nrow(seg), function(i) max(seg[i,CN_lbls]))
for(CN in CNs){
    lbl <- paste('CN', CN, sep="")
    message(paste("  ", lbl))    
    seg[[lbl]] <- mapply(function(LL, maxLL) exp(LL - maxLL), seg[[lbl]], maxLL)
}

# write the final output BED file
message("writing segments BED")
pipe <- paste("cut -f2- | awk 'BEGIN{OFS=\"\\t\"}{$2+=0;$3+=0;print $0}' | bgzip -c > ", HMM_FILE, sep="")
write.table(
    seg,
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)

# write the final output BED file
message("writing bins BED")
d2$START <- d2$START - 1
pipe <- paste("awk 'BEGIN{OFS=\"\\t\"}{$2+=0;$3+=0;print $0}' | bgzip -c > ", BIN_FILE, sep="")
write.table(
    d2[,c('CHROM', 'START', 'END', 'CN', 'READ_COUNT')],
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)

message("R done")

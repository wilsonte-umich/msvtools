
# this script coordinates the HMM of 'refine' and 'segment'

# collect parameters
LIB_DIR      <- Sys.getenv('LIB_DIR')
MODEL_NAME   <- Sys.getenv('MODEL_NAME')  # e.g. a cell line or an individual sample
DATAFILE     <- Sys.getenv('DATAFILE')    # input data, LRR, etc.
PLOT_DIR     <- Sys.getenv('PLOT_DIR')
HMM_FILE     <- Sys.getenv('HMM_FILE')    # output file
TRAIN_RDATA  <- Sys.getenv('TRAIN_RDATA') # for refine, optional data from 'train'
PLOT_PREFIX  <- paste(PLOT_DIR, MODEL_NAME, sep="/")

# set parameters
extreme_zyg    <- 0.85         # more extreme zygosity/BAF is uninformative to copy number
lrr_lim        <- c(-1.5, 1.0) # range of usable LRR values
BAD_PROBE_FREQ <- 1e-4         # assume this many probes are wacky and entirely unpredictable
PERSISTENCE    <- 1 - 1e-4     # reciprocal of HMM transition probability

# load dependencies (don't remember what "ot" stood for! its not the output table...)
ot <- read.table(paste(LIB_DIR, "/state_table.txt", sep=""), header=TRUE, sep="\t")

# load the data
message("loading data")
message(DATAFILE)

####################
d <- read.table(DATAFILE, header=FALSE, sep="\t")
colnames(d) <- c('CHROM','START','POS','NAME','CN_OUT','STRAND','RC','CN_IN','LRR','BAF')

#d <- read.table(DATAFILE, header=TRUE, sep="\t")

isBam <- file.exists(TRAIN_RDATA)
if(isBam) load(TRAIN_RDATA)
hasInf <- "INF" %in% colnames(d)

# analyze the data
message("analyzing the data")
d$ZYG      <- abs(d$BAF - 0.5) + 0.5 # convert B allele frequency to zygosity
if(!hasInf) d$INF <- ifelse(d$ZYG < extreme_zyg, 1, 0) # inherited from refined model for 'segment'
pInf       <- nrow(d[d$INF==1,]) / nrow(d) # fraction of informative probes over entire array
dlrr       <- d$LRR >= lrr_lim[1] & d$LRR <= lrr_lim[2] # probes useful for estimating ranges
CNs        <- sort(unique(d$CN_IN)) # the available starting model copy numbers
chroms     <- sort(unique(d$CHROM)) # the input chromosomes

# set common image properties
width     <- 2.5
height    <- 2.8
units     <- 'in' # w and h in inches
pointsize <- 8
resol     <- 900 #dpi
pch       <- "." #20
cex       <- 1 #0.3

# find zygosity statistics across array
cn_col <- 'CN_IN'
source(paste(LIB_DIR, "model_zygosity.R", sep="/"))

# use zygosity statistics to find LRR statistics across array
source(paste(LIB_DIR, "model_LRR.R", sep="/"))

# plot the ZYG-LRR correlation
boxes  <- TRUE
source(paste(LIB_DIR, "plot_ZYG_LRR_correlation.R", sep="/"))

# bin ZYG and LRR values as observation states
message("binning observation states")
n_zyg_bins <- 100 # per unit, i.e. # bins from 0 to 1
n_lrr_bins <- 50  # 100 gives 51 unique ZYG bins from 0.5 to 1
zyg_inc <- 1/n_zyg_bins / 2
lrr_inc <- 1/n_lrr_bins / 2
d$zyg_bin  <- round(d$ZYG * n_zyg_bins, 0) / n_zyg_bins
d$lrr_bin  <- round(d$LRR * n_lrr_bins, 0) / n_lrr_bins

# calculate all possible emission probabilities
# this first pass is done using fit Gaussian distributions
message("calculating emission probabilities for unique ZYG/LRR combinations")
b          <- unique(d[,c('zyg_bin','lrr_bin')]) # all existing combinations of binned ZYG/LRR
uot        <- ot[ot$zyg_mean!=-88,] # usable output states that exist in model
asts       <- as.vector(unique(uot$allelic_state))
astcns     <- sapply(asts, function(ast) ot[ot$allelic_state==ast & ot$allelism=='monoallelic','copy_number'] )
rc_ep      <- if(isBam) ep # reassign training read count emission probabilities
ep         <- list() # hold calculated probabilities for lookup
for(allelic_state in asts){ # one probability for all output CN states for all input combinations
    message(paste("  ", allelic_state))    
    otrs <- ot[ot$allelic_state==allelic_state,]
    ep[[allelic_state]] <- list()
    for(i in 1:nrow(b)){
        cmb <- paste(b[i,'zyg_bin'], b[i,'lrr_bin'], sep=",")
        ep[[allelic_state]][[cmb]] <- list()   


               
        for (INF in 0:1){                          
                  otr <- otrs[otrs$informative==INF,]  # choose ep based on probe informativity
                  ep[[allelic_state]][[cmb]][[INF]] <- suppressWarnings(
                      max(BAD_PROBE_FREQ, # ensure non-zero probabilities in all states
                          otr$prob_scalar * # edge ZYG distributions are twice the Gaussian probabilities since symmetric
                              diff( pnorm(b[i,'zyg_bin']+c(zyg_inc, -zyg_inc), otr$zyg_mean, otr$zyg_stdev) ) *
                              diff( pnorm(b[i,'lrr_bin']+c(lrr_inc, -lrr_inc), otr$lrr_mean, otr$lrr_stdev) ),
                          na.rm=TRUE)
                  )                   
        }
    }
}


asts
astcns
quit("no")






# associate emission probabilities with probes
message("assigning emission probabilities to individual probes")
setEP <- function(allelic_state, CN, RC, zyg_bin, lrr_bin, INF){     # binned ZYG and LRR values
    cmb <- paste(zyg_bin, lrr_bin, sep=",")
    ep[[allelic_state]][[cmb]][[INF]] * if(isBam){
        rc_ep[[CN]][[round(RC,0)]]
    } else {
        1
    }
}
for(i in 1:length(asts)){
    message(paste("  ", asts[i]))    
    d[[asts[i]]] <- mapply(setEP, asts[i], astcns[i], d$RC, d$zyg_bin, d$lrr_bin, d$INF)
}



# run the HMM on each chromosome
message("solving HMM")
d2 <- data.frame()
d$allelic_state <- 'CN_NA'
EPROB_FILE <- paste(HMM_FILE, ".eprob.txt", sep="")
for (chrom in chroms){
    message(paste("  ", chrom))
    dd <- d[d$CHROM==chrom,]
    write(do.call(paste, c(dd[,asts], sep=",")), file=EPROB_FILE, ncolumns=1)
    pipe=paste("cat ", EPROB_FILE, " | segment -z 0.1 -p ", PERSISTENCE, sep="")  # run Viterbi
    asti <- read.table(pipe(pipe), header=FALSE, sep="\t")[,1] + 1        # get output state calls
    dd$allelic_state <- asts[asti]
    dd$CN_OUT <- astcns[asti]                                             # simplify to copy numbers
    d2 <- rbind(d2, dd)                                                   # assemble the updated table
}
unlink(EPROB_FILE)
d <- d2
rm(d2)

# plot the revised ZYG-LRR correlation (using new CN calls)
cn_col <- 'CN_OUT'
boxes  <- FALSE
source(paste(LIB_DIR, "plot_ZYG_LRR_correlation.R", sep="/"))

# write the final output BED file
pipe <- paste("bgzip -c > ", HMM_FILE, sep="")
write.table(
    
    ##############              
    d[,],
    
    
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)

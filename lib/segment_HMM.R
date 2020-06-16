
# this script coordinates the HMM of 'refine' and 'segment'
# differences between these actions handled inline

# collect parameters
LIB_DIR       <- Sys.getenv('LIB_DIR')
MODEL_NAME    <- Sys.getenv('MODEL_NAME')  # e.g. a cell line or an individual sample
DATAFILE      <- Sys.getenv('DATAFILE')    # input data, LRR, etc.
PLOT_DIR      <- Sys.getenv('PLOT_DIR')
HMM_FILE      <- Sys.getenv('HMM_FILE') 
PROBES_FILE   <- Sys.getenv('PROBES_FILE')   
SEGMENTS_FILE <- Sys.getenv('SEGMENTS_FILE')
SVS_FILE      <- Sys.getenv('SVS_FILE')
TRAIN_RDATA   <- Sys.getenv('TRAIN_RDATA') # for refine, optional data from 'train'
PROBE_COL     <- Sys.getenv('PROBE_COL')
SV_COL        <- Sys.getenv('SV_COL')
PLOT_PREFIX   <- paste(PLOT_DIR, MODEL_NAME, sep="/")
PROBE_COL     <- strsplit(PROBE_COL, ",")[[1]]
SV_COL        <- strsplit(SV_COL, ",")[[1]]

# set parameters
extreme_zyg    <- as.numeric(Sys.getenv('extreme_zyg'))    # more extreme zygosity/BAF is uninformative to copy number
lrr_lim        <- as.numeric(strsplit(Sys.getenv('lrr_lim'), ',')[[1]]) # range of LRR values used for model fitting
zyg_lim        <- c(0.5, 1.0)  # used for plotting
BAD_PROBE_FREQ <- as.numeric(Sys.getenv('BAD_PROBE_FREQ')) # assume this fraction of probes are wacky and entirely unpredictable
PERSISTENCE    <- as.numeric(Sys.getenv('PERSISTENCE'))    # reciprocal of HMM transition probability
PRESERVATION   <- as.numeric(strsplit(Sys.getenv('PRESERVATION'), ',')[[1]])
names(PRESERVATION) <- c('D1','D2','D3','D4')
#PRESERVATION   <- c(D1=0.99,   # penalities (multipliers) for CN changes from model        
#                    D2=0.95,   #     index = CN delta
#                    D3=0.90,   #     gives preference to preserving the input model in the output when segmenting
#                    D4=0.80)

# load dependencies (don't remember what "ot" stood for! its not the output table...)
ot <- read.table(paste(LIB_DIR, "/state_table.txt", sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# load the data
message("loading data")
message(DATAFILE)
d <- read.table(DATAFILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)
nPrbIn <- nrow(d)
message(paste("   ", nPrbIn, "probes in input model"))
isBam <- file.exists(TRAIN_RDATA)
if(isBam) {
    message("loading training emission probabilities")
    #load(TRAIN_RDATA)
    d <- d[!is.na(d$RC),] # only use probes with read count when NGS training is available
    nPrbIn <- nrow(d)
    message(paste("   ", nPrbIn, "probes with read count data will be used"))
}

# sort probes; raw array is NOT necessarily sorted
message("sorting probes")
d <- d[order(d$CHROM, d$POS, d$PRB_NAME),] # yes, there can be duplicate SNP positions over >1 probe

# set up input model parameters
message("identifying segmentation type")
modeling <- !("INF_IN" %in% colnames(d)) # don't already know the presumed probe informativity
if(modeling) {                           # and in that case assume all probes might be informative
    d$INF_IN <- 'NA'
    d$AS_IN  <- 'NA'
}
                
# identify probes with usable data
# must be present in model AND have valid numerical LRR and BAF values
message("determining probe usability")
d$IS_NA  <- ifelse(is.na(d$CN_IN) | is.na(d$LRR) | is.na(d$BAF), 1, 0) # use integer boolean since in BED6 score field       
nUsable  <- nrow(d[d$IS_NA==0,])
message(paste("   ", nUsable, "probes with valid LRR/BAF signal will be used"))

# apply model corrections when segmenting
if(modeling) {
    d$LRR_RAW <- NA
    d$BAF_RAW <- NA    
} else {
    message("applying model corrections")
    source(paste(LIB_DIR, "apply_model_corr.R", sep="/"))
    lrr_lim <- c(-2.5, 2.5)
    zyg_lim <- c(0.5,  1.1)
}

# analyze the data
# d$LRR and d$BAF now contain adjusted values (XXX_RAW contain raw values)
message("calculating derivative values")
chroms     <- sort(unique(d$CHROM))        # the input chromosomes
d$ZYG      <- abs(d$BAF - 0.5) + 0.5       # convert B allele frequency to zygosity
d$INF_OUT  <- ifelse(d$ZYG < extreme_zyg, 'yes', 'no') # informativity based on array data
dlrr       <- d$IS_NA==0 & d$LRR >= lrr_lim[1] & d$LRR <= lrr_lim[2] # probes used for model fitting based on data

# fill the table of informative probabilities
# i.e. the probability of finding an informative probe state in an allelic state
pInf <- nrow(d[d$IS_NA==0 & d$INF_OUT=='yes',]) / nUsable # fraction of informative probes over entire array OR EXTERNAL VALUE?
pINF <- list()
for (allelism in as.character(0:2)){
    pINF[[allelism]] <- list()
    for (INF in c('yes','no')){
        pINF[[allelism]][[INF]] <-
            if(allelism=='0'){
                0.5  
            } else if(allelism=='1'){
                ifelse(INF=='yes', BAD_PROBE_FREQ, 1-BAD_PROBE_FREQ)
            } else {
                ifelse(INF=='yes', pInf,           1-pInf)   
            }  
    }
}

# bin ZYG and LRR values as observation states, using integer bins
message("binning observation states")
n_zyg_bins   <- 150 # per unit, i.e. bins from 0 to 1
n_lrr_bins   <- 30  # 100 gives 51 unique ZYG bins from 0.5 to 1
zyg_bins     <- round(0.5 * n_zyg_bins, 0) : round(max(d$ZYG, na.rm=TRUE) * n_zyg_bins, 0) # ranges of modeled/observed bins
lrr_bins     <- round(min(d$LRR, na.rm=TRUE) * n_lrr_bins, 0) : round(max(d$LRR, na.rm=TRUE) * n_lrr_bins, 0)
d$zyg_bin    <- round(d$ZYG * n_zyg_bins, 0) # the actual bins for each probe observation
d$lrr_bin    <- round(d$LRR * n_lrr_bins, 0)
zyg_bins_obs <- unique(d[d$IS_NA==0,'zyg_bin']) # the list of bins actually observed, for HMM
lrr_bins_obs <- unique(d[d$IS_NA==0,'lrr_bin'])
lrr_zyg_obs  <- unique(d[d$IS_NA==0,c('lrr_bin','zyg_bin')]) # and actually observed pairs
bad_zyg_prob <- 1 / (n_zyg_bins * (1 - 0.5)) # bad probes have equal probability for any state
bad_lrr_prob <- 1 / (n_lrr_bins * (lrr_lim[2] - lrr_lim[1]))
bad_prb_prob <- bad_zyg_prob * bad_lrr_prob
zyg_inc      <- 1/2
lrr_inc      <- 1/2

# fit and plot statistics across array
cn_col <- 'CN_IN'
boxes  <- TRUE
source(paste(LIB_DIR, "model_array.R", sep="/"), local=TRUE, print.eval=TRUE)

#save(lrr_eps, zyg_eps, frc_zyg, ot, file='/treehouse/path-wilsonte-turbo/data/Glover/msvtools/experiments/mosaic.RData')
#quit("no")

#RD_FILE <- '/home/wilsonte_lab/club_house/data/Glover/Irene/samples/msvtools.segment.probes.Scr1-1.bed.bgz.RData'
#load(RD_FILE)
#source(paste(LIB_DIR, "run_HMM.R", sep="/"))

# initialize HMM object
message("initializing HMM")
uot     <- ot[ot$zyg_mean!=-88,] # usable output states that exist in model
asts    <- as.vector(unique(uot$allelic_state))
astcns  <- sapply(asts, function(ast) min(ot[ot$allelic_state==ast,'copy_number']) )
astasms <- sapply(asts, function(ast) min(ot[ot$allelic_state==ast,'allelism']) )
asmscms <- sapply(asts, function(ast)     ot[ot$allelic_state==ast,'mosaic'][1] )
source(paste(LIB_DIR, "run_HMM.R", sep="/"))
if(modeling) {
    ot_cat <- list(AS_IN='NA', INF_IN='NA')    
} else {
    ot_cat <- list(AS_IN=asts[asmscms=="no"], INF_IN=c('yes','no')) # AS_IN can't ever include mosaics 
}
HMM <- HMM_init(
    ot=ot_cat,
    os=list(INF_OUT=c('yes','no'), LRR=lrr_bins_obs, ZYG=zyg_bins_obs),
    hs=list(AS=asts),
    p=PERSISTENCE
)

# calculate all required emission probabilities
# for all combinations of ot, os and hs
message("calculating emission probabilities")
frc_zyg_names <- names(frc_zyg)
os_i <- list() # list of actually observed output state indices
for(INF_OUT in c('yes','no')) {
    os_i[[INF_OUT]] <- HMM$os$i[paste(INF_OUT, lrr_zyg_obs$lrr_bin, lrr_zyg_obs$zyg_bin, sep=",")]
}
# the output hidden states
for(AS_OUT in asts){
    CN_OUT       <- astcns[[AS_OUT]]
    CN_OUT_      <- as.character(CN_OUT)
    allelism_OUT <- astasms[[AS_OUT]]
    mosaic_OUT   <- asmscms[[AS_OUT]]
    message(paste("   ", AS_OUT, " [", CN_OUT_, ":", allelism_OUT, "]", sep=""))
    pinf <- pINF[[as.character(allelism_OUT)]]
    leps <- lrr_eps[[CN_OUT_]]
    leps <- sapply(lrr_zyg_obs$lrr_bin, function(LRR) leps[leps$bin==LRR,'prob'])
    hs_i <- HMM$hs$i[[AS_OUT]]
    frc_zyg_out <- CN_OUT_ %in% frc_zyg_names
    nullallelic <- allelism_OUT == 0
    monoallelic <- allelism_OUT == 1
    # the observation types
    for(AS_IN in ot_cat$AS_IN){ 
        adj_AS_OUT <- 1 # adjust likelihoods of output states based on known modeled states when segmenting
        adj_CN_OUT <- 1
        if(!modeling){
            CN_IN       <- astcns[[AS_IN]]
            allelism_IN <- astasms[[AS_IN]]
            if(allelism_OUT > allelism_IN) adj_AS_OUT <- 0 # known monoallelic can never gain informativity             
            adj_CN_OUT <- if(CN_IN == 0 & CN_OUT >= 1){ # DNA can never become "unmissing" 
                0 
            } else if (CN_OUT == CN_IN){
                if(allelism_OUT == allelism_IN){  # no SV state most favored, weight fully
                    1            
                } else {
                    PRESERVATION[2]               # downgrade CN neutral LOH; double-hit since both a loss and a gain
                }
            } else if (mosaic_OUT == 'yes'){    
                PRESERVATION[1]                   # downgrade mosaic changes similar to CNC=1                
            } else {
                PRESERVATION[abs(CN_OUT - CN_IN)] # downgrade CN changes in a graded fashion
            }
        }
        ep_lev2 <- adj_AS_OUT * adj_CN_OUT         
        for(INF_IN in ot_cat$INF_IN){
            ot_i <- HMM$ot$i[[paste(AS_IN, INF_IN, sep=",")]]       
            zeps <- if(nullallelic){ # output state makes no prediction at all about zygosity
                zyg_eps[["0"]]
            } else if(monoallelic |  # output state only emits to non-informative observation states
                      INF_IN=='no'){ # known uninformative probes only emit non-informative observation states (cannot become informative)
                zyg_eps[["1"]]  
            } else if (frc_zyg_out){ # informative probes in biallelic states (should stay informative when segmenting)
                zyg_eps[[CN_OUT_]]   # the actual zygosity distribution, used for single-peak distributions        
            } else {
                zyg_eps[[AS_OUT]] # modeled zygosity distributions, used for multi-peak distributions and mosaics
            }
            # the observation states 
            for (INF_OUT in c('yes','no')){            
                if(modeling & INF_OUT=='no') zeps <- zyg_eps[["1"]] # first discovering that a probe is uniformative               
                ep_lev3 <- ep_lev2 * # probability of a probe being informative allows modeling monoallelic vs. biallelic states (e.g. ROH)
                           if(modeling) pinf[[INF_OUT]] else 1 # informativity is deterministic when segmenting         
                HMM$ep[ot_i,hs_i,os_i[[INF_OUT]]] <- mapply(function(lep, ZYG){
                    log(ep_lev3 * lep * zeps[zeps$bin==ZYG,'prob'])
                }, leps, lrr_zyg_obs$zyg_bin) 
            }            
        }   
    }
}

# save HMM for future use by compare
message("saving HMM for future use")
if(!modeling) HMM <- HMM_set_median_prob(HMM)
save(HMM, file=HMM_FILE)

# determine observation types and states for each probe
message("setting observation indices for probes")
d$ot <- NA
d$os <- NA
d[,c('ot','os')] <- as.data.frame(HMM_collapse_obs(
    HMM,
    ot=d[,c('AS_IN','INF_IN')],
    os=d[,c('INF_OUT','lrr_bin','zyg_bin')]
))

# run the HMM on each chromosome
message("solving HMM")
d$hs     <- NA
d$AS_OUT <- NA
d$CN_OUT <- NA
for (chrom in chroms){ 
    message(paste("  ", chrom))
    flt <- which(d$IS_NA==0 & d$CHROM==chrom)
    hsi <- HMM_viterbi(HMM, d[flt,c('ot','os')])$hsi
    d[flt,'hs']     <- hsi   
    d[flt,'AS_OUT'] <- asts  [hsi]
    d[flt,'CN_OUT'] <- astcns[hsi]
}

# update informativity based on model output
# probes in monoallelic runs of homozygosity are never considered informative
# probes in a biallelic state have informativity determined by observed zygosity
message("calculating derived output columns")
d$INF_OUT <- mapply(function(IS_NA, AS, INF){
    if(IS_NA==1){
        NA
    } else if(astasms[[AS]]<=1){
        'no'
    } else {
        INF
    }
}, d$IS_NA, d$AS_OUT, d$INF_OUT)

# calculate changes in output relative to input modeled states
d$CNC <- if(modeling) NA else ifelse(d$IS_NA, NA, d$CN_OUT - d$CN_IN)
d$LOH <- if(modeling) NA else ifelse(d$IS_NA, NA, ifelse(d$INF_OUT != d$INF_IN, 'yes', 'no'))

# fit and plot refined statistics across array
cn_col <- 'CN_OUT'
source(paste(LIB_DIR, "model_array.R", sep="/"))

# write the final output probes BED file
message("writing probes BED")
header <- paste(PROBE_COL, collapse="\\t")
header <- paste("#", header, sep="") # comment the header for bgzip
pipe <- paste("awk 'BEGIN{OFS=\"\\t\";print \"", header,"\"}{$2+=0;$3+=0;print $0}' | bgzip -c > ",
              PROBES_FILE, sep="")
write.table(
    d[,PROBE_COL],
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)
nPrbOut <- nrow(d)
message(paste("   ", nPrbOut, "probes written to disk"))

# write the final output structural variants file as needed
if(!modeling) source(paste(LIB_DIR, "find_SVs.R", sep="/"))


#########################
#PROBES_FILE <- '/home/wilsonte_lab/club_house/data/Glover/Irene/samples/msvtools.segment.probes.Scr1-1.bed.bgz'
#RD_FILE <- paste(PROBES_FILE, "RData", sep=".")
##save.image(file=RD_FILE)
##load(RD_FILE)
#source(paste(LIB_DIR, "run_HMM.R", sep="/"))




 
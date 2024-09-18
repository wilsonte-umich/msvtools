# this script coordinates the HMM of 'refine' and 'segment'
# differences between these actions handled inline

#------------------------------------------------------------------
# initialize script
#------------------------------------------------------------------
options(warn=2)
SEGMENTING_ACTION <- Sys.getenv('SEGMENTING_ACTION')
LIB_DIR       <- Sys.getenv('LIB_DIR')
R_LIB_DIR     <- Sys.getenv('R_LIB_DIR')
MODEL_NAME    <- Sys.getenv('MODEL_NAME')   # e.g. a cell line or an individual sample
DATAFILE      <- Sys.getenv('DATAFILE')     # input data, LRR, etc.
GC_BIAS_FILE  <- Sys.getenv('GC_BIAS_FILE') # written by `refine`, read by `segment`
PLOT_DIR      <- Sys.getenv('PLOT_DIR')
HMM_FILE      <- Sys.getenv('HMM_FILE')
PROBES_FILE   <- Sys.getenv('PROBES_FILE')
SEGMENTS_FILE <- Sys.getenv('SEGMENTS_FILE')
SVS_FILE      <- Sys.getenv('SVS_FILE')
PROBE_COL     <- Sys.getenv('PROBE_COL')
SV_COL        <- Sys.getenv('SV_COL')
PLOT_PREFIX   <- file.path(PLOT_DIR, MODEL_NAME)
PROBE_COL     <- strsplit(PROBE_COL, ",")[[1]]
SV_COL        <- strsplit(SV_COL, ",")[[1]]

# load dependencies
suppressMessages(library('data.table', lib.loc=R_LIB_DIR))
suppressMessages(library('R.utils',    lib.loc=R_LIB_DIR))

# set plot shared properties
source(file.path(LIB_DIR, "plot_common.R"))

# set parameters
extreme_zyg    <- as.numeric(Sys.getenv('extreme_zyg'))    # more extreme zygosity/BAF is uninformative to copy number
lrr_lim        <- as.numeric(strsplit(Sys.getenv('lrr_lim'), ',')[[1]]) # range of LRR values used for model fitting
zyg_lim        <- c(0.5, 1.0)  # used for plotting
BAD_PROBE_FREQ <- as.numeric(Sys.getenv('BAD_PROBE_FREQ')) # assume this fraction of probes are wacky and entirely unpredictable
PERSISTENCE    <- as.numeric(Sys.getenv('PERSISTENCE'))    # reciprocal of HMM transition probability
PRESERVATION   <- as.numeric(strsplit(Sys.getenv('PRESERVATION'), ',')[[1]]) # copy number change penalties, CNC1, CNC2, CNC3, CNC4, default 0.99,0.95,0.9,0.8
maxAllowedCN <- 5 # highest CN value available in state_table.txt

#------------------------------------------------------------------
# load and preprocess the array data
#------------------------------------------------------------------
message("loading data")
message(DATAFILE)
d <- fread(DATAFILE, header=TRUE, sep="\t")
nPrbIn <- nrow(d)
message(paste("   ", nPrbIn, "probes in input file"))

# set up input model parameters
message("identifying segmentation type")
IS_SEGMENT_SAMPLE <- SEGMENTING_ACTION == "SEGMENT_SAMPLE"
if(!IS_SEGMENT_SAMPLE) {
    d[, ":="(
        INF_IN = "NA", # don't already know the presumed probe informativity
        AS_IN  = "NA"  # assume all probes might be informative
    )]
}
d[, ZYG := abs(BAF - 0.5) + 0.5]  # convert B allele frequency to zygosity
d[, INF_OUT := ifelse(is.na(ZYG) | ZYG >= extreme_zyg, 'no', 'yes')] # informativity based on array data
d[, ":="(
    LRR_RAW = LRR,
    BAF_RAW = BAF
)]

# setup for probe masking and discarding
# at this point, `refine` probes should all be informative, but `segment` probes may have missing LRR and BAF
# when refining a model, newly unusuable probes (e.g., with extreme GC) are discarded
# when segmenting an array, masked probes are retained in data.table to match other arrays but are not used for analysis
d[, masked := FALSE] 
discardProbes <- function(msg, toDiscard){
    message(msg)
    nBeforeDiscard <- nrow(d)
    nDiscarded <- sum(toDiscard)
    d <<- d[!toDiscard]
    nAfterDiscard <- nrow(d)
    message(paste("   ", nDiscarded, "of", nBeforeDiscard, "probes discarded;", nAfterDiscard, "probes remain"))
}
maskProbes <- function(msg, toMask){
    message(msg)
    nProbes <- nrow(d)
    nMasked <- sum(toMask)
    d[toMask, masked := TRUE]
    nUnmasked <- d[masked == FALSE, .N]
    message(paste("   ", nMasked, "of", nProbes, "probes masked;", nUnmasked, "probes remain unmasked"))
}
rejectProbes <- if(IS_SEGMENT_SAMPLE) maskProbes else discardProbes

# identify probes with usable data
# must be present in model AND have valid numerical LRR and BAF values
message("determining probe usability")
rejectProbes("require CN_IN",     is.na(d$CN_IN))
rejectProbes("require CN_IN > 0", d$CN_IN == 0) # microarrays are useless for mosaic SV detection when the modal state has no input DNA
rejectProbes("require LRR",       is.na(d$LRR))
rejectProbes("require BAF",       is.na(d$BAF))
if(!IS_SEGMENT_SAMPLE) rejectProbes("require FRAC_GC", is.na(d$FRAC_GC))

# sort probes; raw array is NOT necessarily sorted
message("sorting probes")
d <- d[order(CHROM, POS, PRB_NAME)] # yes, there can be duplicate SNP positions over >1 probe

#------------------------------------------------------------------
# apply/calculate model corrections
#------------------------------------------------------------------

# apply GC model corrections when segmenting
if(IS_SEGMENT_SAMPLE){

    # recover refinement metadata for GC bias correction
    message("applying GC bias corrections")
    gcData <- readRDS(GC_BIAS_FILE)
    gcFittableCNs <- gcData$gcFittableCNs
    gcBiasFits    <- gcData$gcBiasFits
    d <- merge(
        d,
        gcData$probe_perc_gc,
        by = 'PRB_NAME',
        all.x = TRUE,
        sort = FALSE
    )
    rejectProbes("require PERC_GC", is.na(d$PERC_GC))
    for(i in 1:length(gcFittableCNs)){
        j <- d[, CN_IN == gcFittableCNs[i]]
        offsets <- d[j, gcBiasFits[[i]]$offset[PERC_GC]]
        d[
            j, 
            ":="(
                LRR_CORR = LRR_RAW - offsets,
                LRR_MED  = LRR_MED - offsets
            )
        ]
    }

# calculate GC bias corrections when modeling
} else {
    message("calculating GC corrections")
    source(file.path(LIB_DIR, "calculate_GC_corr.R"))
}

# allow CN_OUT to range from 0 (homozygous deletion) to +1 gain relative to gcFittableCNs
hmmCNs <- 0:min(max(gcFittableCNs) + 1, maxAllowedCN)
maxHmmCN <- max(hmmCNs)

# only continue with probes from CN_IN states that are well modeled
# i.e., disallow changes INTO hmmCNs from OUTSIDE of that range
# CN_IN==0 previously dicarded above, even though it is always an allowed CN_OUT state
rejectProbes(paste("require CN_IN <=", maxHmmCN), d$CN_IN > maxHmmCN)

# use GC bias correction to calculate LRR_CORR for all copy number states
#   matrix rows are probes 
#   columns are allowed HMM CN states, starting with CN0
#   values are LRR_CORR assuming that probe state was truly that column's CN (true for only one of the columns)
# thus, the correction applied to LRR_RAW varies by the HMM output hidden state
# CN emission probabilites are calculated from lrrCorr relative to the genome-wide LRR_CORR model for each CN
message("modeling GC corrections")
lrrCorr <- sapply(hmmCNs, function(CN){
    # no GC correction applied to CN0
    if(CN == 0) return(d$LRR_RAW) 
    # GC correction for -1 loss or +1 gain relative to gcFittableCNs taken as the nearest actual correction
    cni <-      if (CN > max(gcFittableCNs)) which(gcFittableCNs == max(gcFittableCNs)) 
           else if (CN < min(gcFittableCNs)) which(gcFittableCNs == min(gcFittableCNs)) 
           else                              which(gcFittableCNs == CN)
    d[, LRR_RAW - gcBiasFits[[cni]]$offset[PERC_GC]]
})

# adjust LRR_CORR/lrrCorr and BAF to idealized values before sample modeling
if(IS_SEGMENT_SAMPLE){
    message("rescaling LRR and BAF to idealized values")
    source(file.path(LIB_DIR, "apply_model_corr.R"))
    lrr_lim <- c(-2.5, 2.5)
    zyg_lim <- c(0.5,  1.1)
}

# establish the chromsomes with non-discarded probes (no more probes will be discarded)
message("identifying chromosomes with data")
chroms <- d[, unique(CHROM)]
print(chroms)

#------------------------------------------------------------------
# determine state emissison probabilities
#------------------------------------------------------------------

# bin ZYG and LRR values as observation states, using integer bins
message("binning observation states")
n_zyg_bins   <- 150 # per unit, i.e. bins from 0 to 1, 100 gives 51 unique ZYG bins from 0.5 to 1
zyg_bins     <- round(0.5 * n_zyg_bins, 0) : round(max(d$ZYG, na.rm=TRUE) * n_zyg_bins, 0) # ranges of modeled/observed bins
d$zyg_bin    <- round(d$ZYG * n_zyg_bins, 0) # the actual bins for each probe observation
n_lrr_bins   <- 30 
lrr_bins     <- round(min(lrrCorr, na.rm=TRUE) * n_lrr_bins, 0) : round(max(lrrCorr, na.rm=TRUE) * n_lrr_bins, 0)
d$lrr_bin    <- round(d$LRR_CORR * n_lrr_bins, 0) # as used for state modeling based on input CN
lrrBins      <- sapply(hmmCNs, function(CN) as.character(round(lrrCorr[,CN + 1] * n_lrr_bins, 0))) # as used for HMM based on candidate output CN
bad_zyg_prob <- 1 / (n_zyg_bins * (1 - 0.5)) # bad probes have equal probability for any state
bad_lrr_prob <- 1 / (n_lrr_bins * (lrr_lim[2] - lrr_lim[1]))
bad_prb_prob <- bad_zyg_prob * bad_lrr_prob
zyg_inc      <- 1/2
lrr_inc      <- 1/2

# load pre-defined HMM states
message("load predefined output states")
ot <- fread(file.path(LIB_DIR, "state_table.txt"), header=TRUE)
ot <- ot[copy_number <= maxHmmCN]
if(!IS_SEGMENT_SAMPLE) ot <- ot[mosaic == 'no']
for(col in c('zyg_mean','zyg_stdev','lrr_mean','lrr_stdev')) ot[[col]] <- as.numeric(ot[[col]])

# probes used for model fitting based on data
message("fitting ZYG and LRR emission probablitities by CN_IN")
arrayFittable <- d[, masked == FALSE & !is.na(LRR_CORR) & LRR_CORR >= lrr_lim[1] & LRR_CORR <= lrr_lim[2]] 
message(paste("   ", sum(arrayFittable), "of", sum(d$masked == FALSE), "unmasked probes have LRR within range for fitting"))

# fit and plot statistics across array
# model_array.R fills mean/sd values into ot table and creates zyg_eps and lrr_eps distribution objects
# uses GC-adjusted LRR values, i.e., LRR_CORR, via lrr_bins and d$lrr_bin
cn_col <- 'CN_IN'
boxes  <- TRUE
source(file.path(LIB_DIR, "model_array.R"), local=TRUE, print.eval=TRUE)

# set keys for rapid lookup of binned ZYG and LRR emission probablities
# zyg_eps and lrr_eps have one row per binned state with the associated emission probability
for(key in names(lrr_eps)) {
    lrr_eps[[key]][, bin_ := as.character(bin)]
    setkey(lrr_eps[[key]], bin_)
}
for(key in names(zyg_eps)) {
    zyg_eps[[key]][, bin_ := as.character(bin)]
    setkey(zyg_eps[[key]], bin_)
}

#------------------------------------------------------------------
# build and solve the HMM
#------------------------------------------------------------------

# fill the table of informative probabilities
# i.e. the probability of finding an informative probe state in an allelic state
message("initializing HMM parameters")
pInf <- d[masked == FALSE & INF_OUT=='yes', .N] / d[masked == FALSE, .N] # fraction of informative probes over entire array
pINF <- list()
for (allelism in as.character(0:2)){
    pINF[[allelism]] <- list()
    for (INF in c('yes','no')){
        pINF[[allelism]][[INF]] <-
            if(allelism=='0'){
                0.5  
            } else if(allelism=='1'){
                if(INF=='yes') BAD_PROBE_FREQ else 1-BAD_PROBE_FREQ
            } else {
                if(INF=='yes') pInf           else 1-pInf
            }  
    }
}

# initialize HMM object
uot <- ot[zyg_mean != -88 & lrr_mean != -88] # usable output states that exist in model
asts    <- as.vector(unique(uot$allelic_state)) # the possible HMM output hidden states
astcns  <- sapply(asts, function(ast) uot[allelic_state==ast,min(copy_number)], simplify = FALSE, USE.NAMES = TRUE) # named lists
astasms <- sapply(asts, function(ast) uot[allelic_state==ast,min(allelism)],    simplify = FALSE, USE.NAMES = TRUE)
asmscms <- sapply(asts, function(ast) uot[allelic_state==ast,mosaic[1]],        simplify = FALSE, USE.NAMES = TRUE)

message("calculating probe allelic state emission probabilities")
frc_zyg_names <- names(frc_zyg)
d[, zyg_bin := as.character(zyg_bin)]
source(file.path(LIB_DIR, "simple_HMM.R"))
HMM <- HMM_init(1 - PERSISTENCE, sapply(1:length(asts), function(hs_i){
    AS_OUT       <- asts[hs_i]
    CN_OUT       <- astcns[[AS_OUT]]
    CN_OUT_      <- as.character(CN_OUT)
    allelism_OUT <- astasms[[AS_OUT]]
    mosaic_OUT   <- asmscms[[AS_OUT]]
    message(paste("  ", AS_OUT, " [", CN_OUT_, ":", allelism_OUT, "]", sep=""))
    frc_zyg_out <- CN_OUT_ %in% frc_zyg_names
    leps <- lrr_eps[[CN_OUT_]]
    pinf <- pINF[[as.character(allelism_OUT)]]

    # set probability adjustment based on probe informativity (may be overrident by sample data below)
    zepsKey <- if(allelism_OUT == 0){ # output state makes no prediction at all about zygosity
        "0"
    } else if(allelism_OUT == 1){ # output state only emits to non-informative observation states
        "1"
    } else if (frc_zyg_out){ # informative probes in biallelic states (should stay informative when segmenting)
        CN_OUT_   # the actual zygosity distribution, used for single-peak distributions
    } else {
        AS_OUT    # modeled zygosity distributions, used for multi-peak distributions and mosaics
    }

    # adjust likelihoods of output states based on known modeled states when segmenting
    if(IS_SEGMENT_SAMPLE){
        d[, ":="(
            adj_AS_OUT = ifelse(allelism_OUT > unlist(astasms[AS_IN]), 0, 1), # known monoallelic can never gain informativity
            adj_CN_OUT = ifelse(
                CN_IN == 0 & CN_OUT >= 1, 
                0, # DNA can never become "unmissing"
                ifelse(
                    CN_OUT == CN_IN,
                    ifelse(
                        allelism_OUT == unlist(astasms[AS_IN]), 
                        1,               # no SV state most favored, weight fully
                        PRESERVATION[2]  # downgrade CN neutral LOH; double-hit since both a loss and a gain
                    ),
                    ifelse(
                        mosaic_OUT == 'yes', 
                        PRESERVATION[1],                  # downgrade mosaic changes similar to CNC=1
                        PRESERVATION[abs(CN_OUT - CN_IN)] # downgrade CN changes in a graded fashion, i.e. more penalty at larger CNC
                    )
                )
            ),
            pinf_OUT = 1 # informativity is deterministic when segmenting
        )]
    } else {
        d[, ":="(
            adj_AS_OUT = 1, # we don't know about input pre-genotyping allelic states
            adj_CN_OUT = ifelse(
                CN_OUT == CN_IN,
                1,
                PRESERVATION[abs(CN_OUT - CN_IN)] # downgrade CN changes in a graded fashion, i.e. more penalty at larger CNC
            ),
            pinf_OUT = unlist(pinf[INF_OUT])
        )]
    }

    # return emission probablities reflecting the combination of:
    #   when segmenting, the modeled input state of the reference clone
    #   the sample-specific observed values
    #   the output state, including an output-state-specific GC bias correction
    d[, log(
        adj_AS_OUT * # probability of a probe being informative allows modeling monoallelic vs. biallelic states (e.g. LOH)
        adj_CN_OUT * 
        pinf_OUT * 
        leps[lrrBins[, CN_OUT + 1], prob] *
        ifelse(
            INF_IN == 'no' | # known uninformative probes only emit non-informative observation states (cannot become informative)
            (!IS_SEGMENT_SAMPLE & INF_OUT=='no'), # first discovering that a probe is uniformative
            zyg_eps[["1"]][    zyg_bin, prob],
            zyg_eps[[zepsKey]][zyg_bin, prob]
        )
    )]
}))

# run the HMM on each chromosome
message("solving HMM")
for (chrom in chroms){ 
    message(paste("  ", chrom))
    flt <- d$CHROM==chrom
    hs_i <- HMM_viterbi(HMM, flt)
    d[flt, ":="(
        hs = hs_i,
        AS_OUT = asts  [hs_i],
        CN_OUT = unlist(astcns[hs_i])
    )]
}
print(dcast(d, CHROM ~ CN_OUT, fun.aggregate = length, value.var = "CN_OUT"))
print(dcast(d, CN_IN ~ AS_OUT, fun.aggregate = length, value.var = "AS_OUT"))

# update informativity based on model output
# probes in monoallelic runs of homozygosity are never considered informative
# probes in a biallelic state have informativity determined by observed zygosity
# calculate additional changes in output relative to input modeled states
message("calculating derived output columns")
d[, ":="(
    INF_OUT = ifelse(
        masked == TRUE,
        NA,
        ifelse(
            unlist(astasms[AS_OUT]) <= 1,
            'no',
            INF_OUT
        )
    ),
    CNC = ifelse(masked == TRUE, NA, CN_OUT - CN_IN),
    LOH = if(IS_SEGMENT_SAMPLE) ifelse(masked == TRUE, NA, ifelse(INF_OUT != INF_IN, 'yes', 'no')) else NA
)]

# during refinement, discard probes with CNC of two or more
# these are often CN2_11 LRR outliers falsely called as CN4_22
# this behavior of a sample-set median marks untrustworthy probes that:
#   don't match typical probe performance
#   lead to excessive state breaks during sample HMM segmentation
if(!IS_SEGMENT_SAMPLE){
    refinementRejections <- d$CNC >= 2
    rejectProbes("require CNC <= 1", refinementRejections)
    arrayFittable <- arrayFittable[!refinementRejections] # for use by final model_array.R
    print(dcast(d, CN_IN ~ AS_OUT, fun.aggregate = length, value.var = "AS_OUT"))

    # save HMM for future use by compare
    message("saving HMM for future use")
    HMM$emissProbs <- HMM$emissProbs[!refinementRejections]
    # if(IS_SEGMENT_SAMPLE) HMM <- HMM_set_median_prob(HMM)
    save(HMM, file=HMM_FILE)
}

# fit and plot refined statistics across array
cn_col <- 'CN_OUT'
source(file.path(LIB_DIR, "model_array.R"))

# write the final output probes BED file
message("writing probes BED")
header <- paste(PROBE_COL, collapse="\\t")
header <- paste("#", header, sep="") # comment the header for bgzip
pipe <- paste("awk 'BEGIN{OFS=\"\\t\";print \"", header,"\"}{$2+=0;$3+=0;print $0}' | bgzip -c > ",
              PROBES_FILE, sep="")
write.table(
    d[, .SD, .SDcols = PROBE_COL],
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)
nPrbOut <- nrow(d)
message(paste("   ", nPrbOut, "probes written to disk"))

# write refine GC metadata for use by segment
if(IS_SEGMENT_SAMPLE){
    source(file.path(LIB_DIR, "find_SVs.R"))
} else {
    message("writing gc_bias RDS")
    saveRDS(
        list(
            probe_perc_gc = d[, .(PRB_NAME, PERC_GC)],
            gcFittableCNs = gcFittableCNs,
            gcBiasFits    = gcBiasFits
        ), 
        file = GC_BIAS_FILE
    )
}
message("segment_HMM.R done")

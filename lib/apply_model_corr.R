# adjust sample LRR_CORR (GC bias corrected) and BAF values according to the input model
# occurs in two steps:
# STEP 1 - rescale LRR_CORR and BAF values to their idealized scale
# STEP 2 - move sample LRR_CORR and BAF values toward their idealized input state
#          by the amount that the probe's median differs from the median of medians
# this resists SVs by accounting for probe-specific biases
# and makes values more sensible e.g. log2(1/2), although does result in ZYG > 1.0

# TODO: if linear regression is insufficient
# could use a non-linear regression to rescale LRR_CORR and BAF

# TODO: still concerned that moving the _CNV_ probes by the same amount
# as the non-CNV probes might move them too much toward target?
# scale the movement by the degree of distance from target???

# known deficiency: if the model is almost exclusively CN2 (i.e. normal female)
# cannot perform the linear regression meaningfully!

# LRR STEP 1: re-scale LRR_CORR values to log2 theoretical values
message("    LRR Step 1")
CN_INs      <- sort(unique(c(1:3, d$CN_IN))) # all CN values according to model, and always includes 1:3
CN_INs      <- CN_INs[!is.na(CN_INs) & CN_INs>0] # CN0 creates log2 error and not useful for fitting
lrrX_mom    <- c() # median of input LRR_CORR medians on its arbitrary power-X scale
lrrX_n      <- c() # No. of observations by CN, used as linear regression weights
lrrX_offset <- 0   # amount to shift the lrrX linear relationship to set CN2 to zero
for (i in 1:length(CN_INs)){
    CN <- CN_INs[i]
    dd <- d[CN_IN==CN, LRR_MED] # LRR_MED is GC bias-corrected by caller
    dd <- dd[!is.na(dd)]
    if(length(dd) > 0){
        mom <- median(dd) # the median of medians
        if(CN==2) lrrX_offset <- mom # always use CN2 as LRR reference
        lrrX_mom[i] <- mom 
        lrrX_n[i] <- length(dd)
    } else {
        lrrX_mom[i] <- 0 
        lrrX_n[i]   <- 0 
    }
}
lrrX_mom   <- lrrX_mom   - lrrX_offset    # apply the lrrX shift
d$LRR_MED  <- d$LRR_MED  - lrrX_offset
d$LRR_CORR <- d$LRR_CORR - lrrX_offset # LRR_CORR is GC bias-corrected
for(j in 1:ncol(lrrCorr)) lrrCorr[,j] <- lrrCorr[,j] - lrrX_offset
lrr2_tgt  <- log2(CN_INs/2) # theoretical lrr2 values by CN
lrr_fit   <- lm(lrr2 ~ lrrX,
                data=data.table(lrr2=lrr2_tgt, lrrX=lrrX_mom),
                weights=lrrX_n)
lrr2_mom   <- predict(lrr_fit) # apply the LRR rescaling
d$LRR_MED  <- predict(lrr_fit, data.table(lrrX=d$LRR_MED))
d$LRR_CORR <- predict(lrr_fit, data.table(lrrX=d$LRR_CORR))
for(j in 1:ncol(lrrCorr)) lrrCorr[,j] <- predict(lrr_fit, data.table(lrrX=lrrCorr[,j]))

# LRR STEP 2: calculate and apply lrr2 correction per probe
# CN0 does not get adjusted by this step
message("    LRR Step 2")
lrr2_tgts <- numeric() # used by model_LRR.R
for (i in 1:length(CN_INs)){
    CN  <- CN_INs[i]
    flt <- d$CN_IN==CN
    lrr2_corr <- lrr2_mom[i] - d[flt,LRR_MED]
    d[flt, LRR_CORR := LRR_CORR + lrr2_corr]
    for(j in 1:ncol(lrrCorr)) lrrCorr[flt,j] <- lrrCorr[flt,j] + lrr2_corr
    lrr2_tgts[CN] <- lrr2_mom[i]
}
rm(flt, lrr2_corr)

# BAF STEP 1: re-scale BAF values to theoretical values
message("    BAF Step 1")
tgt_baf  <- c( # low-side BAF target values for informative probes
    CN2_11=0.5,   # for alleleic states that have them
    CN3_12=0.33,
    CN4_13=0.25,
    CN4_22=0.5,
    CN5_14=0.2,
    CN5_23=0.4
)
tgt_baf_in <- ifelse( # returns NA for unusable probes
    d$INF_IN == 'no', # INF determined in order by (i) model allelic state (ii) model median probe zygosity
    0,                # returns the low-side (<=0.5) target BAF based on probe informativity
    tgt_baf[as.character(d$AS_IN)] # as.character prevents use of AS_IN numeric index (even though it is character) WAS PREVIOUSLY A FACTOR!
)
baf_tgts_prb <- as.character(ifelse(d$BAF_MED < 0.5, tgt_baf_in, 1-tgt_baf_in)) # flip to high-side as needed
baf_tgts     <- unique(baf_tgts_prb[!is.na(baf_tgts_prb)]) # target peaks found in model, as strings
baf_tgt_vals <- as.numeric(baf_tgts) # and their numeric equivalents

baf_mom <- numeric() # the median of medians, by BAF peak
baf_n   <- numeric() # No. of observations, used as linear regression weights
for(baf_tgt in baf_tgts){
    dd <- d[masked == FALSE & baf_tgts_prb==baf_tgt,BAF_MED]
    dd <- dd[!is.na(dd)]
    baf_mom[baf_tgt] <- median(dd)
    baf_n[baf_tgt]   <- length(dd)
}

baf_fit   <- lm(bafOut ~ bafIn,
                data=data.table(bafOut=baf_tgt_vals, bafIn=baf_mom),
                weights=baf_n)
baf_mom   <- predict(baf_fit) # apply the BAF rescaling
d$BAF_MED <- predict(baf_fit, data.table(bafIn=d$BAF_MED))
d$BAF     <- predict(baf_fit, data.table(bafIn=d$BAF))
names(baf_mom) <- baf_tgts

# BAF STEP 2: calculate and apply BAF correction per probe
message("    BAF Step 2")
for(baf_tgt in baf_tgts){
    flt <- d$masked==FALSE & baf_tgts_prb == baf_tgt
    baf_corr <- baf_mom[baf_tgt] - d[flt,BAF_MED]
    d[flt, BAF := BAF + baf_corr]
}
rm(tgt_baf_in, baf_tgts_prb, dd, flt, baf_corr)

# adjust additional values based on rescaling
d[, ":="(
    LRR = LRR_CORR, # for eventual plotting
    ZYG = abs(BAF - 0.5) + 0.5
)]
d[, INF_OUT := ifelse(is.na(ZYG) | ZYG >= extreme_zyg, 'no', 'yes')] # informativity based on array data

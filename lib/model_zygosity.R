# step through CN states to progressively collect
# the zygosity distributions present in array data,
# guided by modeled CN for the subject clone's baseline/reference genome

begSD <- 0.05
minSD <- 0.001
maxSD <- 0.1

# general function for fitting multiple gaussians to ZYG plots
fit_zygosity <- function(CN, type, formula, start, lower, upper){

    # set the data based on requested CN
    if(CN == 1) {
        dd <- d[arrayFittable & ZYG>=extreme_zyg]$zyg_bin
    } else {
        dd <- d[arrayFittable & ZYG <extreme_zyg & d[[cn_col]]==CN]$zyg_bin
    }
    dd <- dd[!is.na(dd)]
    
    # establish the zygosity histogram
    h <- aggregate(dd, by=list(dd), length)
    Z <- zyg_bins
    N <- sum(h[[2]])
    n <- as.numeric(sapply(Z, function(ZYG) h[h[[1]]==ZYG,2]))
    n[is.na(n)] <- 0
    F <- n / N

    # initialize the histogram plot
    label <- paste('CN', CN, 'ZYG', sep="_")
    jpgFile <- paste(PLOT_PREFIX, cn_col, label, type, 'jpg', sep=".")
    jpeg(jpgFile, 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    x <- Z/n_zyg_bins
    plot(x, F, type="h", xlim=zyg_lim, xlab="Zygosity", ylab="Frequency")
    lines(x, F, col="forestgreen", lwd=1)
    
    # fit the requested model (one or more gaussians) to the histogram
    start <- binPrm(start, n_zyg_bins)
    lower <- binPrm(lower, n_zyg_bins)
    upper <- binPrm(upper, n_zyg_bins)
    fit <- withWarningsErrors(nls(
        formula,
        data=data.table(Z=Z,F=F),
        algorithm="port",
        start=start,
        lower=lower,
        upper=upper,
        control=list(maxiter=200, warnOnly=TRUE)
    ))
    
    # analyze the fit for errors and warnings
    if(!is.null(fit$error)){
        # singular gradient matrix at initial parameter estimates
        graphics.off()
        unlink(jpgFile)
        return( list(error=TRUE) )
    } else if(!is.null(fit$warning)){
        # Convergence failure: singular convergence (7)
        # Convergence failure: false convergence (8)
        # Convergence failure: function evaluation limit reached without convergence (9)        
        # do nothing, convergence failure is handled later
    }
    fit <- fit[[1]]
    modeled <- predict(fit) # includes EPs for all superimposed gaussians
    
    # add the model lines to the plot
    # and collect modeled EPs for individual gaussians
    c <- coef(fit)
    nDist <- length(start) / 3
    zyg_eps <- list()
    frc_zyg <- list() 
    for(i in 1:nDist){
        L_ <- paste("L", i, sep="")
        M_ <- paste("M", i, sep="")
        S_ <- paste("S", i, sep="")
        lines(x, c[[L_]] * dnorm(Z, c[[M_]], c[[S_]]), col="blue", lwd=1)
        zyg_eps[[paste(CN, i, sep=":")]] <-
            data.table( bin=Z, prob=normalizeEPs(sapply(Z, function(bin) {
                diff( pnorm(bin+c(-zyg_inc, zyg_inc), c[[M_]], c[[S_]]) )
            }), bad_zyg_prob) )
            
    }
    lines(x, modeled, col="red", lwd=1)
  
    # collect the actual EP distribution
    CN <- as.character(CN)
    zyg_eps[[CN]] <- data.table(bin=Z, prob=normalizeEPs(F, bad_zyg_prob))
    if(nDist==1) frc_zyg[[CN]] <- 1 # TODO: prefer actual if only one actual peak?? (regardless of no. of modeled peaks)

    # finish and return
    graphics.off()
    dat <- as.list( binPrm(c, n_zyg_bins, unbin=TRUE) )
    dat$deviance <- deviance(fit)       # boolean
    dat$isConv   <- fit$convInfo$isConv # double
    dat$error    <- FALSE
    dat$jpgFile  <- jpgFile
    dat$zyg_eps  <- zyg_eps
    dat$frc_zyg  <- frc_zyg
    dat
}

# set ZYG parameters for probes with all A or B alleles
# includes CN1 if available, but not required
message("fitting non-informative zygosity")
z1 <- fit_zygosity(
    1, # will use all CN states here, not just CN=1
    'simple',
    F ~ L1 * dnorm(Z, M1, S1),
    start=c(L1=1, M1=1.0,  S1=begSD), # keep the edge distrubutions close to symmetric about limits
    lower=c(L1=1, M1=0.95, S1=minSD),
    upper=c(L1=1, M1=1.05, S1=maxSD)
)
zyg_eps <- c(zyg_eps, z1$zyg_eps)
frc_zyg <- c(frc_zyg, z1$frc_zyg)
ot_flt <- ot$informative=='no'
ot[ot_flt, zyg_mean  := z1$M1]
ot[ot_flt, zyg_stdev := z1$S1]
# CN1 and uninformative always use actual distribution, forced during HMM assembly

# set ZYG parameters for CN2 and CN3 is coordinated manner
# NB: input models are expected to have CN=2 with sufficient probes!
# always ensure there are values set for CN=1,2,3
message("fitting CN2 zygosity")
inf_flt <- d$ZYG<extreme_zyg
z2smp <- fit_zygosity(
    2,
    'simple',
    F ~ L1 * dnorm(Z, M1, S1),
    start=c(L1=1, M1=0.50, S1=begSD),
    lower=c(L1=1, M1=0.49, S1=minSD),
    upper=c(L1=1, M1=0.54, S1=maxSD)
)
isSuff3 <- isSuff(3, inf_flt)
if(isSuff3){
    z2cmp <- fit_zygosity(
        2,
        'complex',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2),
        start=c(L1=0.75, M1=0.50, S1=begSD,
                L2=0.25, M2=0.67, S2=begSD),
        lower=c(L1=0,    M1=0.49, S1=minSD,
                L2=0,    M2=0.62, S2=minSD),
        upper=c(L1=1,    M1=0.54, S1=maxSD,
                L2=1,    M2=0.72, S2=maxSD)
    )
    z2 <- compareFit(z2smp, z2cmp) 
    message("fitting CN3 zygosity")
    z3smp <- fit_zygosity(
        3,
        'simple',
        F ~ L1 * dnorm(Z, M1, S1),
        start=c(L1=1, M1=0.67,  S1=begSD),
        lower=c(L1=1, M1=0.62,  S1=minSD),
        upper=c(L1=1, M1=0.72,  S1=maxSD)
    )
    z3cmp <- fit_zygosity(
        3,
        'complex',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2) +
            L3 * dnorm(Z, M3, S3),
        start=c(L1=0.20, M1=z2$M1, S1=z2$S1,
                L2=0.60, M2=0.67,  S2=begSD,
                L3=0.20, M3=0.75,  S3=begSD),
        lower=c(L1=0,    M1=z2$M1, S1=z2$S1,
                L2=0,    M2=0.62,  S2=minSD,
                L3=0,    M3=0.72,  S3=minSD),
        upper=c(L1=1,    M1=z2$M1, S1=z2$S1,
                L2=1,    M2=0.72,  S2=maxSD,
                L3=1,    M3=0.79,  S3=maxSD)
    )
    z3 <- compareFit(z3smp, z3cmp, list(M2='M1', S2='S1'))
    zyg_eps[["3:2"]] <- if(z3$type=='complex') { zyg_eps[["3:2"]] } else { z3smp$zyg_eps[["3:1"]] }
     
} else {
    # for simple well-behaved arrays, fit CN2 with single gaussian
    z2 <- z2smp
    zyg_eps <- c(zyg_eps, z2$zyg_eps)
    frc_zyg <- c(frc_zyg, z2$frc_zyg)
    message(paste("  ", 'simple', ifelse(z2$isConv, 'convergent', 'NOT convergent')))
    # always allow for duplication even if trained CN lacks CN=3
    message("projecting CN3 zygosity")
    tgt <- if(exists('baf_mom')) {
        ifelse(is.na(baf_mom['CN3_12']), 0.65, baf_mom['CN3_12'])
    } else {
        0.65 # empirically determined
    }
    z3 <- list(L2=1, M2=tgt, S2=z2$S1)
    modeled <- sapply(zyg_bins, function(bin) {
        diff( pnorm(bin+c(-zyg_inc, zyg_inc), z3$M2*n_zyg_bins, z3$S2*n_zyg_bins) )
    })
    zyg_eps[["3"]] <- data.table(bin=zyg_bins, prob=normalizeEPs(modeled, bad_zyg_prob))
    zyg_eps[["3:2"]] <- zyg_eps[["3"]]
    plotModeled(3, 'ZYG', zyg_bins, n_zyg_bins, modeled, "Zygosity", zyg_lim)   
}
PS <- "CN2_11"
ot_flt <- ot$probe_state==PS
ot[ot_flt, zyg_mean  := z2$M1]
ot[ot_flt, zyg_stdev := z2$S1]
zyg_eps[[PS]] <- zyg_eps[["2:1"]] # probe and allelic states synonymous when informative
PS <- "CN3_12"
ot_flt <- ot$probe_state==PS
ot[ot_flt, zyg_mean  := z3$M2]
ot[ot_flt, zyg_stdev := z3$S2]
zyg_eps[[PS]] <- zyg_eps[["3:2"]]

# set ZYG parameters for CN4 if they exist in input model OR if they can be extrapolated from CN3
if(isSuff(4, inf_flt)){
    message("fitting CN4 zygosity")
    z4smp <- fit_zygosity(
        4, 'simple',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2),
        start=c(L1=0.50, M1=z2$M1, S1=z2$S1,
                L2=0.50, M2=0.75,  S2=begSD),
        lower=c(L1=0,    M1=0.49,  S1=minSD,
                L2=0,    M2=0.71,  S2=minSD),
        upper=c(L1=1,    M1=0.54,  S1=maxSD,
                L2=1,    M2=0.79,  S2=maxSD)
    )
    z4cmp <- fit_zygosity(
        4, 'complex',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2) +
            L3 * dnorm(Z, M3, S3),
        start=c(L1=0.40, M1=z2$M1, S1=z2$S1,
                L2=0.20, M2=z3$M2, S2=z3$S2,
                L3=0.40, M3=0.75,  S3=begSD),
        lower=c(L1=0,    M1=0.49,  S1=minSD,
                L2=0,    M2=z3$M2, S2=z3$S2,
                L3=0,    M3=0.71,  S3=minSD),
        upper=c(L1=1,    M1=0.54,  S1=maxSD,
                L2=1,    M2=z3$M2, S2=z3$S2,
                L3=1,    M3=0.79,  S3=maxSD)
    )
    z4 <- compareFit(z4smp, z4cmp, list(M3='M2', S3='S2')) 
    ot_flt <- ot$informative=='no' & ot$copy_number==4
    ot[ot_flt, zyg_mean  := z1$M1]
    ot[ot_flt, zyg_stdev := z1$S1]
    PS <- "CN4_22"
    ot_flt <- ot$probe_state==PS
    ot[ot_flt, zyg_mean  := z4$M1]
    ot[ot_flt, zyg_stdev := z4$S1]
    zyg_eps[[PS]] <- zyg_eps[["4:1"]]
    PS <- "CN4_13"
    ot_flt <- ot$probe_state==PS  
    ot[ot_flt, zyg_mean  := z4$M3]
    ot[ot_flt, zyg_stdev := z4$S3]
    zyg_eps[[PS]] <- if(z4$type=='complex') { zyg_eps[["4:3"]] } else { z4smp$zyg_eps[["4:2"]] }
} else if(isSuff3){
    z3 <- z3smp
    zyg_eps <- c(zyg_eps, z3$zyg_eps)
    frc_zyg <- c(frc_zyg, z3$frc_zyg)
    message(paste("  ", 'simple', ifelse(z3$isConv, 'convergent', 'NOT convergent')))
    message("projecting CN4 zygosity")
    tgt <- if(exists('baf_mom')) {
        ifelse(is.na(baf_mom['CN4_13']), 0.72, baf_mom['CN4_13'])
    } else {
        0.72 # empirically determined
    }
    z4 <- list(L1=1, M1=tgt, S1=z3$S1, isProjected = TRUE)
    modeled <- sapply(zyg_bins, function(bin) {
        diff( pnorm(bin+c(-zyg_inc, zyg_inc), z4$M1*n_zyg_bins, z4$S1*n_zyg_bins) )
    })
    zyg_eps[["4"]] <- data.table(bin=zyg_bins, prob=normalizeEPs(modeled, bad_zyg_prob))
    zyg_eps[["4:1"]] <- zyg_eps[["4"]]
    plotModeled(4, 'ZYG', zyg_bins, n_zyg_bins, modeled, "Zygosity", zyg_lim) 
    PS <- "CN4_13"
    ot_flt <- ot$probe_state==PS  
    ot[ot_flt, zyg_mean  := z4$M1]
    ot[ot_flt, zyg_stdev := z4$S1]
    zyg_eps[[PS]] <- zyg_eps[["4:1"]]
    z2 <- z2smp
    PS <- "CN4_22"
    ot_flt <- ot$probe_state==PS
    ot[ot_flt, zyg_mean  := z2$M1]
    ot[ot_flt, zyg_stdev := z2$S1]
    zyg_eps[[PS]] <- zyg_eps[["CN2_11"]]
}

# set ZYG parameters for CN5, but only if it exists in input model
if(isSuff(5, inf_flt)){
    message("fitting CN5 zygosity")
    z5smp <- fit_zygosity(
        5, 'simple',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2),
        start=c(L1=0.50, M1=0.62,  S1=z4$S3,
                L2=0.50, M2=0.79,  S2=z4$S3),
        lower=c(L1=0,    M1=0.59,  S1=z4$S3*0.5,
                L2=0,    M2=0.75,  S2=z4$S3*0.5),
        upper=c(L1=1,    M1=0.65,  S1=z4$S3*1.5,
                L2=1,    M2=0.82,  S2=z4$S3*1.5)   
    )
    z5cmp <- fit_zygosity(
        5, 'complex',
        F ~ L1 * dnorm(Z, M1, S1) +
            L2 * dnorm(Z, M2, S2) +
            L3 * dnorm(Z, M3, S3) +
            L4 * dnorm(Z, M4, S4),
        start=c(L1=0.10, M1=z4$M1, S1=z4$S1,
                L2=0.40, M2=0.62,  S2=z4$S3,
                L3=0.10, M3=z4$M3, S3=z4$S3,
                L4=0.40, M4=0.79,  S4=z4$S3),
        lower=c(L1=0,    M1=z4$M1, S1=z4$S1,
                L2=0,    M2=0.59,  S2=z4$S3*0.5,
                L3=0,    M3=z4$M3, S3=z4$S3,
                L4=0,    M4=0.75,  S4=z4$S3*0.5),
        upper=c(L1=1,    M1=z4$M1, S1=z4$S1,
                L2=1,    M2=0.65,  S2=z4$S3*1.5,
                L3=1,    M3=z4$M3, S3=z4$S3,
                L4=1,    M4=0.82,  S4=z4$S3*1.5)   
    )
    z5 <- compareFit(z5smp, z5cmp, list(M4='M2', S4='S2',
                                        M2='M1', S2='S1'))
    ot_flt <- ot$informative=='no' & ot$copy_number==5
    ot[ot_flt, zyg_mean  := z1$M1]
    ot[ot_flt, zyg_stdev := z1$S1]
    PS <- "CN5_23"
    ot_flt <- ot$probe_state==PS  
    ot[ot_flt, zyg_mean  := z5$M2]
    ot[ot_flt, zyg_stdev := z5$S2]
    zyg_eps[[PS]] <- if(z5$type=='complex') { zyg_eps[["5:2"]] } else { z5smp$zyg_eps[["5:1"]] }
    PS <- "CN5_14"
    ot_flt <- ot$probe_state==PS  
    ot[ot_flt, zyg_mean  := z5$M4]
    ot[ot_flt, zyg_stdev := z5$S4]
    zyg_eps[[PS]] <- if(z5$type=='complex') { zyg_eps[["5:4"]] } else { z5smp$zyg_eps[["5:2"]] }
}

# set ZYG parameters for CN0; equally weight all outcomes since meaningless
message("projecting CN0 zygosity")
zyg_eps[["0"]] <- data.table(
    bin=zyg_bins,
    prob=1/length(zyg_bins)
)
ot_flt <- ot$copy_number==0
ot[ot_flt, zyg_mean  := 0.5]
ot[ot_flt, zyg_stdev := 0]

# when segmenting, create two mosaic states for CN1-2 and CN2-3
# note that here, low is always CN2, i.e. ZYG=0.5
message("projecting mosaic zygosities")
if(IS_SEGMENT_SAMPLE){
    zyg_eps[["1.5"]] <- data.table(
        bin=zyg_bins,
        prob=project_mosaic(zyg_eps[['2']]$prob, zyg_eps[['1']]$prob, bad_zyg_prob) # will be much wider...
    )
    PS <- "CN1.5_11"
    ot_flt <- ot$probe_state==PS
    ot[ot_flt, zyg_mean  := 0.75]
    ot[ot_flt, zyg_stdev := 0]
    zyg_eps[[PS]] <- zyg_eps[["1.5"]] 
    
    zyg_eps[["2.5"]] <- data.table(
        bin=zyg_bins,
        prob=project_mosaic(zyg_eps[['2']]$prob, zyg_eps[['3']]$prob, bad_zyg_prob) # ...than this
    )    
    PS <- "CN2.5_12"
    ot_flt <- ot$probe_state==PS
    ot[ot_flt, zyg_mean  := 0.585]
    ot[ot_flt, zyg_stdev := 0]
    zyg_eps[[PS]] <- zyg_eps[["2.5"]]
}

# clean  up
rm(inf_flt)

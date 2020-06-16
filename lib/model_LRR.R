
# step through CN states to progressively collect
# the LRR distributions present in array data
# guided by trained CN AND refined ZYG distributions

n_lrr_sd <- 1.5

# general function for fitting multiple gaussians to ZYG plots
fit_LRR <- function(CN, zyg_filter){

    # set the data based on requested CN
    dd <- d[dlrr & d[[cn_col]]==CN & zyg_filter, 'lrr_bin'] # dlrr subsumes usability filter
    dd <- dd[!is.na(dd)]
    
    # establish the LRR histogram and bins
    h <- aggregate(dd, by=list(dd), length)
    L <- lrr_bins
    N <- sum(h[[2]])
    n <- as.numeric(sapply(L, function(LRR) h[h[[1]]==LRR,2]))
    n[is.na(n)] <- 0     
    F <- n / N

    # initialize the histogram plot
    label <- paste('CN', CN, 'LRR', sep="_")    
    jpeg(paste(PLOT_PREFIX, cn_col, label, 'jpg', sep="."),
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    x <- L/n_lrr_bins
    plot(x, F, type="h", xlim=lrr_lim, xlab="LRR", ylab="Frequency")
    lines(x, F, col="forestgreen", lwd=1)
    
    # fit the requested model to the histogram
    fit <- nls(
        F ~ L1 * dnorm(L, M1, S1),
        data=data.frame(L=L,F=F),
        algorithm="port",        
        start=c(L1=1, M1=weighted.mean(L, F), S1=n_lrr_bins * 0.1),
        lower=c(L1=1, M1=min(L),              S1=n_lrr_bins * 0.001),
        upper=c(L1=1, M1=max(L),              S1=n_lrr_bins * 1)
    )
    modeled <- predict(fit)  
    
    # add the model lines to the plot
    c <- coef(fit)    
    lines(x, modeled, col="red", lwd=1)
    
    # LRR HMM always uses actual, not modeled, distributions
    # and are always fitting just a single distribution
    lrr_eps[[as.character(CN)]] <<- data.frame(bin=L, prob=normalizeEPs(F, bad_lrr_prob))
    
    # finish and return
    graphics.off()
    as.list( binPrm(c, n_lrr_bins, unbin=TRUE) )
}

# set LRR parameters for CN2
# use mean +/- n_lrr_sd stdev ZYG to find LRR
# NB: input models are expected to have CN=2 with sufficient probes!
message("fitting CN2 LRR")
inf_flt <- d$ZYG >= 0.5 &
           d$ZYG <= z2$M1 + n_lrr_sd * z2$S1
l2 <- fit_LRR(2, inf_flt)
ot[ot$copy_number==2,'lrr_mean']  <- l2$M1
ot[ot$copy_number==2,'lrr_stdev'] <- l2$S1

# set LRR parameters for CN1
# always ensure there are values set for CN-0,1,2,3
inf_flt <- d$ZYG>=extreme_zyg
if(isSuff(1, inf_flt)){
    message("fitting CN1 LRR")
    l1 <- fit_LRR(1, inf_flt)
} else {
    message("projecting CN1 LRR")   
    tgt <- if(exists('lrr2_tgts')) {
        lrr2_tgts[1]
    } else {
        -0.33 # default LRR mean empirically determined
    }
    l1 <- list(M1=tgt, S1=l2$S1)
    modeled <- sapply(lrr_bins, function(bin) {
        diff( pnorm(bin+c(-lrr_inc, lrr_inc), l1$M1*n_lrr_bins, l1$S1*n_lrr_bins) )
    })
    lrr_eps[["1"]] <- data.frame(bin=lrr_bins, prob=normalizeEPs(modeled, bad_lrr_prob))
    plotModeled(1, 'LRR', lrr_bins, n_lrr_bins, modeled, "LRR", lrr_lim)  
}
ot[ot$copy_number==1,'lrr_mean']  <- l1$M1 # save all critical values
ot[ot$copy_number==1,'lrr_stdev'] <- l1$S1

# set LRR parameters for CN3
inf_flt <- d$ZYG >= z3$M2 - n_lrr_sd * z3$S2 &
           d$ZYG <= z3$M2 + n_lrr_sd * z3$S2
if(isSuff(3, inf_flt)){
    message("fitting CN3 LRR")
    l3 <- fit_LRR(3, inf_flt)
} else {
    message("projecting CN3 LRR")
    tgt <- if(exists('lrr2_tgts')) {
        lrr2_tgts[3]
    } else {
        l2$M1 + (l2$M1-l1$M1)*0.5/0.67
    }
    l3 <- list(M1=tgt, S1=l2$S1)
    modeled <- sapply(lrr_bins, function(bin) {
        diff( pnorm(bin+c(-lrr_inc, lrr_inc), l3$M1*n_lrr_bins, l3$S1*n_lrr_bins) )
    })
    lrr_eps[["3"]] <- data.frame(bin=lrr_bins, prob=normalizeEPs(modeled, bad_lrr_prob))
    plotModeled(3, 'LRR', lrr_bins, n_lrr_bins, modeled, "LRR", lrr_lim) 
}
ot[ot$copy_number==3,'lrr_mean']  <- l3$M1
ot[ot$copy_number==3,'lrr_stdev'] <- l3$S1
    
# set LRR parameters for CN4
# only set values for CN 4 and 5 if they exist in input model
if(exists('z4')){
    inf_flt <- (d$ZYG >= 0.5 &
                d$ZYG <= z4$M1 + n_lrr_sd * z4$S1) |
               (d$ZYG >= z4$M3 - n_lrr_sd * z4$S3 &
                d$ZYG <= z4$M3 + n_lrr_sd * z4$S3)
    if(isSuff(4, inf_flt)){
        message("fitting CN4 LRR")
        l4 <- fit_LRR(4, inf_flt)
        ot[ot$copy_number==4,'lrr_mean']  <- l4$M1
        ot[ot$copy_number==4,'lrr_stdev'] <- l4$S1
    }
}

# set LRR parameters for CN5
if(exists('z5')){
    inf_flt <- d$ZYG >= z5$M2 - n_lrr_sd * z5$S2 &
               d$ZYG <= z5$M2 + n_lrr_sd * z5$S2 # don't use 0.8 ZYG, not well represented in array
    if(isSuff(5, inf_flt)){
        message("fitting CN5 LRR")
        l5 <- fit_LRR(5, inf_flt)
        ot[ot$copy_number==5,'lrr_mean']  <- l5$M1
        ot[ot$copy_number==5,'lrr_stdev'] <- l5$S1
    }        
}

# set LRR parameters for CN0 based on limit value 2SD below CN1 LRR mean
# any LRR above this limit takes the bad probe frequency for CN0
# all LRR bins below this limit are unpredictable and equally weighted
message("projecting CN0 LRR")
minL1 <- l1$M1 - 2*l1$S1
lrr_bin_vals <- lrr_bins / n_lrr_bins
nBin0 <- length(lrr_bin_vals[lrr_bin_vals < minL1])
nullP <- ifelse(lrr_bin_vals < minL1, 1/nBin0, 0)
lrr_eps[["0"]] <- data.frame(
    bin=lrr_bins,
    prob=normalizeEPs(nullP, bad_lrr_prob)
)
ot[ot$copy_number==0,'lrr_mean']  <- minL1
ot[ot$copy_number==0,'lrr_stdev'] <- 0

# when segmenting, create two mosaic states for CN1-2 and CN2-3
message("projecting mosaic LRR")
if(!modeling){
    lrr_eps[["1.5"]] <- data.frame(
        bin=lrr_bins,
        prob=project_mosaic(lrr_eps[['1']]$prob, lrr_eps[['2']]$prob, bad_lrr_prob)
    )
    ot[ot$copy_number==1.5,'lrr_mean']  <- -0.5
    ot[ot$copy_number==1.5,'lrr_stdev'] <- 0
    lrr_eps[["2.5"]] <- data.frame(
        bin=lrr_bins,
        prob=project_mosaic(lrr_eps[['2']]$prob, lrr_eps[['3']]$prob, bad_lrr_prob)
    )
    ot[ot$copy_number==2.5,'lrr_mean']  <- 0.333
    ot[ot$copy_number==2.5,'lrr_stdev'] <- 0
}


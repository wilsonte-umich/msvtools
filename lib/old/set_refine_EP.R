
# calculate all possible 'refine' emission probabilities using fit Gaussian distributions
for(i in 1:length(asts)){ # one probability for all output CN states for all input combinations
    allelic_state <- asts[i]
    allelism      <- astasms[i]
    message(paste("   ", allelic_state, " [", allelism, "]", sep=""))
    otrs <- ot[ot$allelic_state==allelic_state,]
    ep[[allelic_state]] <- list()
    for(i in 1:nrow(b)){
        cmb <- paste(b[i,'zyg_bin'], b[i,'lrr_bin'], sep=",")
        ep[[allelic_state]][[cmb]] <- list()
        for (INF in c('yes','no')){              # informativity is determined by model sometimes, not by observed zygosity   
            otr <- otrs[otrs$informative==INF,]  # choose ep based on probe informativity
            ep[[allelic_state]][[cmb]][[INF]] <- suppressWarnings(
                max(BAD_PROBE_FREQ,           # ensure non-zero probabilities in all states

                    pINF[[allelism]][[INF]] * # use probability of being informative to allow monoallelic vs. biallelic states (e.g. ROH)

                    otr$prob_scalar *         # edge ZYG distributions are twice the Gaussian probabilities since symmetric
                        diff( pnorm(b[i,'zyg_bin']+c(zyg_inc, -zyg_inc), otr$zyg_mean, otr$zyg_stdev) ) *
                        diff( pnorm(b[i,'lrr_bin']+c(lrr_inc, -lrr_inc), otr$lrr_mean, otr$lrr_stdev) ),
                    na.rm=TRUE)
            )                  
        }
    }
}


# set the 'refine' HMM emission probabilities for individual probes
# probe informativity (INF) determined here by zygosity (ZYG) of median probes
setEP <- function(AS, CN, RC, zyg_bin, lrr_bin, INF){ 
    cmb <- paste(zyg_bin, lrr_bin, sep=",")
    ep[[AS]][[cmb]][[INF]] *       # the baseline probability based on array aggregates
    if(isBam){
        rc_ep[[CN]][[round(RC,0)]] # limited by NGS read counts
    } else {
        1                          # unless they weren't used to train the model
    }
}

# individiually process each target allelic state, e.g. CN1_1, CN3_12, etc.
for(i in 1:length(asts)){   
    message(paste("  ", asts[i]))
    # individiually process each probe
    d[[asts[i]]] <- mapply(setEP, asts[i], astcns[i], d$RC, d$zyg_bin, d$lrr_bin, d$INF)
}


# set the 'segment' HMM emission probabilities for individual probes


# probe informativity (INF) determined here by refined input state model for cell line


setEP <- function(AS_OUT, AS_IN, INF_IN, zyg_bin, lrr_bin, INF_OUT){

    # SV can never gain informativity
    if(astasms[[AS_OUT]] > astasms[[AS_IN]]){
        IMPOSSIBLE
        
    # probe should not gain informativity
    } else if(INF_IN  == 'no' &
              INF_OUT == 'yes'){
        BAD_PROBE_FREQ 

    # probe should not lose informativity if allelic state is biallelic
    } else if(astasms[[AS_OUT]] == 2 &
              INF_IN  == 'yes' &
              INF_OUT == 'no'){
        BAD_PROBE_FREQ 

    } else {
        cmb <- paste(zyg_bin, lrr_bin, sep=",")
        
        
    }
    

    
    ep[[AS]][[cmb]][[INF]] *       # the baseline probability based on array aggregates


}

# individiually process each target allelic state, e.g. CN1_1, CN3_12, etc.
for(i in 1:length(asts)){   
    message(paste("  ", asts[i]))
    # individiually process each probe
    d[[asts[i]]] <- mapply(setEP, asts[i], d$AS_IN, d$INF_IN, d$zyg_bin, d$lrr_bin, d$INF_OUT)
}
 
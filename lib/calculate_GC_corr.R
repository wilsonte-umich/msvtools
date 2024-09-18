# perform GC bias correction of probe LRR values, only called during genome model refinement
# each fittable copy number state generates its own set of bias offsets from its modeled median LRR
# offsets are used downstream to correct LRR values before HMM
# GC bias correction became important as array density increased and started to
# have shared/overlapping fragments/GC spans over multiple adjancent probe positions

# determine which input copy number states have sufficient probes for GC bias correction
message("finding fittable copy numbers")
minFitProbes <- 25000
nProbesByCN <- d[, .(nProbes = .N), keyby = .(CN_IN)]
nProbesByCN[, gcFittable := CN_IN > 0 & nProbes >= minFitProbes]
print(nProbesByCN)
gcFittableCNs <- nProbesByCN[gcFittable == TRUE, CN_IN]
d[, gcFittable := CN_IN %in% gcFittableCNs]
message(paste("   ", sum(d$gcFittable), "of", length(d$gcFittable), "probes have a GC-fittable copy number"))

# exclude probes with extreme percentGC
message("finding fittable percent GC")
d[, PERC_GC := as.integer(FRAC_GC * 100)]
gc_lim <- c(0, 0)
gc <- 40
while(TRUE){
    if(d[PERC_GC == gc, .N < minFitProbes]) break
    gc <- gc - 1L
}
gc_lim[1] <- gc
gc <- 40
while(TRUE){
    if(d[PERC_GC == gc, .N < minFitProbes]) break
    gc <- gc + 1L
}
gc_lim[2] <- gc
rejectProbes("discarding extreme GC values", d[, !(PERC_GC %between% gc_lim)])

# fit GC bias by copy number
message(paste("   ", sum(d$gcFittable), "of", length(d$gcFittable), "remaining probes have a GC-fittable copy number"))
lrr_col <- "LRR_RAW"
source(paste(LIB_DIR, "plot_GC_bias.R", sep="/"))
for(i in 1:length(gcFittableCNs)){
    d[
        CN_IN == gcFittableCNs[i], 
        LRR_CORR := LRR_RAW - gcBiasFits[[i]]$offset[PERC_GC]
    ]
}
lrr_col <- "LRR_CORR"
source(paste(LIB_DIR, "plot_GC_bias.R", sep="/"))

# initialize
message("plotting GC bias correlation")

dd <- d[gcFittable == TRUE]
dd[, gc2 := FRAC_GC ^ 2]
gc <- gc_lim[1]:gc_lim[2] / 100
gc2 <- gc ^ 2

nGcPoints <- 25000
sd <- do.call(rbind, lapply(gcFittableCNs, function(CN){
    ddd <- dd[CN_IN == CN]
    if(nrow(ddd) > nGcPoints) ddd[sample(1:nrow(ddd), nGcPoints)] else ddd
}))
sd <- sd[sample(1:nrow(sd))]

# plot the LRR correlation
jpeg(paste(PLOT_PREFIX, 'GC_BIAS_LRR', lrr_col, 'jpg', sep="."),
       width=width*2, height=height*2, unit=units, pointsize=pointsize, res=resol)
plot(jitter(sd$FRAC_GC, a = 0.005), sd[[lrr_col]], pch=".", col=cnColorsAlpha[sd$CN_IN + 1],
     xlab='Fraction GC', ylab='Log R Ratio (LRR)',
     xlim=gc_lim / 100,  ylim=lrr_lim)

# add trendlines for CN 1 to 3
message("plotting trendlines")
gcBiasFits <- lapply(1:length(gcFittableCNs), function(i){
    ddd <- dd[CN_IN == gcFittableCNs[i]]
    ddd[, lrrFit := ddd[[lrr_col]]]
    lrr <- predict(
        lm(lrrFit ~ FRAC_GC + gc2, data = ddd), 
        newdata = data.frame(FRAC_GC = gc, gc2 = gc2)
    )
    lines(gc, lrr, col = cnColors[gcFittableCNs[i] + 1])
    if(lrr_col == "LRR_RAW"){
        rbind(
            data.frame( # for easier lookup downstream
                PERC_GC = 1:(gc_lim[1] - 1L),
                offset = NA
            ),
            data.frame(
                PERC_GC = gc_lim[1]:gc_lim[2],
                offset = lrr - ddd[, median(LRR_RAW, na.rm = TRUE)]
            )
        )
    } else {
        gcBiasFits[[i]]
    }
})

# finish up
graphics.off()

# plot the ZYG correlation, only once for visualization, since zygosity does not show a GC bias (as expected)
if(lrr_col == "LRR_RAW"){
    sd <- do.call(rbind, lapply(2:max(gcFittableCNs), function(CN){
        ddd <- dd[INF_OUT == "yes" & CN_IN == CN]
        if(nrow(ddd) > nGcPoints) ddd[sample(1:nrow(ddd), nGcPoints)] else ddd
    }))
    sd <- sd[sample(1:nrow(sd))]
    jpeg(paste(PLOT_PREFIX, 'GC_BIAS_ZYG', lrr_col, 'jpg', sep="."),
        width=width*2, height=height*2, unit=units, pointsize=pointsize, res=resol)
    plot(jitter(sd$FRAC_GC, a = 0.005), sd$ZYG, pch=".", col=cnColorsAlpha[sd$CN_IN + 1],
        xlab='Fraction GC', ylab='Zygosity',
        xlim=gc_lim / 100,  ylim=zyg_lim)
    graphics.off()
}

# memory management
rm(dd, sd)

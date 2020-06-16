
# plot BAF-LRR correlation, colored by trained or refined CN

# initialize
message("plotting ZYG-LRR correlation")
jpeg(paste(PLOT_PREFIX, cn_col, 'CORRELATION', 'jpg', sep="."), 
       width=width*2, height=height*2, unit=units, pointsize=pointsize, res=resol)

# downsample the data for the plot
# include CN=0, but exclude mosaics
dd <- d[!is.na(d[[cn_col]]) & d[[cn_col]]%%1==0,]
nPoints <- 100000
sd <- if(nrow(dd) > nPoints){
    dd[sample(1:nrow(dd), nPoints),]
} else {
    dd
}

# plot the correlation
plot(sd$ZYG, sd$LRR, pch=".", col=sd[[cn_col]]+1,
     xlab='Zygosity (ZYG)', ylab='Log R Ratio (LRR)',
     xlim=zyg_lim,          ylim=lrr_lim)

# add CN region indicators
if(boxes){
    nsd <- 2
    bx <- unique(ot[ot$zyg_mean>0 & ot$mosaic=="no",c('zyg_mean','zyg_stdev','lrr_mean','lrr_stdev')])
    for(i in 1:nrow(bx)){
        zmin <- bx[i,'zyg_mean'] - nsd * bx[i,'zyg_stdev']
        zmax <- bx[i,'zyg_mean'] + nsd * bx[i,'zyg_stdev']
        lmin <- bx[i,'lrr_mean'] - nsd * bx[i,'lrr_stdev']
        lmax <- bx[i,'lrr_mean'] + nsd * bx[i,'lrr_stdev']
        zmin <- max(zmin, 0.5)
        rect(zmin, lmin, zmax, lmax, border='black', lwd=1)       
    }   
}

# finish up
graphics.off()

# plot the changing CN states resulting from refinement
if(cn_col == 'CN_OUT'){
    label <- 'CN_IN.CN'  
    jpeg(paste(PLOT_PREFIX, cn_col, label, 'jpg', sep="."), 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    xylim <- c(-0.5,max(sd$CN_IN, sd$CN_OUT)+0.5) # be sure to include CN0 space, even if no probes
    plot(jitter(sd$CN_IN, 2), jitter(sd$CN_OUT, 2), pch=".",
         xlab='CN_IN', ylab='CN_OUT',
         xlim=xylim, ylim=xylim) # NAs not plotted, of course
    graphics.off()
    fracChanged <- nrow(d[d$IS_NA==0 & d$CN_OUT != d$CN_IN,]) / nUsable
    message(paste(fracChanged, "= fraction of probes that changed CN (not including any mosaic CNVs)"))
}

# memory management
rm(dd, sd)

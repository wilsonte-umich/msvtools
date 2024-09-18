
# initialize
message("plotting ZYG-LRR correlation")
for(plotCN in 0:5){

# downsample the data for the plot
# include CN=0, but exclude mosaics
# dd <- d[!is.na(d[[cn_col]]) & d[[cn_col]]%%1==0,]
dd <- d[!is.na(d[[cn_col]]) & d[[cn_col]] == plotCN,]
if(nrow(dd) == 0) next
nPoints <- 100000
sd <- if(nrow(dd) > nPoints){
    dd[sample(1:nrow(dd), nPoints),]
} else {
    dd
}

# plot the correlation
jpeg(paste(PLOT_PREFIX, cn_col, paste("cn", plotCN, sep = "_"), 'CORRELATION', 'jpg', sep="."), 
       width=width*2, height=height*2, unit=units, pointsize=pointsize, res=resol)
plot(sd$ZYG, sd$LRR, pch=".", col=cnColorsAlpha[plotCN + 1],
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

# close CN loop
}

# plot the changing CN states resulting from refinement
if(cn_col == 'CN_OUT'){
    dd <- d[!is.na(d[[cn_col]]) & d[[cn_col]]%%1==0] # thus, only integral CN out, not mosaics
    nPoints <- 25000
    sd <- if(nrow(dd) > nPoints){
        dd[sample(1:nrow(dd), nPoints)]
    } else {
        dd
    }
    label <- 'CN_IN.CN'  
    jpeg(paste(PLOT_PREFIX, cn_col, label, 'jpg', sep="."), 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    xylim <- c(-0.5,max(sd$CN_IN, sd$CN_OUT)+0.5) # be sure to include CN0 space, even if no probes
    plot(jitter(sd$CN_IN, a = 0.5), jitter(sd$CN_OUT, a = 0.5), pch=".",
         xlab='CN_IN', ylab='CN_OUT', col=cnColorsAlpha[1],
         xlim=xylim, ylim=xylim) # NAs not plotted, of course
    graphics.off()
    fracChanged <- dd[masked == FALSE & CN_OUT != CN_IN, .N] / dd[masked == FALSE, .N]
    message(paste(fracChanged, "= fraction of probes that changed CN (not including mosaic CNVs)"))
}

# memory management
rm(dd, sd)

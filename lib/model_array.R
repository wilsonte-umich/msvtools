# coordinate zygosity and LRR fitting and modeling

# generic function for plotting modeled probabilities without actual data
plotModeled <- function(CN, type, bins, n_bins, modeled, xlab, xlim){
    label <- paste('CN', CN, type, sep="_")    
    jpeg(paste(PLOT_PREFIX, cn_col, label, 'jpg', sep="."), 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    x <- bins/n_bins
    plot(x, modeled, type="l", col="forestgreen",
         xlab=xlab, ylab="Frequency",
         xlim=xlim)
    graphics.off()
}

# function to convert to integer bins
binPrm <- function(v, n_bins, unbin=FALSE){ 
    f <- ifelse(unbin, 1/n_bins, n_bins)
    v * ifelse(grepl("^L", names(v)), 1, f)
}

# determine data adequacy for modeling a CN state
CNs <- sort(unique(d[[cn_col]])) # the available model copy numbers (NA purged by sort)
isSuff <- function(CN, inf_flt){
    CN %in% CNs & (
        #cn_col == 'CN_OUT' | # always plot the final model
        nrow(d[arrayFittable & d[[cn_col]]==CN & inf_flt]) >= 1e3 # includes usability filter
    )
}

# apply a bad probe allowance (equally likely in all bins)
# and normalize emission probabilities to sum to 1 across all bins
# i.e. any probe must land in one of the target bins for every state
normalizeEPs <- function(eps, bad_prb_prob){
    eps <- (1-BAD_PROBE_FREQ) * eps +    # valid probe probabilities
           BAD_PROBE_FREQ * bad_prb_prob # bad probe probabilities
    eps / sum(eps)  
}

# declare lists to take emission probability lookup tables
# list index=CN, value=data.table of bins and EPs
zyg_eps <- list()
frc_zyg <- list()
lrr_eps <- list()

# function to catch errors and warnings (too complicated in R...)
# used for zygosity nls optimization
withWarningsErrors <- function(expr){
    myWarning <- NULL
    myError   <- NULL
    wHandler  <- function(w){ # warning handler, called by withCallingHandlers
        myWarning <<- w$message
        invokeRestart("muffleWarning")
    }
    list(
        value=withCallingHandlers(
            tryCatch(
                expr,
                error = function(e) myError <<- e$message # error handler, called by tryCatch
            ),
            warning=wHandler
        ),
        warning=myWarning, # returns list with expression value, warning and error
        error=myError      # latter two NULL on success
    )
}

# compare simple and complex nls curve fits (complex assumes input model deficiencies)
# return the best option
chooseFit <- function(dat, type, other, updt=NULL){
    message(paste("  ", type, ifelse(dat$isConv, 'convergent', 'NOT convergent')))
    if(!is.null(updt)){ # if preferring simple, give target values expected names
        for(OUT in names(updt)) dat[OUT] = dat[updt[[OUT]]]     
    }
    unlink(other$jpgFile) # discard the unused fit image
    zyg_eps <<- c(zyg_eps, dat$zyg_eps) # commit the chosen fit data 
    frc_zyg <<- c(frc_zyg, dat$frc_zyg)
    dat$type <- type
    dat
}
compareFit <- function(smp, cmp, updtSmp=NULL){
    if(smp$error & cmp$error){ # so far, this never happens...
        stop("NLS FAILURE: error with both simple and complex models")
    } else if(smp$error){
        chooseFit(cmp, 'complex', smp)
    } else if(cmp$error) {
        chooseFit(smp, 'simple',  cmp, updtSmp)
    } else if(smp$deviance < cmp$deviance){
        chooseFit(smp, 'simple',  cmp, updtSmp)
    } else {
        chooseFit(cmp, 'complex', smp)
    }
}

# function for projecting emission probabilities of mosaic states
# low is the smaller value, high is the large value, EP distribution
# low and high provided as identical length vectors of EP
# result is a EP distribution spread evenly between the low and high peaks
# that eclipses as the boundary distributions rise
#       alternative to a single mosaic state would be to precisely model CN2.1,2.2,etc
project_mosaic <- function(low, high, bad_prb_prob){
    cdf_low  <- cumsum(low)
    cdf_high <- cumsum(high)
    inv_low  <- 1 - cdf_low
    ep       <- pmax(0, 1 - pmax(inv_low, cdf_high)) # pmax(0) prevent tiny negative number overruns via prior math
    normalizeEPs(ep/sum(ep), bad_prb_prob) 
}

# find zygosity statistics across array
source(paste(LIB_DIR, "model_zygosity.R", sep="/"))

# use zygosity statistics to find LRR statistics across array
source(paste(LIB_DIR, "model_LRR.R", sep="/"))

# plot the ZYG-LRR correlation
source(paste(LIB_DIR, "plot_ZYG_LRR_correlation.R", sep="/"))


# this script compares each of a set of SVs pairwise across many grouped samples

# collect parameters
LIB_DIR       <- Sys.getenv('LIB_DIR')
DATAFILE      <- Sys.getenv('DATAFILE')    # input data, LRR, etc.
PLOT_DIR      <- Sys.getenv('PLOT_DIR')
SIZE_THRESHOLD  <- as.numeric(Sys.getenv('SIZE_THRESHOLD'))
SIZE          <- Sys.getenv('SIZE')
GROUP_NAME    <- Sys.getenv('GROUP_NAME')
SAMPLES_DIR   <- Sys.getenv('SAMPLES_DIR')
SAMPLES       <- Sys.getenv('SAMPLES')
RLL_COL       <- Sys.getenv('RLL_COL')
RLL_SMP_COL   <- Sys.getenv('RLL_SMP_COL')
CMP_COL       <- Sys.getenv('CMP_COL')
CMP_FILE      <- Sys.getenv('CMP_FILE')

# set parameters
SAMPLES     <- strsplit(SAMPLES, ",")[[1]]
N_SAMPLES   <- length(SAMPLES)
tgtSmpIs    <- 1:(N_SAMPLES - 1)
RLL_COL     <- strsplit(RLL_COL, ",")[[1]]
RLL_SMP_COL <- strsplit(RLL_SMP_COL, ",")[[1]]
CMP_COL     <- strsplit(CMP_COL, ",")[[1]]
MIN_PROBES_SV  <- 5
N_FLANK_PROBES <- 10 # same as find_SVs.R
N_END_PROBES   <- N_FLANK_PROBES * 2
PLOT_PREFIX    <- paste(PLOT_DIR, GROUP_NAME, sep="/")

# set common image properties
width     <- 2.5
height    <- 2.8
units     <- 'in' # w and h in inches
pointsize <- 8
resol     <- 900 #dpi
pch       <- "." #20
cex       <- 1 #0.3

# load the data
message("loading SV regions")
d <- read.table(DATAFILE, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# collect the source files and load all HMMs and SVs by sample
message("loading sample SVs")
files <- list()
SVs   <- list()
HMMs  <- list()
chroms  <- character()
n_probe <- numeric()
n_inf   <- numeric()
svrll   <- numeric()
addFile <- function(smp, type, ext){
    file <- paste(SAMPLES_DIR, '/msvtools.segment.', type, '.', smp, '.', ext, sep="")
    if(!file.exists(file)[1]) stop(paste("file not found:", file))
    files[[smp]][[type]] <<- file
}
for(smp in SAMPLES){
    files[[smp]] <- list()
    addFile(smp, 'probes', 'bed.bgz') # probes loaded region by region, below     
    addFile(smp, 'HMM', 'RData')
    load(files[[smp]][['HMM']])
    HMMs[[smp]] <- HMM
}
rm(HMM)
SVs[SAMPLES] <- lapply(SAMPLES, function(smp){
    addFile(smp, 'SVs', 'bed.bgz')   
    pipe <- paste("zcat", files[[smp]][['SVs']], "| sed 's/^#//'")
    SVs  <- read.table(pipe(pipe), header=TRUE, sep="\t", stringsAsFactors=FALSE)
    sz   <- SVs$END - SVs$START
    flt  <- if(SIZE=="large") sz > SIZE_THRESHOLD else sz <= SIZE_THRESHOLD
    SVs  <- SVs[flt, ]
    if(nrow(SVs)>0){
        SVs$REGION_ID      <- 0         
        SVs$QRLL           <- 0
        SVs$MIN_TRLL       <- 0
        SVs$MAX_TRLL       <- 0
        SVs$MIN_TRLL_SMP   <- 'NA'
        SVs$MAX_TRLL_SMP   <- 'NA'
        SVs$OVLP_SMP       <- 'NA'
    } else {
        SVs$REGION_ID      <- numeric()         
        SVs$QRLL           <- numeric()
        SVs$MIN_TRLL       <- numeric()
        SVs$MAX_TRLL       <- numeric()
        SVs$MIN_TRLL_SMP   <- character()
        SVs$MAX_TRLL_SMP   <- character()
        SVs$OVLP_SMP       <- character()
    }
    chroms  <<- c(chroms,  SVs$CHROM)
    n_probe <<- c(n_probe, SVs$N_PROBES)
    n_inf   <<- c(n_inf,   SVs$N_INF)
    svrll   <<- c(svrll,   SVs$REL_LL)
    SVs
})
SVs_OUT <- SVs[[1]][FALSE,] # create an empty dataframe for output

# report on the total set of SVs over all samples
print(aggregate(chroms, list(chroms), length))
message(paste(length(chroms), "totals SVs over all samples"))

# make plot correlating N_PROBES to REL_LL
plotCorr <- function(log){
    label <- paste('N_PROBES', 'REL_LL', if(log) 'log' else '', sep="_")
    jpeg(paste(PLOT_PREFIX, 'compare', label, 'jpg', sep="."), 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    xmax <- max(n_probe, na.rm=TRUE)
    ymax <- max(svrll, na.rm=TRUE)    
    plot(1, 1, type="n", pch=19, cex=0.25,
         xlab='No. of Probes',
         ylab='Rel. Log Likelihood',
         xlim=c(if(log) MIN_PROBES_SV else 0, xmax),
         ylim=c(if(log) min(svrll) else 0, ymax),    
         log=if(log) 'xy' else '',
         main=GROUP_NAME)
    for(i in 1:3){ # lines at different levels of Avg_RLL, i.e. RLL per probe
        relL <- 10 ^ i
        np <- MIN_PROBES_SV:xmax
        lines(np, log(relL)*np, lwd=1, col=i+1)
    }
    flt <- n_inf == 0
    points(n_probe[flt], svrll[flt], pch=19, cex=0.25, col="grey")    
    flt <- n_inf > 0
    points(n_probe[flt], svrll[flt], pch=19, cex=0.25)
    graphics.off()      
}
plotCorr(TRUE)

# make plot of histogram of average REL_LL per probe
avgRL <- sapply(exp(svrll / n_probe), function(x) min(x, 250))
label <- 'AVG_REL_L'
jpeg(paste(PLOT_PREFIX, 'compare', label, 'jpg', sep="."), 
       width=width, height=height, unit=units, pointsize=pointsize, res=resol)
hist(avgRL, breaks=25, xlab="Avg. Rel. Likelihood Per Probe")
graphics.off()
  
# functions to calculate relative log likelihood of
#   a specific set of a sample's probes (selected from a query sample's SV)
#   over the allelic states of the reference model
# and to create useful metrics from these values
source(paste(LIB_DIR, "run_HMM.R", sep="/"))
RLL <- function(smp, prs, qhs){              # prs and qhs correspond to same probe series
    flt <- which(prs$IS_NA==0 & !is.na(qhs)) # but some may be unusuable
    N   <- length(flt)
    if(N > 0){
        obs <- prs[flt,c('ot','os')]     
        RLL <- HMM_likelihood_obs(HMMs[[smp]], obs, qhs[flt]) -
               HMM_likelihood_obs(HMMs[[smp]], obs, prs[flt,'AS_IN'])         
    } else {
        NA     
    }
}
get_tgt_smps <- function(tgt_smps, ps, qhs){
    tlls <- sapply(tgt_smps, function(tgt_smp){
        RLL(tgt_smp, prbs[[tgt_smp]][ps,], qhs)
    })
    MIN_RLL <- min(tlls, na.rm=TRUE)   
    MAX_RLL <- max(tlls, na.rm=TRUE)
    MIN_I   <- which(tlls==MIN_RLL)[1]
    MAX_I   <- which(tlls==MAX_RLL)[1]    
    list(min_rll=MIN_RLL, min_smp=tgt_smps[MIN_I],
         max_rll=MAX_RLL, max_smp=tgt_smps[MAX_I])
}

# process one overlapping SV region at a time
message("comparing each SV to all other samples")
message(paste("processing", nrow(d), "SV regions"))

for(i in 1:nrow(d)){
#for(i in 1:5){   
    
    # get the region (may contain one or more SVs in one or more samples)
    rg <- d[i,]
    region <- paste(rg$CHROM, ":", rg$START, "-", rg$END, sep="")
    message(paste('    region ', rg$REGION_ID, ': ', region,
                  ' [', rg$N_SVS, ' SVs in ', rg$N_SAMPLES, ' samples]', sep=''))

    # collect all probes in all samples over the region
    prbs <- list()
    prbs[SAMPLES] <- lapply(SAMPLES, function(smp){
        pipe <- paste("tabix -h", files[[smp]][['probes']], region, "| sed 's/^#//'")
        p <- read.table(pipe(pipe), header=TRUE, sep="\t", stringsAsFactors=FALSE)
        flt <- p$IS_NA==1
        p[flt,'AS_IN']  <- 'NA' # nll_col
        p[flt,'AS_OUT'] <- 'NA' 
        p
    })
    
    # identify the samples (queries) with SVs in the region
    # process one sample at a time
    qry_smps <- strsplit(rg$SAMPLES, ",")[[1]]
    for(qry_smp in qry_smps){

        # get 'query' sample probes in the specific SV region (may include SV, ref, unusable)
        qps   <- prbs[[qry_smp]]
        n_qps <- nrow(qps)
        
        # pull the SVs for this sample in the region
        qsvs  <- SVs[[qry_smp]]
        qsvsi <- which(qsvs$CHROM       == rg$CHROM &
                       qsvs$FLANK_START >= rg$START &
                       qsvs$FLANK_END   <= rg$END) # no SV can ever go past edge of a region 
        for(i in qsvsi){
            qsv  <- qsvs[i,]
            
            # find the usable probes in the query SV span
            # plus some number of flanking probes (presumed to match the model)
            # to aid in finding the truly best query-target match (more effective with smaller SVs)
            qsvps  <- which(qps$POS >= qsv$START &                           
                            qps$POS <= qsv$END)     
            qsvlfsi  <- max(qsvps[1] - N_FLANK_PROBES, 1)        
            qsvrfei  <- min(qsvps[length(qsvps)] + N_FLANK_PROBES, nrow(qps))
            qsvfps   <- qsvlfsi:qsvrfei 
            qsvfprs  <- qps[qsvfps,]              
            qsvfhs   <- qsvfprs$hs # this and prs may still contain unusable flank probes
            QRLL     <- RLL(qry_smp, qsvfprs, qsvfhs)
            
            # assess every other 'target' sample for consistency with the query SV
            # regardless of whether the target had a call or not
            tgt_smps <- SAMPLES[SAMPLES!=qry_smp]
            sv_smps  <- get_tgt_smps(tgt_smps, qsvfps,  qsvfhs)            

            # determine all other samples that had a SV call overlapping the current SV
            ovlp_smps <- sapply(tgt_smps, function(tgt_smp){
                prs <- prbs[[tgt_smp]][qsvps,c('IS_NA','CNC','LOH')]
                prs <- prs[prs$IS_NA==0,]
                if(sum(prs$CNC != 0 | prs$LOH == 'yes') > 0){
                    tgt_smp
                } else {
                    NA    
                }
            })
            ovlp_smps <- paste(ovlp_smps[!is.na(ovlp_smps)], collapse=",")

            # fill and commit the updated SV row
            qsv[,RLL_COL] <-
                c(rg$REGION_ID,
                  sapply(c(QRLL, sv_smps$min_rll, sv_smps$max_rll), round, 3) )
            qsv[,RLL_SMP_COL] <- # do separately to maintain numeric vs. character in columns
                c(sv_smps$min_smp, sv_smps$max_smp, ovlp_smps)
            SVs_OUT <- rbind(SVs_OUT, qsv)
        }
    }
}

# sort the SVs from all samples into one table
message("sorting sample SVs together")
SVs_OUT <- SVs_OUT[order(SVs_OUT$REGION_ID, SVs_OUT$CHROM, SVs_OUT$START, SVs_OUT$END),]

# print final output file
message("print composite SVs file")
header <- paste(CMP_COL, collapse="\t")
header <- paste("#", header, sep="") # comment the header for bgzip
write(header, file=CMP_FILE)         # print it
pipe <- paste("awk 'BEGIN{OFS=\"\\t\"}{$2+=0;$3+=0;print $0}' >> ", CMP_FILE, sep="")    
write.table(                         # then append the data rows
    SVs_OUT,
    file=pipe(pipe),
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE
)


#SVs_OUT <- read.table(CMP_FILE, header=TRUE, sep="\t", comment.char="")


# function for plotting group SV correlations
message("plotting metric correlations by SV")
alab <- c(
    QRLL ='log10[ Rel. Log Likelihood]',   
    MAX_TRLL ='log10[ Max. Other RLL (SV) ]'
)
isLog <- c(
    QRLL =TRUE,
    MAX_TRLL =TRUE
)
alim <- list(
    QRLL =NULL,
    MAX_TRLL =NULL
)
adjust_vals <- function(col){
    v <- SVs_OUT[[col]]
    if(isLog[[col]]) {
        sign(v) * log10(abs(v))
    } else {
        l <- alim[[col]]
        ifelse(v<l[1], l[1], ifelse(v>l[2], l[2], v))
    }   
}
plotRLL <- function(xcol, ycol, v=NULL){
    label <- paste(xcol, ycol, sep="_")
    jpeg(paste(PLOT_PREFIX, 'compare', label, 'jpg', sep="."), 
           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
    x <- adjust_vals(xcol)
    y <- adjust_vals(ycol)
    plot(x, y, type="p", pch=19, cex=0.25,
         xlab=alab[[xcol]], ylab=alab[[ycol]],
         xlim=alim[[xcol]], ylim=alim[[ycol]],
         main=GROUP_NAME)
    abline(h=0)
    if(!is.null(v)) abline(v=v)
    graphics.off() 
}

# make useful correlation plots
plotRLL('QRLL',  'MAX_TRLL')







# the most essential metric for uniqueness of a query SV is:
#   max( [(TLL_THS - TLL_MHS) - (TLL_THS - TLL_QHS)] /
#        [(TLL_THS - TLL_MHS) + (TLL_THS - TLL_QHS)] )  over all target samples
# where:
#   TLL = target (other) sample log likelihood
#   THS = target sample hidden states, i.e. its most likely path (could be model, query, or anything)
#   MHS = model hidden states, i.e. the reference/null hypothesis path
#   QHS = quey sample hidden states, i.e. the candidate SV path
# given that TLL_THS is the most likely target path
# the following are the expected values of the metric:
#   THS == MHS = [0-big]/[0+big] = -1, thus indicating that QHS is unique to query sample
#   THS == QHS = [big-0]/[big+0] =  1, thus positive numbers are similar to query sample, and query SV is not unique
#   THS == ??? = [big-big]/[big+big] ~ 0

#logL <- function(smp, prs, hs){ # prs and hs correspond to same probe series
#    flt <- which(prs$IS_NA==0)  # but some may be unusuable in the subject sample
#    obs <- prs[flt,c('ot','os')]
#    hs  <- hs[flt]
#    round(HMM_likelihood(HMMs[[smp]], obs, hs), 1)
#}
#get_tgt_smps_ <- function(tgt_smps, vals){
#    MIN_VAL <- min(vals)   
#    MAX_VAL <- max(vals)
#    list(min_val=MIN_VAL, min_smp=tgt_smps[which(vals==MIN_VAL)[1]],
#         max_val=MAX_VAL, max_smp=tgt_smps[which(vals==MAX_VAL)[1]])
#}
#get_tgt_smps <- function(tgt_smps, ps, qhs){
#    flt <- which(!is.na(qhs)) # qhs may be from unusable query probes
#    ps  <- ps[flt]
#    qhs <- qhs[flt]
#    
#    
#    tlls <- as.data.frame(t(sapply(tgt_smps, function(tgt_smp){
#        tprs <- prbs[[tgt_smp]][ps,]
#        c( LL_THS = logL(tgt_smp, tprs, tprs$AS_OUT),
#           LL_QHS = logL(tgt_smp, tprs, qhs),
#           LL_MHS = logL(tgt_smp, tprs, tprs$AS_IN) )
#    })))
#    trll_smps <- get_tgt_smps_(tgt_smps, tlls$LL_QHS - tlls$LL_MHS)
#    ths_mhs   <- tlls$LL_THS - tlls$LL_MHS     
#    ths_qhs   <- tlls$LL_THS - tlls$LL_QHS
#    bias_smps <- get_tgt_smps_(tgt_smps, (ths_mhs - ths_qhs) / (ths_mhs + ths_qhs))
#    list(trll=trll_smps, bias=bias_smps)
#}


# the most essential metric for uniqueness of a query SV is:
#   max[ -((TLL_QHS - TLL_MHS)/NT) / ((QLL_QHS - QLL_MHS)/NQ) ] over all target samples
#   i.e. max[ -AVG_TRLL / AVG_QRLL ]
# where:
#   QLL = query sample log likelihood, i.e. the one with the candidate SV
#   TLL = target (other) sample log likelihood
#   QHS = quey sample hidden states, i.e. the candidate SV path
#   MHS = model hidden states, i.e. the reference/null hypothesis path
# using the average RLL adjusts for variable existense of unusable probes
# given that QRLL should always be positive, this negated ratio generally ranges from
#   -1 (equally good SV call in another sample)
#   +1 (specific to query sample), to
# however, it is an _unbounded_ ratio, and becomes volatile at small QRLL

# this metric is calculated for both the entire query SV
# and also for regions bracketing each endpoint

    #    SVs$QRLL           <- 0
    #    SVs$MIN_TRLL       <- 0
    #    SVs$MAX_TRLL       <- 0
    #    SVs$MIN_RATIO      <- 0
    #    SVs$MAX_RATIO      <- 0
    #    SVs$QERLL          <- 0
    #    SVs$MIN_TERLL      <- 0
    #    SVs$MAX_TERLL      <- 0
    #    SVs$MIN_ERATIO     <- 0
    #    SVs$MAX_ERATIO     <- 0
    #    SVs$MIN_TRLL_SMP   <- 'NA'
    #    SVs$MAX_TRLL_SMP   <- 'NA'
    #    SVs$MIN_TERLL_SMP  <- 'NA'
    #    SVs$MAX_TERLL_SMP  <- 'NA'
    #} else {
    #    SVs$REGION_ID      <- numeric()         
    #    SVs$QRLL           <- numeric()
    #    SVs$MIN_TRLL       <- numeric()
    #    SVs$MAX_TRLL       <- numeric()
    #    SVs$MIN_RATIO      <- numeric()
    #    SVs$MAX_RATIO      <- numeric()
    #    SVs$QERLL          <- numeric()
    #    SVs$MIN_TERLL      <- numeric()
    #    SVs$MAX_TERLL      <- numeric()
    #    SVs$MIN_ERATIO     <- numeric()
    #    SVs$MAX_ERATIO     <- numeric()
    #    SVs$MIN_TRLL_SMP   <- character()
    #    SVs$MAX_TRLL_SMP   <- character()
    #    SVs$MIN_TERLL_SMP  <- character()
    #    SVs$MAX_TERLL_SMP  <- character()


            ## find the usable non-reference-model probes in the query SV
            ## (re)-calculate QRLL
            #qsvps  <- which(qps$POS    >= qsv$START &                           
            #                qps$POS    <= qsv$END &                           
            #                qps$IS_NA  == 0 &
            #                qps$AS_OUT != qps$AS_IN)
            #qsvprs <- qps[qsvps,]
            #qsvhs  <- qsvprs$hs
            #QRLL   <- RLL(qry_smp, qsvprs, qsvhs)
            #
            ## calculate the relative log likelihood of the query
            ## over just the probes flanking each SV endpoint (to aid parsing of overlapping SVs)
            #qsvlfsi  <- max(qsvps[1] - N_FLANK_PROBES, 1)
            #qsvlfei  <- min(qsvlfsi+N_END_PROBES-1, n_qps)
            #qsvrfei  <- min(qsvps[length(qsvps)] + N_FLANK_PROBES, nrow(qps))
            #qsvrfsi  <- max(qsvrfei-N_END_PROBES+1, qsvps[1])            
            #qsveps   <- unique(c(qsvlfsi:qsvlfei, qsvrfsi:qsvrfei)) # indices of endpoint region probes
            #qsveprs  <- qps[qsveps,]            
            #qsvehs   <- qsveprs$hs # this and prs may still contain unusable flank probes
            #QERLL    <- RLL(qry_smp, qsveprs, qsvehs)
            #
            #
            #
            #
            #
            ## assess every other 'target' sample for consistency with the query SV
            ## regardless of whether the target had a call or not
            ## keep only the next closest matching SV
            #tgt_smps <- SAMPLES[SAMPLES!=qry_smp]
            #sv_smps  <- get_tgt_smps(tgt_smps, qsvps,  qsvhs)
            #ep_smps  <- get_tgt_smps(tgt_smps, qsveps, qsvehs)
            #
            ## fill and commit the updated SV row
            #qsv[,RLL_COL] <-
            #    c(rg$REGION_ID,
            #      sapply(c(QRLL[['RLL']],
            #               sv_smps$min_rll, sv_smps$max_rll, -sv_smps$min_avg/QRLL[['AVG_RLL']],  -sv_smps$max_avg/QRLL[['AVG_RLL']],
            #               QERLL[['RLL']],
            #               ep_smps$min_rll, ep_smps$max_rll, -ep_smps$min_avg/QERLL[['AVG_RLL']], -ep_smps$max_avg/QERLL[['AVG_RLL']]),
            #             round, 3) )
            #qsv[,RLL_SMP_COL] <- # do separately to maintain numeric vs. character in columns
            #    c(sv_smps$min_smp, sv_smps$max_smp,
            #      ep_smps$min_smp, ep_smps$max_smp)
            #SVs_OUT <- rbind(SVs_OUT, qsv)
#
## function for plotting group SV correlations
#message("plotting metric correlations by SV")
#alab <- c(
#    QRLL ='log10[ Rel. Log Likelihood (SV) ]',
#    QERLL='log10[ Rel. Log Likelihood (Ends) ]',    
#    MAX_TRLL ='log10[ Max. Other RLL (SV) ]',
#    MAX_TERLL='log10[ Max. Other RLL (Ends) ]',
#    MAX_RATIO ='Uniqueness (SV)',
#    MAX_ERATIO='Uniqueness (Ends)'
#)
#isLog <- c(
#    QRLL =TRUE,
#    QERLL=TRUE,
#    MAX_TRLL =TRUE,
#    MAX_TERLL=TRUE,
#    MAX_RATIO =FALSE,
#    MAX_ERATIO=FALSE
#)
#alim <- list(
#    QRLL =NULL,
#    QERLL=NULL,
#    MAX_TRLL =NULL,
#    MAX_TERLL=NULL,
#    MAX_RATIO =c(-2.5,2.5),
#    MAX_ERATIO=c(-2.5,2.5)   
#)
#adjust_vals <- function(col){
#    v <- SVs_OUT[[col]]
#    if(isLog[[col]]) {
#        sign(v) * log10(abs(v))
#    } else {
#        l <- alim[[col]]
#        ifelse(v<l[1], l[1], ifelse(v>l[2], l[2], v))
#    }   
#}
#plotRLL <- function(xcol, ycol, v=NULL){
#    label <- paste(xcol, ycol, sep="_")
#    bitmap(file=paste(PLOT_PREFIX, 'compare', label, 'jpg', sep="."), type='jpeg',
#           width=width, height=height, unit=units, pointsize=pointsize, res=resol)
#    x <- adjust_vals(xcol)
#    y <- adjust_vals(ycol)
#    plot(x, y, type="p", pch=19, cex=0.25,
#         xlab=alab[[xcol]], ylab=alab[[ycol]],
#         xlim=alim[[xcol]], ylim=alim[[ycol]],
#         main=GROUP_NAME)
#    abline(h=0)
#    if(!is.null(v)) abline(v=v)
#    graphics.off() 
#}
#
## make useful correlation plots
#plotRLL('QRLL',      'QERLL')
#plotRLL('QRLL',      'MAX_TRLL')
#plotRLL('QERLL',     'MAX_TERLL')
#plotRLL('QRLL',      'MAX_RATIO',   v=0)
#plotRLL('QERLL',     'MAX_ERATIO',  v=0)
#plotRLL('MAX_RATIO',  'MAX_ERATIO', v=0)
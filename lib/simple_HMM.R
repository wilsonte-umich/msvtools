# build and solve a Hidden Markov model for an array sample or set of clonally related samples
# these are generic functions, not array specific

HMM_init <- function(transProb, emissProbs){
    T <- nrow(emissProbs) # i.e., number of observations
    N <- ncol(emissProbs) # i.e., number of hidden states
    transProbs <- matrix(transProb, nrow = N, ncol = N)
    diag(transProbs) <- 1 - transProb * (N - 1)
    list(
        emissProbs = emissProbs,
        transProbs = log(transProbs)
    )
}
HMM_viterbi <- function(hmm, observations = TRUE){
    ep <- hmm$emissProbs[observations,] # is a matrix, not a data.table
    tp <- hmm$transProbs
    ep[is.na(ep)] <- -Inf # block unusable paths as having zero probability
    ep[apply(ep, 1, max) == -Inf, ] <- log(1) # mask unusable bins, i.e., those with no usable paths

    # 1. initialization (observation t=1)
    T          <- nrow(ep) # length of the sequence of observations
    N          <- ncol(ep) # number of states
    delta      <- log(matrix(0, nrow = T, ncol = N))
    delta[1, ] <- sapply(1:N, function(i) log(1 / N) + ep[1, i])
    phi        <- matrix(NA, nrow = T, ncol = N)

    # 2. recursion;
    # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
    for (t in 2:T){
        pt <- t - 1
        for (j in 1:N){     # j = this hs
            ep_j <- ep[t, j]
            for (i in 1:N){ # i = prev hs
                delta_ <- delta[pt, i] + tp[i, j] + ep_j
                if(delta[t, j] < delta_){
                    delta[t, j] <- delta_
                    phi[pt, j]  <- i
                }
            }
        }
    }
    
    # 3. termination
    prob <- -Inf
    hsi  <- rep(1, T)
    for (j in 1:N){
        if(prob < delta[T, j]){
            prob <- delta[T, j]
            hsi[T] <- j
        }
    }
    
    # 4. reconstruction and return the hidden state indices
    for (t in (T - 1):1) hsi[t] <- phi[t, hsi[t + 1]]
    hsi
}


HMM_weighted.median <- function (x, w) {
    flt <- !is.na(x) & !is.na(w)
    x   <- x[flt]
    w   <- w[flt]
    if(length(w)==0 | sum(w)==0){
        0
    } else {    
        o   <- order(x)
        x   <- x[o]
        w   <- w[o]
        p   <- cumsum(w)/sum(w)
        n   <- sum(p < 0.5)
        if (p[n + 1] > 0.5){
            x[n + 1]
        } else {
            (x[n + 1] + x[n + 2])/2
        }         
    }
}
HMM_set_median_prob <- function(HMM){ # used to estimate expected likelihoods
    N <- ncol(ep) # number of states
    HMM$wm <- log(array(0,dim=c(HMM$ot$N,N,N)))
    for(ot in 1:HMM$ot$N){
        for(hs_in in 1:N){ # weight by the input hs
            w <- exp(HMM$emissProbs[ot,hs_in,])
            for(hs_out in 1:N){ # use ep from the output hs
                p <- exp(HMM$emissProbs[ot,hs_out,])
                HMM$wm[ot,hs_in,hs_out] <- log(HMM_weighted.median(p, w))
            }
        }
    }
    HMM
}

#-------------------------------------------------------------------------------
# likelihood of specific path, using all parameters of the model
#-------------------------------------------------------------------------------
# find the likelihood of a specific path given a set of observations and a HMM
HMM_likelihood <- function(HMM, obs, hs){ # obs provided as indices, not categories
    T <- length(obs$ot) 
    if(length(obs$os) != T) stop("different numbers of observations provided for ot vs. os")
    if(length(hs)     != T) stop("different numbers of observations provided for ot vs. hs")
    if(is.character(hs)) hs <- HMM$hs$i[hs] # hs may be indices or categories
    t <- 1 # initialize first observation based on starting probabilities 
    prob <- HMM$sp[hs[t]] + HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]
    for (t in 2:T){
        prob <- prob + HMM$tp[hs[t-1],hs[t]] + HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]
    }
    prob
}
#-------------------------------------------------------------------------------
# observed and expected likelihoods of a specific set of hs
# using only ep, not tp
#-------------------------------------------------------------------------------
HMM_likelihood_init <- function(env){
    T <- length(env$obs$ot) 
    if(length(env$obs$os) != T) stop("different numbers of observations provided for ot vs. os")
    if(length(env$hs)     != T) stop("different numbers of observations provided for ot vs. hs")
    if(is.character(env$hs)) env$hs <- env$HMM$hs$i[env$hs] # hs may be indices or categories
    T
}
HMM_likelihood_obs <- function(HMM, obs, hs){ # likelihood IGNORING transition probabilities
    T <- HMM_likelihood_init(environment())
    sum(sapply(1:T, function(t) HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]))
}
HMM_likelihood_exp <- function(HMM, obs, hs_in, hs){
    T <- HMM_likelihood_init(environment())
    if(is.character(hs_in)) hs_in <- HMM$hs$i[hs_in]
    lvl <- aggregate(1:T, list(ot=obs$ot, hs_in=hs_in, hs_out=hs), length)
    sum(sapply(1:nrow(lvl), function(i){
        HMM$wm[lvl[i,'ot'],lvl[i,'hs_in'],lvl[i,'hs_out']] * lvl[i,'x']        
    }))
}
#===============================================================================

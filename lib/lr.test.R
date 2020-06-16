
# adapted from https://cran.r-project.org/web/packages/extRemes
lr.test <- function(l1, l2, alpha=0.05, df=1, ...) {
    STATISTIC        <- 2 * (l2 - l1)
    #CRITVAL          <- qchisq(alpha, df=df, lower.tail=FALSE)    
    PVAL             <- pchisq(STATISTIC, df=df, lower.tail=FALSE)
    PVAL
    #list(alternative="greater", p.value=PVAL)
}

#null == l1 should be the model with fewer parameters = FORWARD
#alt  == l2 should be the mode with more parameters, has the same or greater log-likelihood = SV PATH
#df = n probes

#The test is only valid for comparing nested models. That is, the parameters of
#one model must be a subset of the parameters of the second model.
#
#Suppose the base model, m0, is nested within the model m1. Let x be the
#negative log-likelihood for m0 and y for m1. Then the likelihood-ratio statistic
#(or deviance statistic) is given by (Coles, 2001, p 35; Reiss and Thomas, 2007, p 118):
#
#D = -2*(y - x).
#
#Letting c.alpha be the (1 - alpha) quantile of the chi-square distribution with
#degrees of freedom equal to the difference in the number of model parameters,
#the null hypothesis that D = 0 is rejected if D > c.alpha (i.e., in favor of model m1).
# https://github.com/RemkoDuursma/bootpredictlme4/blob/master/R/predict.merMod.R

#' A predict method for merMod models with the bootstrap
#'@description When this package is loaded after loading \code{lme4}, it replaces the predict method for linear mixed-effects models (\code{merMod} objects) with this function, which provides the argument \code{se.fit} (just like in \code{\link{predict.lm}}). It uses the (semi-)parametric bootstrap, as implemented in \code{\link{bootMer}} to estimate standard errors and confidence intervals for the predictions. Also useful in conjunction with the \code{\link{visreg}} package, which by default does not plot confidence intervals for \code{merMod} models. 
#'
#'Note that this package cannot be used in conjunction with \code{lmerTest}, as that package also replaces some \code{lme4} methods. 
#'
#'Note the argument \code{re.form} in \code{\link{predict.merMod}}, which controls how random effects are incorporated in the predictions. This argument can be set via this function as well (see Examples).
#'@param object An object as returned by \code{\link{lmer}}
#'@param nsim The number of bootstrap replicates. The default is a small number to allow quick testing. A warning is printed when nsim < 100, and it can be set via \code{\link{options(bootnsim=1000)}} (see Examples).
#'@param se.fit If TRUE, returns standard error (se.fit) and confidence interval (ci.fit) for the predictions, as components of the returned list.
#'@param alpha Controls the coverage of the confidence interval. 
#'@param \dots Further arguments passed to \code{\link{bootMer}}.
#'@details Two standard errors are calculated. The first (\code{se.fit}) is the effective standard error, consistent with the asymptotic confidence interval. This standard error is calculated because other methods may estimate a confidence interval based on a standard error, which is inaccurate when using the simple standard error (\code{se.boot}). For example, when the desired coverage of the confidence interval is 0.95, giving a normal quantile of 1.96, the \code{se.boot} is calculated as half confidence interval width divided by 1.96.This method is applied especially because \code{\link{visreg}} calculates a confidence interval assuming normality. 
#'
#'The second ('se.boot') is imply the standard deviation of the bootstrap replicates (i.e. the standard bootstrapped standard error). 
#'
#'@value When \code{se.fit=FALSE} (the default), simply calls \code{\link{predict.merMod}}. If \code{se.fit=TRUE}, invokes the bootstrap, and returns a list with components \code{fit} (fitted value), \code{se.fit} (bootstrapped standard error), \code{se.boot} (effective standard error, reproducing the confidence interval with the assumption of normality), \code{ci.fit} (the bootstrapped confidence interval).
#'@export
#'@examples 
#'
#'# Fit a linear mixed-effects model 
#'fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#'# Predictions without standard error (fixed effects only)
#'predict(fm1, newdata=data.frame(Days=5), re.form=NA)
#'
#'# Predictions with standard error and confidence interval
#'# Set number of bootstrap replicates first (you should use a larger number)
#'options(bootnsim=20)
#'predict(fm1, newdata=data.frame(Days=5), re.form=NA, se.fit=TRUE)
#'
#'# Add confidence intervals to visreg
#'library(visreg)
#'visreg(fm1, "Days")
#'
#'# Also works with an overlay plot, using visreg
#'# First add an artificial group to the sleepstudy data
#'high <- with(sleepstudy, levels(reorder(Subject,Reaction,mean)))[1:9]
#'sleepstudy$Group <- factor(ifelse(sleepstudy$Subject %in% high, "A", "B"))
#'fm2 <- lmer(Reaction ~ Days*Group + (Days | Subject), sleepstudy)
#'visreg(fm2, "Days", by="Group", overlay=TRUE)
#'
predict_lmer <- function(object, nsim=getOption("bootnsim"), se.fit=FALSE, alpha=0.05, ...){
  
  if(se.fit){
    if(is.null(nsim))nsim <- 10
    if(nsim < 100)message("Number of bootstrap replicates very low. \n Set to higher value with e.g. options(bootnsim = 500)")
    
    b <- bootMer(object, FUN=function(x)lme4:::predict.merMod(x, ...), nsim=nsim)
    serr <- apply(b$t, 2, sd)
    ci <- apply(b$t,2, quantile,probs=c(alpha/2, 1 - alpha/2))
    
    qn <- qnorm(1 - alpha/2)
    se_eff <- apply(sweep(ci, 2, b$t0) / qn, 2, function(x)mean(abs(x)))
    
    return(list(fit=b$t0, se.fit=se_eff, se.boot=serr, ci.fit=ci))
  } else {
    return(lme4:::predict.merMod(object, ...))
  }
  
}

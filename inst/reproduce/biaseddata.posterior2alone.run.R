library(bettertogether)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
setmytheme()
registerDoParallel(cores = detectCores())
set.seed(16)

## data
load("biaseddata_dataset.RData")

y1 <- matrix(synthetic_dataset$y1, ncol = 1)
y2 <- matrix(synthetic_dataset$y2, ncol = 1)
#
### Model
## Module 1
module1 <- biaseddata_module1()
hyper1 <- module1$hyper1
## Module 2
module2 <- biaseddata_module2()
hyper2 <- module2$hyper2
#
hyper1$theta1_prior_sd <- 1
hyper2$theta2_prior_sd <- 0.1
## Algorithmic parameters
param_algo <- list(nthetas = 1024, minimum_diversity = 0.8, nmoves = 10, proposal = mixture_rmixmod())
#
doRun <- TRUE
###
### Posterior in Module 2 alone (ignoring Y1)
target_2 <- list(thetadim = 2, ydim = 1,
                 rprior = function(n, ...){
                   th <- matrix(ncol = 2, nrow = n)
                   th[,1] <- rnorm(n, mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd)
                   th[,2] <- rnorm(n, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd)
                   return(th)
                 },
                 dprior = function(thetas, ...){
                   dnorm(thetas[,1], mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd, log = TRUE) +
                     dnorm(thetas[,2], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)
                 },
                 fullloglikelihood = function(thetas, ys, ...){
                   return(module2$loglik(thetas[,1], thetas[,2], ys[,1]))
                 },
                 conditionallikelihood = function(thetas, ys, idata, ...){
                   return(module2$loglik(thetas[,1], thetas[,2], ys[idata,1]))
                 },
                 parameters = list())
rep <- 5
filename <- paste0("biaseddata_module2alone.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
if (doRun){
  results2alone <- foreach(i = 1:rep) %dorng% {
    res <- smcsampler(y2, target_2, param_algo)
    list(thetas = res$thetas_history[[n2+1]], normw = res$normw_history[[n2+1]], logevidence = res$logevidence)
  }
  save(results2alone, file = filename)
} else {
  load(filename)
}


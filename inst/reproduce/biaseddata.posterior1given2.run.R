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
# whether to compute things or load pre-computed things
doRun <- TRUE
rep <- 5
###
### Posterior in Module 1 given Module 2

filename <- paste0("biaseddata_module2alone.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
load(filename)

## then initialize from the particles of the previous SMC, and assimilate Module 1
target_1given2 <- list(thetadim = 2, ydim = 1,
                       dprior = function(thetas, parameters){
                         dnorm(thetas[,1], mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd, log = TRUE) +
                           dnorm(thetas[,2], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE) +
                           module2$loglik(thetas[,1], thetas[,2], parameters$y2)
                       },
                       fullloglikelihood = function(thetas, ys, ...){
                         return(module1$loglik(thetas[,1], ys[,1]))
                       },
                       conditionallikelihood = function(thetas, ys, idata, ...) module1$loglik(thetas[,1], ys[idata,1]),
                       parameters = list(y2 = y2))


filename <- paste0("biaseddata_module1givenmodule2.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
if (doRun){
  results1given2 <- foreach(i = 1:rep) %dorng% {
    thetas_given2 <- results2alone[[i]]$thetas
    normw_given2 <- results2alone[[i]]$normw
    res <- smcsampler(y1, target_1given2, param_algo, thetas_given2, normw_given2)
    list(thetas = res$thetas_history[[n1+1]], normw = res$normw_history[[n1+1]], logevidence = res$logevidence)
  }
  save(results1given2, file = filename)
} else {
  load(filename)
}

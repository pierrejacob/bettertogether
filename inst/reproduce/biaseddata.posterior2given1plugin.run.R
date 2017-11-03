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
### Posterior in Module 1 given Module 1
filename <- paste0("biaseddata_module1.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
load(filename)
### Plug-in approach
target_plugin <- list(thetadim = 1, ydim = 1,
                      rprior = function(n, ...){
                        return(matrix(rnorm(n, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd), ncol = 1))
                      },
                      dprior = function(thetas, ...){
                        dnorm(thetas[,1], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)
                      },
                      fullloglikelihood = function(thetas, ys, parameters = NULL){
                        return(module2$loglik(parameters$theta1hat, thetas[,1], ys[,1]))
                      },
                      conditionallikelihood = function(thetas, ys, idata, parameters = NULL){
                        return(module2$loglik(parameters$theta1hat, thetas[,1], ys[idata,1]))
                      },
                      parameters = list())

filename <- paste0("biaseddata_module2givenmodule1plugin.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
if (doRun){
  results2given1plugin <- foreach(i = 1:rep) %dorng% {
    thetas1 <- results1[[i]]$thetas
    normw1 <- results1[[i]]$normw
    target_plugin$parameters$theta1hat <- wmean(thetas1, normw1)
    res <- smcsampler(y2, target_plugin, param_algo)
    list(thetas2 = res$thetas_history[[n2+1]], normw2 = res$normw_history[[n2+1]], sumlogevidence = sum(res$logevidence),
         theta1hat = target_plugin$parameters$theta1hat)
  }
  save(results2given1plugin, file = filename)
} else {
  load(file = filename)
}

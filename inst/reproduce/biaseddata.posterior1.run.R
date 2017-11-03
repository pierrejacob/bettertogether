library(bettertogether)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
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
### Posterior in Module 1
target_module1 <- list(thetadim = 1, ydim = 1,
                       rprior = function(n, ...) matrix(rnorm(n, mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd), ncol = 1),
                       dprior = function(thetas, ...) dnorm(thetas[,1], mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd, log = TRUE),
                       conditionallikelihood = function(thetas, ys, idata, ...) module1$loglik(thetas[,1], ys[idata,1]),
                       fullloglikelihood = function(thetas, ys, ...) module1$loglik(thetas[,1], ys[,1]),
                       parameters = list())

# whether to compute things or load pre-computed things
doRun <- TRUE
rep <- 5
filename <- paste0("biaseddata_module1.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")

if (doRun){
  results1 <- foreach(i = 1:rep) %dorng% {
    res <- smcsampler(y1, target_module1, param_algo)
    list(thetas = res$thetas_history[[n1+1]], normw = res$normw_history[[n1+1]], logevidence = res$logevidence)
  }
  save(results1, file = filename)
} else {
  load(filename)
}


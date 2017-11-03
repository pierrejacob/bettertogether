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
filename <- paste0("biaseddata_module1.N", param_algo$nthetas, ".K10.rep", rep, ".RData")
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


# based on N1 samples
N1 <- 1024
filename <- paste0("biaseddata_module2givenmodule1cut.N", param_algo$nthetas, ".K", param_algo$nmoves, ".M", N1, ".rep", rep, ".RData")

if (doRun){
  results2given1cut <- list()
  scores2given1cut <- list()
  for (irep in 1:rep){
    print(irep)
    thetas1 <- results1[[irep]]$thetas
    normw1 <- results1[[irep]]$normw
    thetas1reduced <- thetas1[systematic_resampling_n(normw1, N1),,drop=F]
    results_cut <- foreach (i = 1:N1) %dorng% {
      target_plugin$parameters$theta1hat <- thetas1reduced[i,]
      res <- smcsampler(y2, target_plugin, param_algo)
      thetas2 = res$thetas_history[[n2+1]]
      normw2 = res$normw_history[[n2+1]]
      list(thetas2 = thetas2,
           normw2 = normw2,
           theta1 = thetas1reduced[i,], logevidence = res$logevidence)
    }
    log_score <- rep(0, n2)
    for (time in 1:n2){
      logevid <- rep(0, N1)
      for (j in 1:N1){
        logevid[j] <-   results_cut[[j]]$logevidence[time]
      }
      maxlogevid <- max(logevid)
      log_score[time] <- maxlogevid + log(mean(exp(logevid - maxlogevid)))
    }
    for (j in 1:N1){
      results_cut[[j]]$logevidence <- NULL # (to save memory space)
    }
    scores2given1cut[[irep]] <- sum(log_score)
    results2given1cut[[irep]] <- results_cut
    save(results2given1cut, scores2given1cut, file = filename)
  }
} else {
  load(file = filename)
}


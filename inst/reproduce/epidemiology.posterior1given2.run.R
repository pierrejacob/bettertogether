library(bettertogether)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
setmytheme()
registerDoParallel(cores = detectCores())
set.seed(16)
# nhpv considered as Y
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
y1 <- matrix(nhpv, nrow = 1)
# Npart is put in the parameters
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
           143, 229, 696, 93)
J <- 13
# For module 2, ncases considered data
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
y2 <- matrix(ncases, nrow = 1)
# Npop considered parameters
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)
Npop_normalized <- log(10**(-3) * Npop)
#
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv

# Find parameters given Y2
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))

dprior <- function(thetas,...){
  theta1s <- thetas[,1:13]
  theta2s <- thetas[,14:15]
  evals <- sapply(1:nrow(theta2s), function(index)
    sum(dnorm(theta2s[index,], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
  for (j in 1:J){
    evals <- evals + dbeta(theta1s[,j], shape1 = 1, shape2 = 1, log = TRUE)
  }
  which.finiteprior <- which(is.finite(evals))
  if (length(which.finiteprior) > 0){
    evals[which.finiteprior] <- evals[which.finiteprior] +
      plummer_module2_loglikelihood(theta1s[which.finiteprior,,drop=F], theta2s[which.finiteprior,,drop=F], ncases, Npop_normalized)
  }
  return(evals)
}

# loglikelihood
loglik1 <- function(thetas, ys, parameters){
  if (is.null(nrow(thetas))) thetas <- matrix(thetas, nrow = 1)
  theta1s <- thetas[,1:13]
  return(sapply(1:nrow(theta1s), function(index) sum(dbinom(ys[1,], size = parameters$Npart, prob = theta1s[index,], log = TRUE))))
}

### Posterior in Module 1
target_module1given2 <- list(thetadim = 15, ydim = 1,
                             dprior = dprior,
                             fullloglikelihood = loglik1,
                             conditionallikelihood = function(thetas, ys, idata, parameters){
                               return(loglik1(thetas, ys, parameters))
                             },
                             parameters = list(Npart = Npart))

filename <- "epidemiology_module2alone.N1024.K100.rep5.RData"
load(file = filename)
param_algo <- list(nthetas = 1024, minimum_diversity = 0.8, nmoves = 100, proposal = mixture_rmixmod())
rep <- 5
filename <- paste0("epidemiology_module1givenmodule2.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
results1given2 <- foreach (irep = 1:rep) %dorng% {
  thetas_given2 <- results2alone[[irep]]$thetas
  normw_given2 <- results2alone[[irep]]$normw
  res <- smcsampler(y1, target_module1given2, param_algo, thetas_given2, normw_given2)
  list(thetas = res$thetas_history[[2]], normw = res$normw_history[[2]], logevidence = res$logevidence)
}
save(results1given2, file = filename)
load(file = filename)


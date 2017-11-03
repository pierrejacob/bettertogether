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
# generate from prior
rprior1 <- function(n,...){
  return(matrix(rbeta(n*J, shape1 = 1, shape2 = 1), ncol = 13))
}
# logdensity of prior
dprior1 <- function(thetas,...){
  return(sapply(1:nrow(thetas), function(index) sum(dbeta(thetas[index,], shape1 = 1, shape2 = 1, log = TRUE))))
}
# loglikelihood
loglik1 <- function(thetas, ys, parameters){
  if (is.null(nrow(thetas))) thetas <- matrix(thetas, nrow = 1)
  return(sapply(1:nrow(thetas), function(index) sum(dbinom(ys[1,], size = parameters$Npart, prob = thetas[index,], log = TRUE))))
}

#
param_algo <- list(nthetas = 2^10, minimum_diversity = 0.8, nmoves = 100, proposal = mixture_rmixmod())
#
### Posterior in Module 1
target_module1 <- list(thetadim = J, ydim = 1,
                       rprior = rprior1,
                       dprior = dprior1,
                       fullloglikelihood = loglik1,
                       conditionallikelihood = function(thetas, ys, idata, parameters) loglik1(thetas, ys[idata,,drop=F], parameters),
                       parameters = list(Npart = Npart))
# number of repeats
rep <- 5
filename <- paste0("epidemiology_module1altogether.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
results1 <- foreach(irep = 1:rep) %dorng% {
  res <- smcsampler(y1, target_module1, param_algo)
  list(thetas = res$thetas_history[[2]], normw = res$normw_history[[2]], logevidence = res$logevidence)
}
save(results1, file = filename)




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
target2alone <- list(thetadim = 15, ydim = 1,
                     rprior = function(n, ...){
                       theta1s <- matrix(0, nrow = n, ncol = 13)
                       for (j in 1:13){
                         theta1s[,j] <- rbeta(n, 1, 1)
                       }
                       theta2s <- matrix(0, nrow = n, ncol = 2)
                       for (j in 1:2){
                         theta2s[,j] <- rnorm(n, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd)
                       }
                       return(cbind(theta1s, theta2s))
                     },
                     dprior = function(thetas, ...){
                       theta1s <- thetas[,1:13]
                       theta2s <- thetas[,14:15]
                       evals <- sapply(1:nrow(theta2s), function(index)
                         sum(dnorm(theta2s[index,], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
                       for (j in 1:J){
                         evals <- evals + dbeta(theta1s[,j], shape1 = 1, shape2 = 1, log = TRUE)
                       }
                       return(evals)
                     },
                     fullloglikelihood = function(thetas, ys, ...){
                       return(plummer_module2_loglikelihood(thetas[,1:13,drop=F], thetas[,14:15,drop=F], ncases, Npop_normalized))
                     },
                     # here it's dummy; we pretend there's only one observation
                     conditionallikelihood = function(thetas, ys, idata, ...){
                       return(plummer_module2_loglikelihood(thetas[,1:13,drop=F], thetas[,14:15,drop=F], ncases, Npop_normalized))
                     },
                     parameters = list())

param_algo <- list(nthetas = 2^10, minimum_diversity = 0.8, nmoves = 100, proposal = mixture_rmixmod())
rep <- 5
filename <- paste0("epidemiology_module2alone.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
results2alone <- foreach(irep = 1:rep) %dorng% {
  res <- smcsampler(y2, target2alone, param_algo)
  list(thetas = res$thetas_history[[2]], normw = res$normw_history[[2]], logevidence = res$logevidence)
}
save(results2alone, file = filename)
load(file = filename)

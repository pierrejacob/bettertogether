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

dprior2 <- function(theta2s, hyper2){
  return(sapply(1:nrow(theta2s), function(index)
    sum(dnorm(theta2s[index,], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE))))
}
rprior2 <- function(n, hyper2){
  theta2s <- matrix(0, nrow = n, ncol = 2)
  for (j in 1:2){
    theta2s[,j] <- rnorm(n, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd)
  }
  return(theta2s)
}

#
# posterior in the first module is given by
# 13 Beta distributions with parameters
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv
##
param_algo <- list(nthetas = 2^10, minimum_diversity = 0.8, nmoves = 100, proposal = mixture_rmixmod())
target2given1 <- list(thetadim = 15, ydim = 1,
                      dprior = function(thetas, ...){
                        theta1s <- thetas[,1:13]
                        theta2s <- thetas[,14:15]
                        evals <- dprior2(theta2s, hyper2)
                        for (j in 1:J){
                          evals <- evals + dbeta(theta1s[,j], shape1 = posterior_phi_alpha[j], shape2 = posterior_phi_beta[j], log = TRUE)
                        }
                        return(evals)
                      },
                      fullloglikelihood = function(thetas, ys, ...){
                        return(plummer_module2_loglikelihood(thetas[,1:13], thetas[,14:15], ncases, Npop_normalized))
                      },
                      conditionallikelihood = function(thetas, ys, idata, ...){
                        return(plummer_module2_loglikelihood(thetas[,1:13], thetas[,14:15], ncases, Npop_normalized))
                      },
                      parameters = list())

rep <- 5
filename <- paste0("epidemiology_module2givenmodule1.N", param_algo$nthetas, ".K", param_algo$nmoves, ".rep", rep, ".RData")
results2given1 <- foreach(irep = 1:rep) %dorng% {
  theta1 <- matrix(0, nrow = param_algo$nthetas, ncol = J)
  for (j in 1:J){
    theta1[,j] <- rbeta(n = param_algo$nthetas, shape1 = posterior_phi_alpha[j], shape2 = posterior_phi_beta[j])
  }
  theta2 <- rprior2(param_algo$nthetas, hyper2)
  # augment the particles with theta2
  thetas_prior <- cbind(theta1, theta2)
  res <- smcsampler(y2, target2given1, param_algo, thetas_prior, rep(1/param_algo$nthetas, param_algo$nthetas))
  list(thetas = res$thetas_history[[2]], normw = res$normw_history[[2]], logevidence = res$logevidence)
}
save(results2given1, file = filename)

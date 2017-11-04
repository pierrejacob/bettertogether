# this script samples naively from the cut distribution
# in Plummer's epidemiology example
library(bettertogether)
rm(list = ls())
setmytheme()
registerDoParallel(cores = detectCores())
set.seed(16)

## This is the data
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

## The first posterior is conjugate, and is a Beta distribution
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv

# sample from the posterior in first module
nsamples <- 2^8
thetas1 <- foreach(isample = 1:nsamples, .combine = rbind) %dorng% {
  rbeta(J, posterior_phi_alpha, posterior_phi_beta)
}

# Now, we can approximate theta2 given each theta1, using any Monte Carlo method we like
# the following code does it with SMC samplers, but it could be a simpler Metropolis-Hastings algorithm
# as in naivecut-with-mcmc.R in the inst/tutorial/ folder

# prior in second module
# hyper parameters
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
# log prior density
dprior2 <- function(theta2s, hyper2){
  return(sapply(1:nrow(theta2s), function(index)
    sum(dnorm(theta2s[index,], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE))))
}
# sample from prior
rprior2 <- function(n, hyper2){
  theta2s <- matrix(0, nrow = n, ncol = 2)
  for (j in 1:2){
    theta2s[,j] <- rnorm(n, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd)
  }
  return(theta2s)
}
# define "target 2", i.e. posterior of theta2 given Y2 and theta1
target2given1plugin <- list(thetadim = 2, ydim = 1,
                            rprior = function(n, ...){
                              return(rprior2(n, hyper2))
                            },
                            dprior = function(thetas, ...){
                              evals <- dprior2(thetas, hyper2)
                              return(evals)
                            },
                            fullloglikelihood = function(thetas, ys, parameters){
                              return(plummer_module2_conditional(parameters$theta1hat, thetas, ncases, Npop_normalized))
                            },
                            conditionallikelihood = function(thetas, ys, idata, parameters){
                              return(plummer_module2_conditional(parameters$theta1hat, thetas, ncases, Npop_normalized))
                            },
                            parameters = list(theta1hat = NULL))

# algorithmic parameters for the SMC samplers
param_algo <- list(nthetas = 2^8, minimum_diversity = 0.5, nmoves = 1, proposal = mixture_rmixmod())
# if you want more accurate results, try this setting, which takes > 50 times longer to run:
# param_algo <- list(nthetas = 2^10, minimum_diversity = 0.8, nmoves = 10, proposal = mixture_rmixmod())

# then run SMC sampler for each theta1
# and get one theta2 out of each SMC sampler
thetas2cut <- foreach (isample = 1:nsamples, .combine = rbind) %dorng% {
  target2given1plugin$parameters$theta1hat <- thetas1[isample,]
  res <- smcsampler(y2, target2given1plugin, param_algo)
  thetas2 = res$thetas_history[[2]]
  normw2 = res$normw_history[[2]]
  index <- sample(1:nrow(thetas2), size = 1, prob = normw2)
  thetas2[index,]
}

#
plot(thetas2cut[,1], thetas2cut[,2])
hist(thetas2cut[,1])
hist(thetas2cut[,2])

# See the file unbiasedcut.R  in inst/tutorial 
# for a more sophisticated approach to sampling from the cut distribution

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

# prior for module 2
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
dprior2 <- function(theta2, hyper2){
  return(sum(dnorm(theta2, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
}

# Metropolis-Hastings kernel targeting the posterior of
# theta2 given theta1 and Y2, using a random walk proposal
# with covariance "Sigma_proposal"
get_MH_module2 <- function(theta1, Sigma_proposal = diag(c(.1,.1), 2, 2)){
  target <- function(x) plummer_module2_conditional(theta1, matrix(x, nrow = 1), ncases, Npop_normalized) + dprior2(x, hyper2)
  dimension <- 2
  # Markov kernel of MH
  MH_kernel <- function(chain_state){
    proposal_value <- chain_state + fast_rmvnorm(1, rep(0, dimension), Sigma_proposal)[1,]
    proposal_pdf <- target(proposal_value)
    current_pdf <- target(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(proposal_value)
    } else {
      return(chain_state)
    }
  }
  return(MH_kernel)
}

# First, we run the chain conditional on the posterior mean in the first module
# denoted by theta1hat
theta1hat <- colMeans(thetas1)
# define Metropolis Hastings kernel
kernel <- get_MH_module2(theta1hat, Sigma_proposal = diag(c(.1,.1), 2, 2))
# initial state of the chain
chain_state <- fast_rmvnorm(1, mean = c(0, 0), diag(1, 2, 2))[1,]
# number of MH iterations
niterations <- 10000
# run MH
chain <- matrix(nrow = niterations, ncol = 2)
chain[1,] <- chain_state
for (iteration in 2:niterations){
  chain_state <- kernel(chain_state)
  chain[iteration,] <- chain_state
}
# chain is of dimension niterations x 2
matplot(chain, type = "l")

# this chain approximates the distribution of theta2 given theta1hat and Y2
# i.e. a plug-in approximation
# let's compute its mean and variance, discarding first 1000 iterations as burn-in
theta2plugin_mean <- colMeans(chain[1001:niterations,])
theta2plugin_variance <- cov(chain[1001:niterations,])

# Now we are going to approximate the cut, by running MH chains of length 2000
# for each theta1 from the first posterior
# and retaining only the last state of each chain
# we'll use theta2plugin_variance for the covariance of the MH proposal
# and theta2plugin_mean for the initial position of the chain

niterations <- 2000
theta2s <- foreach(itheta = 1:nrow(thetas1), .combine = rbind) %dorng% {
  theta1 <- thetas1[itheta,]
  kernel <- get_MH_module2(theta1 = theta1, Sigma_proposal = theta2plugin_variance)
  chain_state <- theta2plugin_mean
  for (iteration in 1:niterations){
    chain_state <- kernel(chain_state)
  }
  chain_state  
}
#
plot(theta2s[,1], theta2s[,2])
hist(theta2s[,1])
hist(theta2s[,2])

# this naive approach is wasteful, as 2000 iterations are probably too many 
# but on the other hand it's hard to know exactly how many iterations are necessary
# See the file unbiasedcut.R  in inst/tutorial 
# for a more sophisticated approach to sampling from the cut distribution
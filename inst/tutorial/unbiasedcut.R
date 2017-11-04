library(bettertogether)
# install package: debiasedmcmc, available at https://github.com/pierrejacob/debiasedmcmc
library(debiasedmcmc)
setmytheme()
rm(list = ls())
set.seed(21)
registerDoParallel(cores = detectCores())
#
# nhpv considered as Y
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
y1 <- matrix(nhpv, nrow = 1)
# Npart is put in the parameters
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
           143, 229, 696, 93)
J <- 13
# posterior is beta in each study, with parameters
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv

###
# For module 2, ncases considered data
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
y2 <- matrix(ncases, nrow = 1)
# Npop considered parameters
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)
Npop_normalized <- log(10**(-3) * Npop)
# Find parameters given Y2
hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
dprior2 <- function(theta2, hyper2){
  return(sum(dnorm(theta2, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE)))
}

nsamples <- 2^8
thetas1 <- foreach(isample = 1:nsamples, .combine = rbind) %dorng% {
  rbeta(J, posterior_phi_alpha, posterior_phi_beta)
}


# get Metropolis-Hastings kernel targeting the posterior of
# theta2 given theta1 and Y2
get_kernels <- function(theta1, Sigma_proposal){
  target <- function(x) plummer_module2_conditional(theta1, matrix(x, nrow = 1), ncases, Npop_normalized) + dprior2(x, hyper2)
  ##
  # Markov kernel of the chain
  single_kernel <- function(chain_state){
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
  # Markov kernel of the coupled chain
  Sigma_chol <- t(chol(Sigma_proposal))
  Sigma_chol_inv <- t(solve(chol(Sigma_proposal)))
  coupled_kernel <- function(chain_state1, chain_state2){
    distance_ <- mean((chain_state1 - chain_state2)^2)
    proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    proposal_pdf1 <- target(proposal1)
    proposal_pdf2 <- target(proposal2)
    current_pdf1 <- target(chain_state1)
    current_pdf2 <- target(chain_state2)
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(target = target, coupled_kernel = coupled_kernel, single_kernel = single_kernel))
}

# compute mean of first posterior
theta1hat <- colMeans(thetas1)
# we are going to tune the MH algorithm conditional on theta1hat
# and cross fingers that this will give a decent tuning for other draws theta1 from the first posterior
dimension <- 2
Sigma_proposal <- diag(c(1,1), dimension, dimension)
rinit <- function() fast_rmvnorm(1, mean = c(0, 0), diag(1, dimension, dimension))[1,]
# first, get a sense of the meeting times
kernels <- get_kernels(theta1hat, Sigma_proposal)
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
# choose values for k and m
k <- as.numeric(floor(quantile(meetingtime, probs = 0.95)))
m <- 5*k

c_chains_continued_ <-  foreach(irep = 1:nsamples) %dorng% {
  continue_coupled_chains(c_chains_[[irep]], kernels$single_kernel, K = m)
}
# then compute mean and variance of plug-in distribution  
mean_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], k = k, K = m)
}
square_estimators <-  foreach(irep = 1:nsamples) %dorng% {
  H_bar(c_chains_continued_[[irep]], h = function(x) x^2, k = k, K = m)
}
est_mean <- rep(0, dimension)
est_var <- rep(0, dimension)
for (component in 1:dimension){
  estimators <- sapply(mean_estimators, function(x) x[component])
  est_mean[component] <- mean(estimators)
  s_estimators <- sapply(square_estimators, function(x) x[component])
  est_var[component] <- mean(s_estimators) - est_mean[component]^2
}
Sigma_proposal <- diag(est_var, dimension, dimension)
# this yields the following variance estimator
print(Sigma_proposal)
# now we can restart with a better tuned MH
rinit <- function() fast_rmvnorm(1, mean = est_mean, Sigma_proposal)[1,]
c_chains_ <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- thetas1[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit)
}
meetingtime <- sapply(c_chains_, function(x) x$meetingtime)
summary(meetingtime)
hist(meetingtime)
# meeting time is now much smaller!
# now we can choose a new value for k and m
k <- as.numeric(floor(quantile(meetingtime, probs = 0.95)))
m <- 10*k
# ##
c_chains_final <-  foreach(irep = 1:nsamples) %dorng% {
  theta1 <- thetas1[irep,]
  kernels <- get_kernels(theta1, Sigma_proposal)
  coupled_chains(kernels$single_kernel, kernels$coupled_kernel, rinit, K = m)
}
## histogram
histogram1 <- histogram_c_chains(c_chains_final, 1, k, m, nclass = 30)
g1 <- plot_histogram(histogram1, with_bar = T) + xlab(expression(theta[2.1])) + ylab("density")
g1
histogram2 <- histogram_c_chains(c_chains_final, 2, k, m, nclass = 30)
g2 <- plot_histogram(histogram2, with_bar = T) + xlab(expression(theta[2.2])) + ylab("density")
g2

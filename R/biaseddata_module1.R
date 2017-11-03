#'@rdname biaseddata_module1
#'@title Biased Data - module 1
#'@description Generate a module object, in the form of a list, that contains all
#' we need from module 1: a function to evaluate the prior, a function to evaluate the log-likelihood,
#' ...
#'@export

# theta1 ~ Normal(theta1_mean_prior, theta1_prior_sd)
# Y_1 ~ Normal(theta1, 1)

biaseddata_module1 <- function(){
  # some hyperparameters
  hyper1 <- list(theta1_mean_prior = 0, theta1_prior_sd = 1)
  # logdensity of prior
  dprior <- function(theta1s, hyper1){
    return(dnorm(theta1s, mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd, log = TRUE))
  }
  # loglikelihood
  loglik <- function(theta1s, y1){
    return(sapply(1:length(theta1s), function(index) sum(dnorm(y1, mean = theta1s[index], sd = 1, log = TRUE))))
  }
  # log density of the posterior
  dposterior <- function(theta1s, y1, hyper1){
    return(dprior(theta1s, hyper1) + loglik(theta1s, y1))
  }
  return(list(dprior = dprior, loglik = loglik, dposterior = dposterior,
              hyper1 = hyper1, theta_dimension = 1, y_dimension = 1))
}

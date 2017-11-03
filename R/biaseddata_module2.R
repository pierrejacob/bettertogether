#'@rdname biaseddata_module2
#'@title Biased Data - module 2
#'@description Generate a module object, in the form of a list, that contains all
#' we need from module 2: a function to evaluate the prior, a function to evaluate the log-likelihood,
#' ...
#'@export

# theta2 ~ Normal(theta2_mean_prior, theta2_prior_sd)
# Y_2    ~ Normal(theta1 + theta2, 1)

biaseddata_module2 <- function(){
  # some hyperparameters
  hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = 1)
  # logdensity of prior
  dprior <- function(theta2s, hyper2){
    return(dnorm(theta2s, mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE))
  }
  # loglikelihood; note that we want to use this function for both the joint posterior,
  # and the conditional posterior of theta 2 given one value of theta 1
  loglik <- function(theta1s, theta2s, y2){
    # in other model, the next line should be replaced, e.g. by length(theta1s) == nrow(theta2s) or something
    if (length(theta1s) == length(theta2s)){
      # for pairs (theta1, theta2)
      return(sapply(1:length(theta2s),
                    function(index) sum(dnorm(y2, mean = theta1s[index] + theta2s[index], sd = 1, log = TRUE))))
    } else {
      # theta1s might be made of only one element
      return(sapply(1:length(theta2s),
                    function(index) sum(dnorm(y2, mean = theta1s + theta2s[index], sd = 1, log = TRUE))))
    }
  }
  # log density of the posterior
  dposterior <- function(theta1s, theta2s, y2, hyper2){
    return(dprior(theta2s, hyper2) + loglik(theta1s, theta2s, y2))
  }
  return(list(dprior = dprior, loglik = loglik, dposterior = dposterior,
              hyper2 = hyper2, theta_dimension = 1, y_dimension = 1))
}

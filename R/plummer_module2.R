#'@rdname plummer_module2
#'@title Plummer Model - module 2
#'@description Generate a module object, in the form of a list, that contains all
#' we need from module 1: a function to evaluate the prior, a function to evaluate the log-likelihood,
#' ...
#'@export

# The second module is characterized by the following distributions.
# Parameter: $\theta = (\theta_1, \theta_2) \sim \mathcal{N}(\mu = 0, \sigma^2 = 10^{3})$, independently.
# Likelihood: $\forall i \in \{1,\ldots, J\} \quad  Y_i \sim \mathcal{P}\text{oisson}(\mu_i)$,
# where $Y_i$ is ncases[i] and
# \[\forall i \in \{1,\ldots, J\} \quad  \mu_i = \exp\left(\theta_1 + \phi_i \theta_2 + D_i \right)\]
# where $D_i = \log (10^{-3} \text{Npop}[i])$.

plummer_module2 <- function(){
  ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
  Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
            26751, 75815, 150302, 354993, 3683043, 507218)
  Npop_normalized <- log(10**(-3) * Npop)
  J <- 13
  # some hyperparameters
  hyper2 <- list(theta2_mean_prior = 0, theta2_prior_sd = sqrt(1000))
  # logdensity of prior
  dprior <- function(theta2s, hyper2){
    return(sapply(1:nrow(theta2s), function(index)
      sum(dnorm(theta2s[index,], mean = hyper2$theta2_mean_prior, sd = hyper2$theta2_prior_sd, log = TRUE))))
  }
  # loglikelihood
  loglik <- function(theta1s, theta2s, ncases, Npop_normalized){
    if (is.null(nrow(theta1s))) theta1s <- matrix(theta1s, nrow = 1)
    if (is.null(nrow(theta2s))) theta2s <- matrix(theta2s, nrow = 1)
    # evals <- rep(0, nrow(theta2s))
    evals <- sapply(1:nrow(theta2s), function(index) {
      mu <- exp(theta2s[index, 1] + theta1s[index,] * theta2s[index, 2] + Npop_normalized)
      return(sum(dpois(x = ncases, lambda = mu, log = TRUE)))
    })
    return(evals)
  }
  # log density of the posterior
  dposterior <- function(theta1s, theta2s, ncases, Npop_normalized, hyper2){
    evals <- dprior(theta2s, hyper2)
    evals <- evals + loglik(theta1s, theta2s, ncases, Npop_normalized)
    return(evals)
  }
  return(list(dprior = dprior, loglik = loglik, dposterior = dposterior,
              hyper2 = hyper2,
              theta_dimension = 2, ncases = ncases, J = J, Npop_normalized = Npop_normalized, Npop = Npop))
}

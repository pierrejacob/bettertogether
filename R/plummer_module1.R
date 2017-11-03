#'@rdname plummer_module1
#'@title Plummer Model - module 1
#'@description Generate a module object, in the form of a list, that contains all
#' we need from module 1: a function to evaluate the prior, a function to evaluate the log-likelihood,
#' ...
#'@export

# The first module is characterized by the following distributions.
#
# Parameter: $\forall i \in \{1,\ldots, J\} \quad \phi_i \sim \mathcal{U}([0,1]) \equiv Beta(1, 1)$.
# Likelihood: $\forall i \in \{1,\ldots, J\} \quad  Z_i \sim \mathcal{B}\text{inom}(N_i, \phi_i)$
#  where $N_i$ is Npart[i] and $Z_i$ is nhpv[i].
# By conjugacy, posterior distribution : $\forall i \in \{1,\ldots, J\} \quad \phi_i \sim Beta(1 + Z_i, 1 + N_i - Z_i)$.

plummer_module1 <- function(){
  nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
  Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
             143, 229, 696, 93)
  J <- 13
  # some hyperparameters
  hyper1 <- list()
  # logdensity of prior
  dprior <- function(theta1s){
    return(sapply(1:nrow(theta1s), function(index) sum(dbeta(theta1s[index,], shape1 = 1, shape2 = 1, log = TRUE))))
  }
  # loglikelihood
  loglik <- function(theta1s, nhpv, Npart){
    if (is.null(nrow(theta1s))) theta1s <- matrix(theta1s, nrow = 1)
    return(sapply(1:nrow(theta1s), function(index) sum(dbinom(nhpv, Npart, prob = theta1s[index,], log = TRUE))))
  }
  # log density of the posterior
  dposterior <- function(theta1s, nhpv, Npart){
    evals <- dprior(theta1s)
    if (any(is.finite(evals))){
      evals[is.finite(evals)] <- evals[is.finite(evals)] + loglik(theta1s[is.finite(evals),], nhpv, Npart)
    }
    return(evals)
  }
  return(list(dprior = dprior, loglik = loglik, dposterior = dposterior,
              theta_dimension = J, nhpv = nhpv, Npart = Npart))
}

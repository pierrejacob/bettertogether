#'@rdname biaseddata_exact
#'@title Biased Data - exact calculations
#'@description Returns a list with convenient functions
#'@export

biaseddata_exact <- function(){
  # calculate posterior moments of theta1 in the first module (returns each partial posterior moments)
  posterior_moments_1 <- function(hyper1, y1){
    n1 <- length(y1)
    precision <- 1:n1 + 1/(hyper1$theta1_prior_sd^2)
    expectation <-  cumsum(y1) / precision
    return(list(expectation = expectation, precision = precision, sd = 1/sqrt(precision)))
  }
  # calculate posterior moments of theta1 and theta2 in the full model (returns each partial posterior
  # during the assimilation of y2, after having assimilated all of y1)
  posterior_moments_full <- function(hyper1, hyper2, y1, y2){
    n1 <- length(y1)
    n2 <- length(y2)
    # prior hyperparameters
    lambda1 <- 1/(hyper1$theta1_prior_sd^2)
    lambda2 <- 1/(hyper2$theta2_prior_sd^2)

    precision1 <- n1 + lambda1 + (1:n2 * lambda2 / (1:n2 + lambda2))
    expectation1 <- (mean(y1) * n1 + cumsum(y2) / (1:n2) * (1:n2 * lambda2) / (1:n2 + lambda2)) / precision1
    precision2 <- lambda2 + 1:n2 * (n1 + lambda1) / (lambda1 + n1 + 1:n2)
    expectation2 <- (cumsum(y2) / (1:n2) * (1:n2 * (n1 + lambda1) / (lambda1 + n1 + 1:n2)) - mean(y1) * (1:n2 * n1 / (lambda1 + n1 + 1:n2))) / precision2
    # covariance
    covar12 <- (cumsum(y2)) / (1:n2 + lambda2) * expectation1 -
      (1:n2 / (1:n2 + lambda2)) * (1/precision1 + (expectation1^2))
    covar12 <- covar12 - (expectation1 * expectation2)
    # theta_sum = theta1 + theta2
    expectation_sum <- expectation1 + expectation2
    variance_sum <- 1/precision1 + 1/precision2 + 2*covar12
    return(list(expectation1 = expectation1, precision1 = precision1, sd1 = 1/sqrt(precision1),
                expectation2 = expectation2, precision2 = precision2, sd2 = 1/sqrt(precision2), covar12 = covar12,
                expectation_sum = expectation_sum, variance_sum = variance_sum))
  }

  # calculate posterior moments of theta1 and theta2, by plugging in the posterior mean of theta1 in the second module
  # (returns each partial posterior during the assimilation of y2, after having assimilated all of y1)
  posterior_moments_plugin <- function(hyper1, hyper2, y1, y2){
    n1 <- length(y1)
    n2 <- length(y2)
    # prior hyperparameters
    lambda1 <- 1/(hyper1$theta1_prior_sd^2)
    lambda2 <- 1/(hyper2$theta2_prior_sd^2)
    #
    precision1 <- n1 + lambda1
    expectation1 <-  sum(y1) / precision1
    #
    expectation2 <- (1:n2) / ((1:n2) + lambda2) * (cumsum(y2) / (1:n2) - expectation1)
    precision2 <- ((1:n2) + lambda2)
    # theta_sum = theta1 + theta2
    expectation_sum <- expectation1 + expectation2
    variance_sum <- 1/precision2
    return(list(expectation1 = expectation1, precision1 = precision1, sd1 = 1/sqrt(precision1),
                expectation2 = expectation2, precision2 = precision2, sd2 = 1/sqrt(precision2), covar12 = 0,
                expectation_sum = expectation_sum, variance_sum = variance_sum))
  }

  # calculate posterior moments of theta1 and theta2 in the cut model
  # (returns each partial posterior during the assimilation of y2, after having assimilated all of y1)
  posterior_moments_cut <- function(hyper1, hyper2, y1, y2){
    n1 <- length(y1)
    n2 <- length(y2)
    # prior hyperparameters
    lambda1 <- 1/(hyper1$theta1_prior_sd^2)
    lambda2 <- 1/(hyper2$theta2_prior_sd^2)
    #
    precision1 <- n1 + lambda1
    expectation1 <-  sum(y1) / precision1
    #
    expectation2 <- (1:n2) / ((1:n2) + lambda2) * (cumsum(y2) / (1:n2) - mean(y1) * n1 / (n1 + lambda1))
    precision2 <- (n1 + lambda1) * ((1:n2) + lambda2) / (n1 + lambda1 + ((1:n2)^2)/((1:n2) + lambda2))
    covar12 <- - (1:n2) / ((1:n2 + lambda2) * (n1 + lambda1))

    # theta_sum = theta1 + theta2
    expectation_sum <- expectation1 + expectation2
    variance_sum <- 1/precision1 + 1/precision2 + 2*covar12
    return(list(expectation1 = expectation1, precision1 = precision1, sd1 = 1/sqrt(precision1),
                expectation2 = expectation2, precision2 = precision2, sd2 = 1/sqrt(precision2), covar12 = covar12,
                expectation_sum = expectation_sum, variance_sum = variance_sum))
  }

  # calculate expected logarithmic scoring rule for predicting Y1,
  # using a candidate representation on theta1 with mean theta1mean and precision theta1precision
  expect_utility_module1 <- function(theta1_star, theta1mean, theta1precision){
    # Y_1 ~ Normal(theta1, 1)
    # so if theta1 ~ Normal(theta1mean, theta1precision)
    # we have the following predictive mean and precision
    predmu <- theta1mean
    predprecision <- 1/(1 + 1/theta1precision)
    if (theta1precision < 1e-100){
      predprecision <- 1
    }
    # for a predictive distribution p(y) = Normal(y; mu, lambda)
    # computes int log p(y) p_star(dy)
    # where p_star(y) = Normal(y, theta1_star, 1)
    #-0.5*log(2*pi) = -0.9189385
    return(-0.9189385 + 0.5 * log(predprecision) - 0.5 * predprecision * ((theta1_star - predmu)^2 + 1))
  }
  # calculate expected logarithmic scoring rule for predicting Y2,
  # using a candidate representation on eta defined as theta1 + theta2
  expect_utility_module2 <- function(theta1_star, theta2_star, eta_mean, eta_precision){
    predmu <- eta_mean
    predprecision <- 1 / (1 + 1/eta_precision)
    if (eta_precision < 1e-100){
      predprecision <- 1
    }
    return(-0.9189385 + 0.5 * log(predprecision) - 0.5 * predprecision *
             (((theta1_star + theta2_star) - predmu)^2 + 1))
    return(exp_util_full)
  }

  return(list(posterior_moments_1 = posterior_moments_1,
              posterior_moments_full = posterior_moments_full,
              posterior_moments_cut = posterior_moments_cut,
              posterior_moments_plugin = posterior_moments_plugin,
              expect_utility_module1 = expect_utility_module1,
              expect_utility_module2 = expect_utility_module2))
}

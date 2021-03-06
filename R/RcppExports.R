# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rmvnorm <- function(nsamples, mean, covariance) {
    .Call('_bettertogether_rmvnorm', PACKAGE = 'bettertogether', nsamples, mean, covariance)
}

dmvnorm <- function(x, mean, covariance) {
    .Call('_bettertogether_dmvnorm', PACKAGE = 'bettertogether', x, mean, covariance)
}

plummer_full_posterior_ <- function(theta1s, theta2s, nhpv, Npart, ncases, Npop_normalized, theta2_mean_prior, theta2_sd_var) {
    .Call('_bettertogether_plummer_full_posterior_', PACKAGE = 'bettertogether', theta1s, theta2s, nhpv, Npart, ncases, Npop_normalized, theta2_mean_prior, theta2_sd_var)
}

plummer_module2_conditional_ <- function(theta1, theta2s, ncases, Npop_normalized) {
    .Call('_bettertogether_plummer_module2_conditional_', PACKAGE = 'bettertogether', theta1, theta2s, ncases, Npop_normalized)
}

plummer_module2_loglikelihood_ <- function(thetas1, theta2s, ncases, Npop_normalized) {
    .Call('_bettertogether_plummer_module2_loglikelihood_', PACKAGE = 'bettertogether', thetas1, theta2s, ncases, Npop_normalized)
}

systematic_resampling_ <- function(weights) {
    .Call('_bettertogether_systematic_resampling_', PACKAGE = 'bettertogether', weights)
}

systematic_resampling_n_ <- function(weights, ndraws) {
    .Call('_bettertogether_systematic_resampling_n_', PACKAGE = 'bettertogether', weights, ndraws)
}

systematic_resampling_n_u_ <- function(weights, ndraws, u) {
    .Call('_bettertogether_systematic_resampling_n_u_', PACKAGE = 'bettertogether', weights, ndraws, u)
}

wmean_ <- function(x, unnormalized_w) {
    .Call('_bettertogether_wmean_', PACKAGE = 'bettertogether', x, unnormalized_w)
}

wcovariance_ <- function(x, unnormalized_w, xbar) {
    .Call('_bettertogether_wcovariance_', PACKAGE = 'bettertogether', x, unnormalized_w, xbar)
}


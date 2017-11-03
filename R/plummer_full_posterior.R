#'@rdname plummer_full_posterior
#'@title Plummer Model - full posterior in c++
#'@export
plummer_full_posterior <- function(thetas){
  nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
  Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
             143, 229, 696, 93)
  ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
  Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
            26751, 75815, 150302, 354993, 3683043, 507218)
  Npop_normalized <- log(10**(-3) * Npop)
  theta2_mean_prior <- 0
  theta2_sd_var <- 1000
  return(plummer_full_posterior_(thetas[,1:13], thetas[,14:15], nhpv, Npart, ncases, Npop_normalized,
                                 theta2_mean_prior, theta2_sd_var))
}

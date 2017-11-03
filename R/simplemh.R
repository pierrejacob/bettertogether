#' #'@rdname simpleMH
#' #'@title Simple Metropolis Hastings
#' #'@description Could be useful.
#' #'@export
#' simpleMH <- function(niterations, nparameters, logtarget,
#'                proposalcovariance = diag(0.001, nparameters), initvalue = rep(0, nparameters),
#'                verbose = TRUE){
#'   chain <- matrix(nrow = niterations, ncol = nparameters)
#'   currenttheta <- initvalue
#'   chain[1,] <- currenttheta
#'   currentlogtarget <- logtarget(chain[1,])
#'   naccepts <- 0
#'   for (iteration in 2:niterations){
#'     if ((iteration %% 1000 == 1 && verbose)){
#'       cat("iteration ", iteration, "/", niterations, "\n")
#'       cat("Acceptance rate so far:", 100 * naccepts / iteration, "%\n")
#'     }
#'     proposaltheta <- normalmodule$rmvnorm(1, mean = currenttheta, sigma = proposalcovariance)
#'     proposallogtarget <- logtarget(proposaltheta)
#'     if (log(runif(1)) < (proposallogtarget - currentlogtarget)){
#'       naccepts <- naccepts + 1
#'       currenttheta <- proposaltheta
#'       currentlogtarget <- proposallogtarget
#'     }
#'     chain[iteration,] <- currenttheta
#'   }
#'   if (verbose) cat("Acceptance rate of M-H:", 100 * naccepts / niterations, "%\n")
#'   return(chain)
#' }

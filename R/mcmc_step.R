#' #'@export
#' mcmc_step <- function(thetas, logtargetvalues, ys, target, param_algo){
#'   nmoves <- param_algo$nmoves
#'   if (is.null(nmoves)){
#'     nmoves <- 1
#'   }
#'   proposal <- param_algo$proposal
#'   nthetas <- nrow(thetas)
#'   #
#'   acceptrate <- 0
#'   if (nmoves > 0){
#'     for (imove in 1:nmoves){
#'       thetas_prop <- proposal$r(thetas, proposal$param_prop)
#'       prior_proposed <- target$dprior(thetas_prop, target$parameters)
#'       prop_proposed <- proposal$d(thetas_prop, proposal$param_prop)
#'       prop_current <- proposal$d(thetas, proposal$param_prop)
#'       #
#'       logratios <- prior_proposed + (prop_current - prop_proposed)
#'       # only compute loglikelihood if prior is not -Inf
#'       loglikelihood_proposed <- rep(-Inf, nthetas)
#'       potential_indices <- which(!is.infinite(logratios))
#'       # for (ipot in potential_indices){
#'       loglikelihood_proposed[potential_indices] <- target$fullloglikelihood(matrix(thetas_prop[potential_indices,], ncol = target$thetadim), ys)
#'       # }
#'       logratios <- logratios + loglikelihood_proposed - logtargetvalues
#'       accepted <- (log(runif(nthetas)) < logratios)
#'       #
#'       thetas[accepted] <- thetas_prop[accepted]
#'       logtargetvalues[accepted] <- loglikelihood_proposed[accepted] + prior_proposed[accepted]
#'       acceptrate <- acceptrate + mean(accepted)
#'     }
#'     acceptrate <- acceptrate/nmoves
#'   } else {
#'     print("no rejuvenation move")
#'   }
#'   return(list(acceptrate = acceptrate, thetas = thetas, logtargetvalues = logtargetvalues))
#' }

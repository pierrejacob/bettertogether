# SMC sampler, recursively assimilating all data
#'@export
smcsampler <- function(y, target, param_algo, thetas = NULL, normw = NULL){
  nthetas <- param_algo$nthetas
  param_algo$original_proposal <- param_algo$proposal
  if (is.null(thetas) && is.null(normw)){
    print("will sample from the prior at the initial step")
    # sample from prior
    thetas <- target$rprior(nthetas, target$parameters)
    logtargetdensities <- target$dprior(thetas, target$parameters)
    normw <- rep(1/nthetas, nthetas)
    logw <- rep(0, nthetas)
  } else {
    if (nrow(thetas) != param_algo$nthetas){
      print("mismatch between desired number of particles and number of particles given as prior samples")
    }
    # assumes that thetas, weighted with normw, already comes from
    # the distribution with density target$dprior
    thetas <- thetas
    logtargetdensities <- target$dprior(thetas, target$parameters)
    normw <- normw
    logw <- log(normw)
  }
  #
  thetas_history <- list()
  thetas_history[[1]] <- thetas
  normw_history <- list()
  normw_history[[1]] <- normw
  #
  logevidence <- rep(0, nrow(y))
  for (idata in 1:nrow(y)){
    results <- assimilate_one(y, idata, target, param_algo, thetas, logtargetdensities, logw, normw)
    param_algo <- results$param_algo
    thetas <- results$thetas
    normw <- results$normw
    logw <- results$logw
    logtargetdensities <- results$logtargetdensities
    logevidence[idata] <- results$logcst
    thetas_history[[idata+1]] <- thetas
    normw_history[[idata+1]] <- normw
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, logevidence = logevidence,
              logtargetdensities = logtargetdensities, param_algo = param_algo))
}
#

# move step
#'@export
move_step <- function(thetas, logtargets, target_dlog, param_algo){
  nmoves <- param_algo$nmoves
  if (is.null(nmoves)){
    nmoves <- 1
  }
  proposal <- param_algo$proposal
  nthetas <- nrow(thetas)
  #
  acceptrate <- 0
  if (nmoves > 0){
    for (imove in 1:nmoves){
      # sample proposals
      param_algo <- update_proposal(param_algo, thetas, rep(1/nthetas, nthetas))
      proposal <- param_algo$proposal
      # proposal$param_prop <- proposal$param_update(thetas, rep(1/nthetas, nthetas))
      thetas_prop <- proposal$r(thetas, proposal$param_prop)
      # evaluate proposals' log-density
      prop_proposed <- proposal$d(thetas_prop, proposal$param_prop)
      prop_current <- proposal$d(thetas, proposal$param_prop)
      # evaluate proposals' target log-density
      targets_prop <- target_dlog(thetas_prop)
      # acceptance log probability
      logratios <- targets_prop - logtargets + (prop_current - prop_proposed)
      # accept or not
      accepted <- (log(runif(nthetas)) < logratios)
      # perform appropriate replacement
      thetas[accepted,] <- thetas_prop[accepted,]
      logtargets[accepted] <- targets_prop[accepted]
      acceptrate <- acceptrate + mean(accepted)
    }
    acceptrate <- acceptrate/nmoves
  } else {
    print("no move")
  }
  return(list(acceptrate = acceptrate, thetas = thetas, logtargets = logtargets, param_algo = param_algo))
}

# assimilate one data point, using tempering
#'@export
assimilate_one <- function(y, idata, target, param_algo, thetas, logtargetdensities, logw, normw){
  nthetas <- param_algo$nthetas
  minimum_diversity <- param_algo$minimum_diversity
  current_gamma <- 0
  logcst <- 0
  while (current_gamma < 1){
    logw_incremental <- target$conditionallikelihood(thetas, y, idata, target$parameters)
    one_uniform <- runif(1)
    diversity_given_gamma <- function(gamma){
      logw_ <- logw + (gamma - current_gamma) * logw_incremental
      maxlogw <- max(logw_)
      w <- exp(logw_ - maxlogw)
      normw <- w / sum(w)
      a <- systematic_resampling_given_u(normw, nthetas, one_uniform)
      th <- thetas[a,1]
      return(length(unique(th))/nthetas)
    }
    # try gamma = 1 first
    if (diversity_given_gamma(1) > minimum_diversity){
      gamma <- 1
    } else {
      if (diversity_given_gamma(current_gamma) < minimum_diversity){
        gamma <- current_gamma
      } else {
        optimal_gamma <- optim(par = current_gamma, fn = function(x) (diversity_given_gamma(x) - minimum_diversity)^2,
                               method = "L-BFGS-B", lower = current_gamma, upper = 1)
        gamma <- optimal_gamma$par
        current_diversity <- diversity_given_gamma(gamma)
        if (current_diversity < (0.9*minimum_diversity) || gamma == current_gamma){
          # then something went really wrong in the above optimization program, so we turn to a
          # manual way of binary searching gamma
          maxtrials <- 100
          gamma_increment <- 1 - current_gamma
          current_diversity <- diversity_given_gamma(gamma + gamma_increment)
          trial <- 1
          while ((current_diversity < minimum_diversity) && (trial < maxtrials)){
            gamma_increment <- gamma_increment / 2
            gamma <- current_gamma + gamma_increment
            current_diversity <- diversity_given_gamma(gamma)
            trial <- trial + 1
          }
          gamma <- current_gamma + gamma_increment
          cat("using gamma = ", gamma, " with corresponding diversity of ", diversity_given_gamma(gamma), "\n")
        }
      }
    }
    # now we've found our gamma
    logw_incremental_gamma <- (gamma - current_gamma) * logw_incremental
    logtargetdensities <- logtargetdensities + logw_incremental_gamma
    current_gamma <- gamma
    # compute increment to the normalizing constant
    maxlogw <- max(logw_incremental_gamma)
    w <- exp(logw_incremental_gamma - maxlogw)
    logcst <- logcst + log(sum(normw * w)) + maxlogw
    # normalize weights
    logw <- logw + logw_incremental_gamma
    w <- exp(logw - max(logw))
    normw <- w / sum(w)
    ##
    cat("Step", idata, ", gamma = ", gamma, ", ESS = ", ESSfunction(normw), "\n")
    if (gamma<1){
      # we need to resample and move
      # resampling step
      ancestors <- systematic_resampling(normw)
      thetas <- matrix(thetas[ancestors,], ncol = target$thetadim)
      logtargetdensities <- logtargetdensities[ancestors]
      logw <- rep(0, nthetas)
      normw <- rep(1/nthetas, nthetas)
      # move step
      dtarget <- function(thetas){
        eval <- target$dprior(thetas, target$parameters)
        # don't evaluate full log likelihood if the prior precludes the value
        which.finiteprior <- which(is.finite(eval))
        if ((length(which.finiteprior)) > 0){
          if (idata > 1){
            eval[which.finiteprior] <- eval[which.finiteprior] + target$fullloglikelihood(thetas[which.finiteprior,,drop=FALSE], y[1:(idata-1),,drop=FALSE], parameters = target$parameters)
          }
          eval[which.finiteprior] <- eval[which.finiteprior] + gamma * target$conditionallikelihood(thetas[which.finiteprior,,drop=FALSE], y, idata, parameters = target$parameters)
        }
        return(eval)
      }
      # moves with independent proposal
      move_results <- move_step(thetas, logtargets = logtargetdensities, target_dlog =  dtarget, param_algo = param_algo)
      thetas <- move_results$thetas
      logtargetdensities <- move_results$logtargets
      param_algo <- move_results$param_algo
      cat("Acceptance rate (independent proposal): ", 100*move_results$acceptrate, "%\n")
    }
  }
  return(list(thetas = thetas, normw = normw, logw = logw, logtargetdensities = logtargetdensities, logcst = logcst, param_algo = param_algo))
}


#'@export
update_proposal <- function(param_algo, thetas, weights){
  # reverse back to original proposal, in case it's been replaced by something else
  param_algo$proposal <- param_algo$original_proposal
  param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
  a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
  if (inherits(a, "try-error")){
    print("error in fitting MCMC proposal, trying Rmixmod 5 times..")
    param_algo$proposal <- mixture_rmixmod()
    param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
    a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
    nattempts <- 5
    attempt <- 0
    while ((inherits(a, "try-error")) && (attempt < nattempts)){
      attempt <- attempt + 1
      param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
      a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
    }
    if (inherits(a, "try-error")){
      print("error in fitting MCMC proposal, switching to independent Gaussian")
      param_algo$proposal <- independent_proposal()
      param_algo$proposal$param_prop <- param_algo$proposal$param_update(thetas, weights)
      a <- try(param_algo$proposal$r(thetas[1:2,,drop=F], param_prop = param_algo$proposal$param_prop))
      if (inherits(a, "try-error")){
        print("independent Gaussian proposal failed")
        tmpfilename <- "~/tmpfitfail.RData"
        save(thetas, weights, param_algo, file = tmpfilename)
        cat("dumping thetas, weights, param_algo in ", tmpfilename,
            "; try executing 'load('", tmpfilename, "'); param_algo$proposal$param_update(thetas, weights)'\n")
      }
    }
  }
  return(param_algo)
}


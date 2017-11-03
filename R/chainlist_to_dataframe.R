#'@rdname chainlist_to_dataframe
#'@title Convert list of chains to data.frame
#'@description to handle the output of metropolishastings
#'@export
chainlist_to_dataframe <- function(chains_list){
  nchains <- length(chains_list)
  niterations <- nrow(chains_list[[1]])
  chaindf <- foreach (i = 1:nchains, .combine = rbind) %do% {
    data.frame(ichain = rep(i, niterations), iteration = 1:niterations, X = chains_list[[i]])
  }
  return(chaindf)
}

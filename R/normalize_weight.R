#'@rdname normalize_weight
#'@title normalize_weight
#'@description normalize log weights
#'@param w vector of log weights
#'@export
normalize_weight <- function(w){
  w <- exp(w - max(w))
  w <- w / sum(w)
  return(w)
}
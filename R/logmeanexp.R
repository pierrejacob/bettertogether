#'@export
logmeanexp <- function(logw){
  ml <- max(logw)
  l <- logw - ml
  return(ml + log(mean(exp(l))))
}

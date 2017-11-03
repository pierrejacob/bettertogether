#'@rdname plummer_module2_conditional
#'@title Plummer Model - module 2 posterior conditional on theta1 in c++
#'@export
plummer_module2_conditional <- function(theta1, theta2s, ncases, Npop_normalized){
  return(plummer_module2_conditional_(theta1, theta2s, ncases, Npop_normalized))
}

# following function is for when theta1s and theta2s have same number of rows
#'@export
plummer_module2_loglikelihood <- function(theta1s, theta2s, ncases, Npop_normalized){
  return(plummer_module2_loglikelihood_(theta1s, theta2s, ncases, Npop_normalized))
}

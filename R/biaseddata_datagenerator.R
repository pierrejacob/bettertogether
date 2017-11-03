#'@rdname biaseddata_datagenerator
#'@title Biased Data - generator
#'@description Generate synthetic dataset from the biased data model
#'@export

# Y_1 ~ Normal(theta1_star, 1)
# Y_2 ~ Normal(theta1_star + theta2_star, 1)


biaseddata_datagenerator <- function(n1, n2, theta1_star, theta2_star){
  y1 <- rnorm(n1, mean=theta1_star, sd = 1)
  y2 <- rnorm(n2, mean=theta1_star + theta2_star, sd = 1)
  return(list(y1 = y1, y2 = y2))
}

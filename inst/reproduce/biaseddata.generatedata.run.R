library(bettertogether)
set.seed(16)

## data-generating process
n1 <- 100
n2 <- 1000
theta1_star <- 0
theta2_star <- 1
synthetic_dataset <- biaseddata_datagenerator(n1, n2, theta1_star, theta2_star)
save(n1, n2, theta1_star, theta2_star, synthetic_dataset, file = "biaseddata_dataset.RData")

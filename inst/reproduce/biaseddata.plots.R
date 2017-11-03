##
library(bettertogether)
library(xtable)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
setmytheme()
set.seed(16)

## data
load("biaseddata_dataset.RData")

y1 <- matrix(synthetic_dataset$y1, ncol = 1)
y2 <- matrix(synthetic_dataset$y2, ncol = 1)

K <- 10

#
rep <- 5
## Module 1
load(paste0("biaseddata_module1.N1024.K", K, ".rep5.RData"))
#
theta.df1 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results1[[i]]$thetas,
  weight = results1[[i]]$normw,
  irep = rep(i, nrow(results1[[i]]$thetas)))
}

load(paste0("biaseddata_module2givenmodule1.N1024.K", K, ".rep5.RData"))
rep <- length(results2given1)
theta.df2given1 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results2given1[[i]]$thetas[,1],
             theta2 = results2given1[[i]]$thetas[,2],
             weight = results2given1[[i]]$normw,
             irep = rep(i, nrow(results2given1[[i]]$thetas)))
}
load(paste0("biaseddata_module1givenmodule2.N1024.K", K, ".rep5.RData"))
theta.df1given2 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results1given2[[i]]$thetas[,1],
             theta2 = results1given2[[i]]$thetas[,2],
             weight = results1given2[[i]]$normw,
             irep = rep(i, nrow(results1given2[[i]]$thetas)))
}
### Plug-in approach
load(paste0("biaseddata_module2givenmodule1plugin.N1024.K", K, ".rep5.RData"))
theta.df2given1plugin <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta2 = results2given1plugin[[i]]$thetas[,1],
             weight = results2given1plugin[[i]]$normw,
             irep = rep(i, nrow(results2given1plugin[[i]]$thetas)),
             s = results2given1plugin[[i]]$theta1hat + results2given1plugin[[i]]$thetas[,1])
}

## Cut approach
load(paste0("biaseddata_module2givenmodule1cut.N1024.K10.M1024.rep5.RData"))
N1 <- 1024
theta.df2given1cut <- data.frame()
for (irep in 1:rep){
  theta.df2given1cut_ <- foreach(i = 1:N1, .combine = rbind) %dopar% {
    Nreduced <- 10
    theta2 <-  results2given1cut[[irep]][[i]]$thetas2
    normw2 <- results2given1cut[[irep]][[i]]$normw2
    theta2 <- theta2[systematic_resampling_n(normw2, Nreduced),]
    theta1s <- matrix(0, nrow = Nreduced, ncol = 1)
    for (j in 1:Nreduced){
      theta1s[j,] <- results2given1cut[[irep]][[i]]$theta1
    }
    data.frame(theta1 = theta1s, theta2 = theta2)
  }
  theta.df2given1cut_$irep <- irep
  theta.df2given1cut <- rbind(theta.df2given1cut, theta.df2given1cut_)
}


load(paste0("biaseddata_module2alone.N1024.K", K, ".rep5.RData"))
theta.df2alone <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results2alone[[i]]$thetas[,1],
             theta2 = results2alone[[i]]$thetas[,2],
             weight = results2alone[[i]]$normw,
             irep = rep(i, nrow(results2alone[[i]]$thetas)))
}

#### estimation of theta1
g <- ggplot(theta.df1, aes(x = theta1, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "module 1  "))
# g <- g + geom_density(data = theta.df1given2, aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2given1, aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + xlab(expression(theta[1]))
g <- g + geom_vline(xintercept = 0) + stat_function(fun = dnorm, aes(colour = "prior  "))
g <- g + scale_colour_manual(name = "", values = c("orange", "black", "skyblue", "brown"))
g <- g + xlim(-2,2)
g
ggsave(plot = g, filename = "biaseddata.theta1.pdf", width = 10, height = 7)

g <- ggplot(theta.df1given2, aes(x = theta2, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + geom_density(data = theta.df2given1cut, aes(y = ..density.., weight = NULL, colour = "cut  "))
g <- g + geom_density(data = theta.df2given1plugin, aes(y = ..density.., colour = "plugin  "))
g <- g + xlab(expression(theta[2]))
g <- g + geom_vline(xintercept = 1) + stat_function(fun = function(x) dnorm(x, 0, 0.1), aes(colour = "prior  "))
g <- g + scale_colour_manual(name = "", values = c("black", "orange", "skyblue", "red", "brown"))
g
ggsave(plot = g, filename = "biaseddata.theta2.pdf", width = 10, height = 7)

#### overlay multiple posteriors using 2d contours
g <- ggplot(theta.df2given1cut, aes(x = theta1, y = theta2, group = irep))
g <- g + stat_density_2d(geom = "polygon", alpha = 0.2, aes(fill = "cut  "))
g <- g + stat_density_2d(data = theta.df2alone, alpha = 0.2, geom = "polygon", aes(fill = "module 2  ", weight = weight))
g <- g + stat_density_2d(data = theta.df2given1, alpha = 0.2, geom = "polygon", aes(fill = "full posterior  ", weight = weight))
g <- g + xlab(expression(theta[1])) + ylab(expression(theta[2]))
g <- g + geom_vline(xintercept = 0) + geom_hline(yintercept = 1)
g <- g + geom_abline(intercept = 1, slope = -1, linetype = 2)
g <- g + scale_x_continuous(breaks = 0:5/5)  + scale_y_continuous(breaks = 0:5/5)
g <- g + scale_fill_manual(name = "", values = c("black", "orange", "skyblue", "red", "brown"))
g
ggsave(plot = g, filename = "biaseddata.thetas.pdf", width = 10, height = 7)

theta.df2given1cut <- theta.df2given1cut %>% mutate(s = theta1 + theta2)
theta.df2alone <- theta.df2alone %>% mutate(s = theta1 + theta2)
theta.df2given1 <- theta.df2given1 %>% mutate(s = theta1 + theta2)

g <- ggplot(theta.df2given1, aes(x = s, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + geom_density(data = theta.df2given1plugin, aes(y = ..density.., colour = "plugin  "))
g <- g + geom_density(data = theta.df2given1cut, aes(y = ..density.., weight = NULL, colour = "cut  "))
g <- g + xlab(expression(theta[1]+theta[2]))
g <- g + geom_vline(xintercept = 1)
g <- g + scale_colour_manual(name = "", values = c("black", "orange", "skyblue", "red", "brown"))
g
ggsave(plot = g, filename = "biaseddata.thetasum.pdf", width = 10, height = 7)


## What about the prior predictive performance?
module1 <- biaseddata_module1()
hyper1 <- module1$hyper1
## Module 2
module2 <- biaseddata_module2()
hyper2 <- module2$hyper2

nthetas <- 1024
logevidence1prior <- as.numeric(foreach(i = 1:rep, .combine = c) %dorng% {
  theta_prior <- matrix(rnorm(nthetas, mean = hyper1$theta1_mean_prior, sd = hyper1$theta1_prior_sd), ncol = 1)
  prior_logscore <- rep(0, n1)
  for (i in 1:n1){
    evals <- module1$loglik(theta_prior, y1 = y1[i])
    prior_logscore[i] <- max(evals) + log(mean(exp(evals - max(evals))))
  }
  sum(prior_logscore)
})

logevidence2alone <- as.numeric(foreach(i = 1:rep, .combine = c) %dorng% {
  theta1 <- results2alone[[i]]$thetas[,1]
  normw1 <- results2alone[[i]]$normw
  prior_logscore <- rep(0, n1)
  for (i in 1:n1){
    evals <- module1$loglik(theta1, y1 = y1[i])
    prior_logscore[i] <- max(evals) + log(sum(normw1 * exp(evals - max(evals))))
  }
  sum(prior_logscore)
})

logevidences1 <- c()
logevidences1given2 <- c()
for (i in 1:rep){
  logevidences1 <- c(logevidences1, sum(results1[[i]]$logevidence))
  logevidences1given2 <- c(logevidences1given2, sum(results1given2[[i]]$logevidence))
}

## prepare table
f <- function(logevid, ro = 1){
  m <- round(mean(logevid), ro)
  mmin <- round(min(logevid), ro)
  mmax <- round(max(logevid), ro)
  print_evid <- paste0(m, " [", mmin, ", ", mmax, "]")
  return(print_evid)
}

tab <- data.frame(module1 = f(logevidences1),
                  prior = f(logevidence1prior),
                  full = f(logevidences1given2),
                  module2 = f(logevidence2alone))
tab <- t(tab)
# tab <- rbind(c("predictive score on $Y_1$"), tab)
rownames(tab) <- c("module 1", "prior", "full model", "module 2")
colnames(tab) <- "predictive score on $Y_1$"
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T)
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T,
      file = "biaseddata.predictY1.tex")


## prepare for table
N1 = length(results2given1cut[[1]])
#

logevidences2given1 <- c()
logevidences2alone <- c()
logevidences2plugin <- c()
logevidences2cut <- c()
for (i in 1:5){
  logevidences2given1 <- c(logevidences2given1, sum(results2given1[[i]]$logevidence))
  logevidences2alone <- c(logevidences2alone, sum(results2alone[[i]]$logevidence))
  logevidences2plugin <- c(logevidences2plugin, results2given1plugin[[i]]$sumlogevidence)
  logevidences2cut <- c(logevidences2cut, scores2given1cut[[i]])
  # logevid <- rep(0, N1)
  # for (j in 1:N1){
  #   logevid[j] <- results2given1cut[[i]][[j]]$sumlogevidence
  # }
  # maxlogevid <- max(logevid)
  # logevidences2cut <- c(logevidences2cut, maxlogevid + log(mean(exp(logevid - maxlogevid))))
}

logevidences2cut

#
f <- function(logevid, ro = 1){
  m <- round(mean(logevid), ro)
  mmin <- round(min(logevid), ro)
  mmax <- round(max(logevid), ro)
  print_evid <- paste0(m, " [", mmin, ", ", mmax, "]")
  return(print_evid)
}

tab <- data.frame(module2 = f(logevidences2alone),
                  full = f(logevidences2given1),
                  cut = f(logevidences2cut),
                  plugin = f(logevidences2plugin))
tab <- t(tab)
xtable(tab)

rownames(tab) <- c("module 2", "full model", "cut approach", "plug-in approach")
colnames(tab) <- "predictive score on $Y_2$"
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T)
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T,
      file = "biaseddata.predictY2.tex")

library(bettertogether)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
setmytheme()
registerDoParallel(cores = detectCores())
set.seed(16)
#
# nhpv considered as Y
nhpv <- c(7, 6, 10, 10, 1, 1, 10, 4, 35, 0, 10, 8, 4)
y1 <- matrix(nhpv, nrow = 1)
# Npart is put in the parameters
Npart <- c(111, 71, 162, 188, 145, 215, 166, 37, 173,
           143, 229, 696, 93)
J <- 13
# For module 2, ncases considered data
ncases <- c(16, 215, 362, 97, 76, 62, 710, 56, 133,28, 62, 413, 194)
y2 <- matrix(ncases, nrow = 1)
# Npop considered parameters
Npop <- c(26983, 250930, 829348, 157775, 150467, 352445, 553066,
          26751, 75815, 150302, 354993, 3683043, 507218)
Npop_normalized <- log(10**(-3) * Npop)
#
posterior_phi_alpha <- 1 + nhpv
posterior_phi_beta <- 1 + Npart - nhpv


### Module 1 only
load("epidemiology_module1altogether.N1024.K100.rep5.RData")

rep <- length(results1)
theta.df1 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results1[[i]]$thetas,
             weight = results1[[i]]$normw,
             irep = rep(i, nrow(results1[[i]]$thetas)))
}
theta.df1 %>% head

### Full
load("epidemiology_module2givenmodule1.N1024.K100.rep5.RData")
results2given1[[1]]$thetas %>% head

theta.df2given1 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results2given1[[i]]$thetas[,1:13],
             theta2 = results2given1[[i]]$thetas[,14:15],
             weight = results2given1[[i]]$normw,
             irep = rep(i, nrow(results2given1[[i]]$thetas)))
}
theta.df2given1 %>% head

load("epidemiology_module1givenmodule2.N1024.K100.rep5.RData")
theta.df1given2 <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results1given2[[i]]$thetas[,1:13],
             theta2 = results1given2[[i]]$thetas[,14:15],
             weight = results1given2[[i]]$normw,
             irep = rep(i, nrow(results1given2[[i]]$thetas)))
}
### Plug-in approach
load("epidemiology_module2givenmodule1plugin.N1024.K100.rep5.RData")
theta.df2given1plugin <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta2 = results2given1plugin[[i]]$thetas,
             weight = results2given1plugin[[i]]$normw,
             irep = rep(i, nrow(results2given1plugin[[i]]$thetas)))
}
theta.df2given1plugin %>% head

## Cut approach
load("epidemiology_module2givenmodule1cut.N1024.K100.M1024.rep5.RData")
theta.df2given1cut <- data.frame()
for (irep in 1:rep){
  theta.df2given1cut_ <- foreach(i = 1:1024, .combine = rbind) %dopar% {
    Nreduced <- 10
    theta2 <-  results2given1cut[[irep]][[i]]$thetas2
    normw2 <- results2given1cut[[irep]][[i]]$normw2
    theta2 <- theta2[systematic_resampling_n(normw2, Nreduced),]
    theta1s <- matrix(0, nrow = Nreduced, ncol = 13)
    for (j in 1:Nreduced){
      theta1s[j,] <- results2given1cut[[irep]][[i]]$theta1
    }
    data.frame(theta1 = theta1s, theta2 = theta2)
  }
  theta.df2given1cut_$irep <- irep
  theta.df2given1cut <- rbind(theta.df2given1cut, theta.df2given1cut_)
}

theta.df2given1cut %>% head

load("epidemiology_module2alone.N1024.K100.rep5.RData")

rep <- length(results2alone)
theta.df2alone <- foreach(i = 1:rep, .combine = rbind) %dopar% {
  data.frame(theta1 = results2alone[[i]]$thetas[,1:13],
             theta2 = results2alone[[i]]$thetas[,14:15],
             weight = results2alone[[i]]$normw,
             irep = rep(i, nrow(results2alone[[i]]$thetas)))
}

theta.df1given2 %>% head

g <- ggplot(theta.df1, aes(x = theta1.1, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "module 1  "))
g <- g + geom_density(data = theta.df1given2, aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + xlab(expression(theta[1.1]))
g <- g + scale_colour_manual(name = "", values = c("orange", "black", "skyblue", "brown"))
g <- g + xlim(0,1)
g
ggsave(plot = g, filename = "epidemiology.theta1.1.pdf", width = 10, height = 7)

g <- ggplot(theta.df1, aes(x = theta1.9, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "module 1  "))
g <- g + geom_density(data = theta.df1given2, aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + xlab(expression(theta[1.9]))
g <- g + scale_colour_manual(name = "", values = c("orange", "black", "skyblue", "brown"))
g <- g + xlim(0,1)
g
ggsave(plot = g, filename = "epidemiology.theta1.9.pdf", width = 10, height = 7)


#### overlay multiple posteriors using 2d contours
g <- ggplot(theta.df2given1cut, aes(x = theta2.1, y = theta2.2, group = irep))
g <- g + stat_density_2d(geom = "polygon", alpha = 0.2, aes(fill = "cut  "))
g <- g + stat_density_2d(data = theta.df2alone, alpha = 0.2, geom = "polygon", aes(fill = "module 2  ", weight = weight))
g <- g + stat_density_2d(data = theta.df2given1, alpha = 0.2, geom = "polygon", aes(fill = "full posterior  ", weight = weight))
g <- g + xlab(expression(theta[2.1])) + ylab(expression(theta[2.2]))
g <- g + scale_fill_manual(name = "", values = c("black", "orange", "skyblue", "red", "brown"))
g <- g + xlim(-4,3) + ylim(-10,35)
g
ggsave(plot = g, filename = "epidemiology.theta2.1.pdf", width = 10, height = 7)

g <- ggplot(theta.df1given2, aes(x = theta2.2, weight = weight, group = irep)) + geom_density(aes(y = ..density.., colour = "full posterior  "))
g <- g + geom_density(data = theta.df2alone, aes(y = ..density.., colour = "module 2  "))
g <- g + geom_density(data = theta.df2given1cut, aes(y = ..density.., weight = NULL, colour = "cut  "))
g <- g + geom_density(data = theta.df2given1plugin, aes(y = ..density.., colour = "plugin  "))
g <- g + xlab(expression(theta[2.2]))
g <- g + stat_function(fun = function(x) dnorm(x, 0, sqrt(1000)), aes(colour = "prior  "))
g <- g + scale_colour_manual(name = "", values = c("black", "orange", "skyblue", "red", "brown"))
g <- g + xlim(-10, 35)
g
ggsave(plot = g, filename = "epidemiology.theta2.2.pdf", width = 10, height = 7)

### tables
### Predict Y1
# Module 1 only
load("epidemiology_module1altogether.N1024.K100.rep5.RData")
names(results1[[1]])
logevidences1 <- c()
for (i in 1:length(results1)){
  logevidences1 <- c(logevidences1, results1[[i]]$logevidence)
}

logevidences1
# in fact we can get the exact evidence here, by conjugacy
-sum(log(Npart+1))

### Module 1 given module 2
load("epidemiology_module1givenmodule2.N1024.K100.rep5.RData")
logevidences1given2 <- c()
for (i in 1:length(results1given2)){
  logevidences1given2 <- c(logevidences1given2, results1given2[[i]]$logevidence)
}
logevidences1given2
#
##  Prior in module 1
## What about the prior predictive performance?
sumlogevidprior <- 0
for (index in 1:J){
  probs <- rbeta(1e5, 1, 1)
  evid <- log(mean(dbinom(y1[index], size = Npart[index], prob = probs)))
  sumlogevidprior <- sumlogevidprior + evid
}
print(sumlogevidprior)

load("epidemiology_module2alone.N1024.K100.rep5.RData")
logevidences2alone <- c()
for (i in 1:length(results2alone)){
  theta1 <- results2alone[[i]]$thetas[,1:13]
  normw1 <- results2alone[[i]]$normw
  prior_logscore <- rep(0, J)
  for (j in 1:J){
    evals <- dbinom(y1[j], size = Npart[j], prob = theta1[,j], log = TRUE)
    prior_logscore[j] <- max(evals) + log(sum(normw1 * exp(evals - max(evals))))
  }
  logevidences2alone <- c(logevidences2alone, sum(prior_logscore))
}
logevidences2alone
#### TO SUMMARIZE
### module 1 predictive score
-sum(log(Npart+1))
### full model predictive score
logevidences1given2

# Prepare table
## prepare table
f <- function(logevid, ro = 1){
  m <- round(mean(logevid), ro)
  mmin <- round(min(logevid), ro)
  mmax <- round(max(logevid), ro)
  print_evid <- paste0(m, " [", mmin, ", ", mmax, "]")
  return(print_evid)
}

library(xtable)
tab <- data.frame(module1 = f(logevidences1),
                  full = f(logevidences1given2),
                  module2 = f(logevidences2alone))
tab <- t(tab)

rownames(tab) <- c("module 1", "full model", "module 2")
colnames(tab) <- "predictive score on $Y_1$"
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T)
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T,
      file = "epidemiology.predictY1.tex")



### Predict Y2
## full model predictive
load("epidemiology_module2givenmodule1.N1024.K100.rep5.RData")
logevidences2given1 <- c()
for (i in 1:length(results2given1)){
  logevidences2given1 <- c(logevidences2given1, results2given1[[i]]$logevidence)
}
logevidences2given1

## plug in predictive
load("epidemiology_module2givenmodule1plugin.N1024.K100.rep5.RData")
logevidences2given1plugin <- c()
for (i in 1:length(results2given1plugin)){
  logevidences2given1plugin <- c(logevidences2given1plugin, results2given1plugin[[i]]$logevidence)
}
logevidences2given1plugin

## predictive ignoring Y1 altogether
load("epidemiology_module2alone.N1024.K100.rep5.RData")
logevidences2alone <- c()
for (i in 1:length(results2given1plugin)){
  logevidences2alone <- c(logevidences2alone, results2alone[[i]]$logevidence)
}
logevidences2alone

## Cut approach
## predictive ignoring Y1 altogether
load("epidemiology_module2givenmodule1cut.N1024.K100.M1024.rep5.RData")
N1 = length(results2given1cut[[1]])
logevidences2cut <- rep(0, 5)
for (irep in 1:5){
  evid <- rep(0, N1)
  for (i in 1:N1){
    evid[i] <- results2given1cut[[irep]][[i]]$sumlogevidence
  }
  logevidences2cut[irep] <- max(evid) + log(mean(exp(evid - max(evid))))
}

logevidences2cut

## prepare table
f <- function(logevid, ro = 1){
  m <- round(mean(logevid), ro)
  mmin <- round(min(logevid), ro)
  mmax <- round(max(logevid), ro)
  print_evid <- paste0(m, " [", mmin, ", ", mmax, "]")
  return(print_evid)
}

print2given1 <- f(logevidences2given1)
print2given1plugin <- f(logevidences2given1plugin)
print2alone <- f(logevidences2alone)
print2cut <- f(logevidences2cut)

library(xtable)
tab <- data.frame(module2 = print2alone,
                  full = print2given1,
                  cut = print2cut,
                  plugin = print2given1plugin)
tab <- t(tab)
xtable(tab)

rownames(tab) <- c("module 2", "full model", "cut approach", "plug-in approach")
colnames(tab) <- "predictive score on $Y_2$"
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T)
print(xtable(tab, align = "ll"), floating = F, sanitize.text.function = function(x) x, include.colnames = T,
      file = "epidemiology.predictY2.tex")



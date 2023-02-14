# Title: Parturition site selection
# Subtitle: 08 - Used-habitat calibration plots
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 29 Sep 2021
# Date completed: 4 Oct 2021
# Date modified: 11 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)       # modeling
library(mefa4)         # %notin% function
library(mvtnorm)       # multivariate normal distribution
library(KernSmooth)    # calculate simulation envelopes
library(cowplot)       # arranging multiple plots

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

load("partsite_models.RData")

#_____________________________________________________________________________________________________________
# 3. Split strata into test and training sets (3 folds) ----
#_____________________________________________________________________________________________________________

fold.size <- length(unique(scaled.data$Animal)) / 3

# fold 1
test.data.1 <- scaled.data %>% filter(Case == 1) %>% 
               slice_sample(n = fold.size)

fold1.test <- scaled.data %>% filter(Animal %in% unique(test.data.1$Animal))

fold1.train <- scaled.data %>% filter(Animal %notin% unique(test.data.1$Animal))

# define animal names
fold1.animals <- unique(fold1.test$Animal)

# fold 2
# extract used locations
test.data.2 <- scaled.data %>% filter(Animal %notin% fold1.animals,
                                      Case == 1) %>% 
                               slice_sample(n = fold.size)

# add random locations
fold2.test <- scaled.data %>% filter(Animal %in% unique(test.data.2$Animal))

# define animal names
fold2.animals <- unique(fold2.test$Animal)

# concatenate animal names from folds 1 and 2
fold1.2.animals <- c(fold1.animals, fold2.animals)

# train data
fold2.train <- scaled.data %>% filter(Animal %notin% fold2.animals)

# fold 3
# extract used locations
test.data.3 <- scaled.data %>% filter(Animal %notin% fold1.2.animals,
                                      Case == 1)

# add random locations
fold3.test <- scaled.data %>% filter(Animal %in% unique(test.data.3$Animal))

# define animal names
fold3.animals <- unique(fold3.test$Animal)

# train data
fold3.train <- scaled.data %>% filter(Animal %notin% fold3.animals)

#_____________________________________________________________________________________________________________
# 4. Fit models to training data ----
#_____________________________________________________________________________________________________________
# 4a. Fold 1 ----
#_____________________________________________________________________________________________________________

# fit model 1
train.model1.1struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 IJI.250 + IJI.250.Q +
                                 ED.500 +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold1.train,
                               doFit = FALSE)

train.model1.1struc$parameters$theta[1] = log(1e3)
train.model1.1struc$mapArg = list(theta = factor(c(NA)))

train.model1.1 <- glmmTMB:::fitTMB(train.model1.1struc)

# fit model 2
train.model1.2struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold1.train,
                               doFit = FALSE)

train.model1.2struc$parameters$theta[1] = log(1e3)
train.model1.2struc$mapArg = list(theta = factor(c(NA)))

train.model1.2 <- glmmTMB:::fitTMB(train.model1.2struc)

#_____________________________________________________________________________________________________________
# 4b. Fold 2 ----
#_____________________________________________________________________________________________________________

# fit model 1
train.model2.1struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 IJI.250 + IJI.250.Q +
                                 ED.500 +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold2.train,
                               doFit = FALSE)

train.model2.1struc$parameters$theta[1] = log(1e3)
train.model2.1struc$mapArg = list(theta = factor(c(NA)))

train.model2.1 <- glmmTMB:::fitTMB(train.model2.1struc)

# fit model 2
train.model2.2struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold2.train,
                               doFit = FALSE)

train.model2.2struc$parameters$theta[1] = log(1e3)
train.model2.2struc$mapArg = list(theta = factor(c(NA)))

train.model2.2 <- glmmTMB:::fitTMB(train.model2.2struc)

#_____________________________________________________________________________________________________________
# 4c. Fold 3 ----
#_____________________________________________________________________________________________________________

# fit model 1
train.model3.1struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 IJI.250 + IJI.250.Q +
                                 ED.500 +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold3.train,
                               doFit = FALSE)

train.model3.1struc$parameters$theta[1] = log(1e3)
train.model3.1struc$mapArg = list(theta = factor(c(NA)))

train.model3.1 <- glmmTMB:::fitTMB(train.model3.1struc)

# fit model 2
train.model3.2struc <- glmmTMB(Case ~ EVI.2pix +
                                 canopy + canopy.Q +
                                 TRI +
                                 TPI.250 +
                                 SRI.500.P +
                                 dRoad.1000 +
                                 (1 | Animal),
                               family = poisson,
                               data = fold3.train,
                               doFit = FALSE)

train.model3.2struc$parameters$theta[1] = log(1e3)
train.model3.2struc$mapArg = list(theta = factor(c(NA)))

train.model3.2 <- glmmTMB:::fitTMB(train.model3.2struc)

#_____________________________________________________________________________________________________________
# 5. Draw possible parameter values from a normal distribution ----
#_____________________________________________________________________________________________________________
# 5a. Fold 1 ----
#_____________________________________________________________________________________________________________

# model 1
param.est1.1 <- train.model1.1$fit$par

vcov.matrix1.1 <- vcov(train.model1.1)[[1]]

# multivariate normal distribution
mvt.norm.values1.1 <- as.data.frame(rmvnorm(n = 500, mean = param.est1.1, sigma = vcov.matrix1.1))

# model 2
param.est1.2 <- train.model1.2$fit$par

vcov.matrix1.2 <- vcov(train.model1.2)[[1]]

# multivariate normal distribution
mvt.norm.values1.2 <- as.data.frame(rmvnorm(n = 500, mean = param.est1.2, sigma = vcov.matrix1.2))

#_____________________________________________________________________________________________________________
# 5b. Fold 2 ----
#_____________________________________________________________________________________________________________

# model 1
param.est2.1 <- train.model2.1$fit$par

vcov.matrix2.1 <- vcov(train.model2.1)[[1]]

# multivariate normal distribution
mvt.norm.values2.1 <- as.data.frame(rmvnorm(n = 500, mean = param.est2.1, sigma = vcov.matrix2.1))

# model 2
param.est2.2 <- train.model2.2$fit$par

vcov.matrix2.2 <- vcov(train.model2.2)[[1]]

# multivariate normal distribution
mvt.norm.values2.2 <- as.data.frame(rmvnorm(n = 500, mean = param.est2.2, sigma = vcov.matrix2.2))

#_____________________________________________________________________________________________________________
# 5c. Fold 3 ----
#_____________________________________________________________________________________________________________

# model 1
param.est3.1 <- train.model3.1$fit$par

vcov.matrix3.1 <- vcov(train.model3.1)[[1]]

# multivariate normal distribution
mvt.norm.values3.1 <- as.data.frame(rmvnorm(n = 500, mean = param.est3.1, sigma = vcov.matrix3.1))

# model 2
param.est3.2 <- train.model3.2$fit$par

vcov.matrix3.2 <- vcov(train.model3.2)[[1]]

# multivariate normal distribution
mvt.norm.values3.2 <- as.data.frame(rmvnorm(n = 500, mean = param.est3.2, sigma = vcov.matrix3.2))

#_____________________________________________________________________________________________________________
# 6. Calculate RSF scores ----
#_____________________________________________________________________________________________________________
# 6a. Fold 1 ----
#_____________________________________________________________________________________________________________

# model 1
used.RSF.scores1.1 <- data.frame("Animal" = fold1.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.1.1 <- data.frame(exp(mvt.norm.values1.1[y, 2]*fold1.test$EVI.2pix +
                                         mvt.norm.values1.1[y, 3]*fold1.test$canopy +
                                         mvt.norm.values1.1[y, 4]*fold1.test$canopy.Q +
                                         mvt.norm.values1.1[y, 5]*fold1.test$IJI.250 +
                                         mvt.norm.values1.1[y, 6]*fold1.test$IJI.250.Q +
                                         mvt.norm.values1.1[y, 7]*fold1.test$ED.500 +
                                         mvt.norm.values1.1[y, 8]*fold1.test$TRI +
                                         mvt.norm.values1.1[y, 9]*fold1.test$TPI.250 +
                                         mvt.norm.values1.1[y, 10]*fold1.test$SRI.500.P +
                                         mvt.norm.values1.1[y, 11]*fold1.test$dRoad.1000))
  
  used.RSF.scores1.1 <- cbind(used.RSF.scores1.1, used.RSF.score.1.1)
  
}

# change column names
colnames(used.RSF.scores1.1) <- c("Animal", 1:500)

# model 2
used.RSF.scores1.2 <- data.frame("Animal" = fold1.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.1.2 <- data.frame(exp(mvt.norm.values1.2[y, 2]*fold1.test$EVI.2pix +
                                         mvt.norm.values1.2[y, 3]*fold1.test$canopy +
                                         mvt.norm.values1.2[y, 4]*fold1.test$canopy.Q +
                                         mvt.norm.values1.2[y, 5]*fold1.test$TRI +
                                         mvt.norm.values1.2[y, 6]*fold1.test$TPI.250 +
                                         mvt.norm.values1.2[y, 7]*fold1.test$SRI.500.P +
                                         mvt.norm.values1.2[y, 8]*fold1.test$dRoad.1000))
  
  used.RSF.scores1.2 <- cbind(used.RSF.scores1.2, used.RSF.score.1.2)
  
}

# change column names
colnames(used.RSF.scores1.2) <- c("Animal", 1:500)

#_____________________________________________________________________________________________________________
# 6b. Fold 2 ----
#_____________________________________________________________________________________________________________

# model 1
used.RSF.scores2.1 <- data.frame("Animal" = fold2.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.2.1 <- data.frame(exp(mvt.norm.values2.1[y, 2]*fold2.test$EVI.2pix +
                                         mvt.norm.values2.1[y, 3]*fold2.test$canopy +
                                         mvt.norm.values2.1[y, 4]*fold2.test$canopy.Q +
                                         mvt.norm.values2.1[y, 5]*fold2.test$IJI.250 +
                                         mvt.norm.values2.1[y, 6]*fold2.test$IJI.250.Q +
                                         mvt.norm.values2.1[y, 7]*fold2.test$ED.500 +
                                         mvt.norm.values2.1[y, 8]*fold2.test$TRI +
                                         mvt.norm.values2.1[y, 9]*fold2.test$TPI.250 +
                                         mvt.norm.values2.1[y, 10]*fold2.test$SRI.500.P +
                                         mvt.norm.values2.1[y, 11]*fold2.test$dRoad.1000))
  
  used.RSF.scores2.1 <- cbind(used.RSF.scores2.1, used.RSF.score.2.1)
  
}

# change column names
colnames(used.RSF.scores2.1) <- c("Animal", 1:500)

# model 2
used.RSF.scores2.2 <- data.frame("Animal" = fold2.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.2.2 <- data.frame(exp(mvt.norm.values2.2[y, 2]*fold2.test$EVI.2pix +
                                         mvt.norm.values2.2[y, 3]*fold2.test$canopy +
                                         mvt.norm.values2.2[y, 4]*fold2.test$canopy.Q +
                                         mvt.norm.values2.2[y, 5]*fold2.test$TRI +
                                         mvt.norm.values2.2[y, 6]*fold2.test$TPI.250 +
                                         mvt.norm.values2.2[y, 7]*fold2.test$SRI.500.P +
                                         mvt.norm.values2.2[y, 8]*fold2.test$dRoad.1000))
  
  used.RSF.scores2.2 <- cbind(used.RSF.scores2.2, used.RSF.score.2.2)
  
}

# change column names
colnames(used.RSF.scores2.2) <- c("Animal", 1:500)

#_____________________________________________________________________________________________________________
# 6c. Fold 3 ----
#_____________________________________________________________________________________________________________

# model 1
used.RSF.scores3.1 <- data.frame("Animal" = fold3.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.3.1 <- data.frame(exp(mvt.norm.values3.1[y, 2]*fold3.test$EVI.2pix +
                                         mvt.norm.values3.1[y, 3]*fold3.test$canopy +
                                         mvt.norm.values3.1[y, 4]*fold3.test$canopy.Q +
                                         mvt.norm.values3.1[y, 5]*fold3.test$IJI.250 +
                                         mvt.norm.values3.1[y, 6]*fold3.test$IJI.250.Q +
                                         mvt.norm.values3.1[y, 7]*fold3.test$ED.500 +
                                         mvt.norm.values3.1[y, 8]*fold3.test$TRI +
                                         mvt.norm.values3.1[y, 9]*fold3.test$TPI.250 +
                                         mvt.norm.values3.1[y, 10]*fold3.test$SRI.500.P +
                                         mvt.norm.values3.1[y, 11]*fold3.test$dRoad.1000))
  
  used.RSF.scores3.1 <- cbind(used.RSF.scores3.1, used.RSF.score.3.1)
  
}

# change column names
colnames(used.RSF.scores3.1) <- c("Animal", 1:500)

# model 2
used.RSF.scores3.2 <- data.frame("Animal" = fold3.test$Animal)

for (y in 1:500) {
  
  used.RSF.score.3.2 <- data.frame(exp(mvt.norm.values3.2[y, 2]*fold3.test$EVI.2pix +
                                         mvt.norm.values3.2[y, 3]*fold3.test$canopy +
                                         mvt.norm.values3.2[y, 4]*fold3.test$canopy.Q +
                                         mvt.norm.values3.2[y, 5]*fold3.test$TRI +
                                         mvt.norm.values3.2[y, 6]*fold3.test$TPI.250 +
                                         mvt.norm.values3.2[y, 7]*fold3.test$SRI.500.P +
                                         mvt.norm.values3.2[y, 8]*fold3.test$dRoad.1000))
  
  used.RSF.scores3.2 <- cbind(used.RSF.scores3.2, used.RSF.score.3.2)
  
}

# change column names
colnames(used.RSF.scores3.2) <- c("Animal", 1:500)

#_____________________________________________________________________________________________________________
# 7. Sample from test data (one from each stratum) ----
#_____________________________________________________________________________________________________________
# 7a. Fold 1 ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows1.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold1.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold1.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores1.1[used.RSF.scores1.1$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows1.1 <- rbind(sampled.rows1.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows1.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold1.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold1.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores1.2[used.RSF.scores1.2$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows1.2 <- rbind(sampled.rows1.2, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7b. Fold 2 ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows2.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold2.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold2.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores2.1[used.RSF.scores2.1$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows2.1 <- rbind(sampled.rows2.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows2.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold2.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold2.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores2.2[used.RSF.scores2.2$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows2.2 <- rbind(sampled.rows2.2, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7c. Fold 3 ----
#_____________________________________________________________________________________________________________

# model 1
sampled.rows3.1 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold3.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold3.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores3.1[used.RSF.scores3.1$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows3.1 <- rbind(sampled.rows3.1, sampled.row)
    
  }
  
}

# model 2
sampled.rows3.2 <- data.frame()

for (y in 1:500) {
  
  value <- y
  
  for (x in unique(fold3.test$Animal)) {
    
    CollarID <- x
    
    stratum.data <- fold3.test %>% filter(Animal == CollarID) %>%
      bind_cols(RSF.score = used.RSF.scores3.2[used.RSF.scores3.2$Animal == CollarID, y + 1])
    
    # randomly sample a row based on probability
    sampled.row <- stratum.data %>% slice_sample(n = 1, weight_by = RSF.score) %>%
      mutate("run" = value)
    
    # bind to master data frame
    sampled.rows3.2 <- rbind(sampled.rows3.2, sampled.row)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 8. Create simulation envelopes ----
#_____________________________________________________________________________________________________________
# 8a. Write function ----
#_____________________________________________________________________________________________________________

sim_env <- function(variable, x.val = 500) {
  
  # create blank df
  sim.env <- data.frame()
  
  # iterate by fold and model
  # fold 1 
  # model 1
  # define range of values
  range.val1.1 <- c(min(sampled.rows1.1[variable]), max(sampled.rows1.1[variable]))
  
  # how many x values?
  dens1.1 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows1.1 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val1.1, gridsize = 500)
    
    dens1.1[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f1.1 <- apply(dens1.1, 2, mean, na.rm = TRUE)
  
  low.f1.1 <- apply(dens1.1, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f1.1 <- apply(dens1.1, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env1.1 <- data.frame("var" = variable,
                           "fold" = "Fold 1",
                           "model" = "Global",
                           "mean" = mean.f1.1,
                           "low" = low.f1.1,
                           "high" = high.f1.1,
                           "x" = seq(min(sampled.rows1.1[variable]), 
                                     max(sampled.rows1.1[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env1.1)
  
  # model 2
  # define range of values
  range.val1.2 <- c(min(sampled.rows1.2[variable]), max(sampled.rows1.2[variable]))
  
  # how many x values?
  dens1.2 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows1.2 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val1.2, gridsize = 500)
    
    dens1.2[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f1.2 <- apply(dens1.2, 2, mean, na.rm = TRUE)
  
  low.f1.2 <- apply(dens1.2, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f1.2 <- apply(dens1.2, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env1.2 <- data.frame("var" = variable,
                           "fold" = "Fold 1",
                           "model" = "FORCOV + TOPO + HUMAN" ,
                           "mean" = mean.f1.2,
                           "low" = low.f1.2,
                           "high" = high.f1.2,
                           "x" = seq(min(sampled.rows1.2[variable]), 
                                     max(sampled.rows1.2[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env1.2)
  
  # fold 2 
  # model 1
  # define range of values
  range.val2.1 <- c(min(sampled.rows2.1[variable]), max(sampled.rows2.1[variable]))
  
  # how many x values?
  dens2.1 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows2.1 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val2.1, gridsize = 500)
    
    dens2.1[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f2.1 <- apply(dens2.1, 2, mean, na.rm = TRUE)
  
  low.f2.1 <- apply(dens2.1, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f2.1 <- apply(dens2.1, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env2.1 <- data.frame("var" = variable,
                           "fold" = "Fold 2",
                           "model" = "Global",
                           "mean" = mean.f2.1,
                           "low" = low.f2.1,
                           "high" = high.f2.1,
                           "x" = seq(min(sampled.rows2.1[variable]), 
                                     max(sampled.rows2.1[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env2.1)
  
  # model 2
  # define range of values
  range.val2.2 <- c(min(sampled.rows2.2[variable]), max(sampled.rows2.2[variable]))
  
  # how many x values?
  dens2.2 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows2.2 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val2.2, gridsize = 500)
    
    dens2.2[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f2.2 <- apply(dens2.2, 2, mean, na.rm = TRUE)
  
  low.f2.2 <- apply(dens2.2, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f2.2 <- apply(dens2.2, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env2.2 <- data.frame("var" = variable,
                           "fold" = "Fold 2",
                           "model" = "FORCOV + TOPO + HUMAN" ,
                           "mean" = mean.f2.2,
                           "low" = low.f2.2,
                           "high" = high.f2.2,
                           "x" = seq(min(sampled.rows2.2[variable]), 
                                     max(sampled.rows2.2[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env2.2)
  
  # fold 3 
  # model 1
  # define range of values
  range.val3.1 <- c(min(sampled.rows3.1[variable]), max(sampled.rows3.1[variable]))
  
  # how many x values?
  dens3.1 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows3.1 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val3.1, gridsize = 500)
    
    dens3.1[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f3.1 <- apply(dens3.1, 2, mean, na.rm = TRUE)
  
  low.f3.1 <- apply(dens3.1, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f3.1 <- apply(dens3.1, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env3.1 <- data.frame("var" = variable,
                           "fold" = "Fold 3",
                           "model" = "Global",
                           "mean" = mean.f3.1,
                           "low" = low.f3.1,
                           "high" = high.f3.1,
                           "x" = seq(min(sampled.rows3.1[variable]), 
                                     max(sampled.rows3.1[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env3.1)
  
  # model 2
  # define range of values
  range.val3.2 <- c(min(sampled.rows3.2[variable]), max(sampled.rows3.2[variable]))
  
  # how many x values?
  dens3.2 <- matrix(NA, x.val, x.val)
  
  # for loop to iterate over 500 possible values
  for (z in 1:500) {
    
    which.run <- z
    
    run.rows <- sampled.rows3.2 %>% filter(run == which.run)
    
    run.dens <- bkde(c(run.rows[variable][,1]), range.x = range.val3.2, gridsize = 500)
    
    dens3.2[z, ] <- run.dens$y
    
  }
  
  # compute 95% simulation envelope
  mean.f3.2 <- apply(dens3.2, 2, mean, na.rm = TRUE)
  
  low.f3.2 <- apply(dens3.2, 2, quantile, prob = 0.025, na.rm = TRUE)
  
  high.f3.2 <- apply(dens3.2, 2, quantile, prob = 0.975, na.rm = TRUE)
  
  # create a df
  sim.env3.2 <- data.frame("var" = variable,
                           "fold" = "Fold 3",
                           "model" = "FORCOV + TOPO + HUMAN" ,
                           "mean" = mean.f3.2,
                           "low" = low.f3.2,
                           "high" = high.f3.2,
                           "x" = seq(min(sampled.rows3.2[variable]), 
                                     max(sampled.rows3.2[variable]), 
                                     length.out = x.val))
  
  # bind to sim.env df
  sim.env <- rbind(sim.env, sim.env3.2)
  
  # factor levels
  sim.env$model <- factor(sim.env$model,
                          levels = c("Global", "FORCOV + TOPO + HUMAN"))
  
  # return the final df
  return(sim.env)
  
}

#_____________________________________________________________________________________________________________
# 8b. Simulation envelopes per variable ----
#_____________________________________________________________________________________________________________

sim.EVI <- sim_env(variable = "EVI.2pix")

sim.canopy <- sim_env(variable = "canopy")

sim.IJI <- sim_env(variable = "IJI.250")

sim.ED <- sim_env(variable = "ED.250")

sim.TRI <- sim_env(variable = "TRI")

sim.TPI <- sim_env(variable = "TPI.250")

sim.SRI <- sim_env(variable = "SRI.500.P")

sim.dRoad <- sim_env(variable = "dRoad.1000")

#_____________________________________________________________________________________________________________
# 9. Used-habitat calibration plots ----

# resolution: 410 x 543

# combine test data to make facetting possible
fold1.test$fold <- "Fold 1"

fold2.test$fold <- "Fold 2"

fold3.test$fold <- "Fold 3"

all.folds.test <- rbind(fold1.test, fold2.test, fold3.test)
#_____________________________________________________________________________________________________________
# 9a. EVI ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.EVI) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(EVI.2pix),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(EVI.2pix),
                    size = 1.25,
                    color = "#66CC33") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(EVI[1000])))

#_____________________________________________________________________________________________________________
# 9b. Canopy ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.canopy) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(canopy),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(canopy),
                    size = 1.25,
                    color = "#66CC33") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(Canopy)))

#_____________________________________________________________________________________________________________
# 9c. IJI ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.IJI) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(IJI.250),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(IJI.250),
                    size = 1.25,
                    color = "#3399FF") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(IJI[250])))

#_____________________________________________________________________________________________________________
# 9d. ED ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.ED) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(ED.500),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(ED.500),
                    size = 1.25,
                    color = "#3399FF") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(ED[500])))

#_____________________________________________________________________________________________________________
# 9e. TRI ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.TRI) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(TRI),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(TRI),
                    size = 1.25,
                    color = "#006633") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(TRI)))

#_____________________________________________________________________________________________________________
# 9f. TPI ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.TPI) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(TPI.250),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(TPI.250),
                    size = 1.25,
                    color = "#006633") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(TPI[250])))

#_____________________________________________________________________________________________________________
# 9g. SRI ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.SRI) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(SRI.500.P),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(SRI.500.P),
                    size = 1.25,
                    color = "#006633") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(SRI[500]))) +
       coord_cartesian(xlim = c(-5, 5))

#_____________________________________________________________________________________________________________
# 9h. dRoad ----
#_____________________________________________________________________________________________________________

# facetted plot
ggplot(data = sim.dRoad) +
       theme_bw() +
       facet_grid(rows = vars(fold),
                  cols = vars(model)) + 
       geom_ribbon(aes(y = mean, 
                       ymin = low, 
                       ymax = high,
                       x = x),
                   fill = "lightgray",
                   alpha = 0.5) +
       geom_density(data = all.folds.test[all.folds.test$Case == 0,], 
                    aes(dRoad.1000),
                    linetype = "dashed") +
       geom_density(data = all.folds.test[all.folds.test$Case == 1,], 
                    aes(dRoad.1000),
                    size = 1.25,
                    color = "#FF3300") +
       theme(panel.grid = element_blank(),
             axis.title.y = element_blank()) +
       xlab(expression(scale(dRoad[1000])))

#_____________________________________________________________________________________________________________
# 10. Save image ----
#_____________________________________________________________________________________________________________

save.image("UHC_3_fold.RData")


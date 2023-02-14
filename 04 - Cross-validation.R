# Title: Parturition site selection
# Subtitle: 04 - Cross-validation
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 15 Apr 2021
# Date completed: 15 Apr 2021
# Date modified: 7 Dec 2021
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)     # modeling

#_____________________________________________________________________________________________________________
# 2. Read in .Rdata ----
#_____________________________________________________________________________________________________________

load("partsite_models.RData")

#_____________________________________________________________________________________________________________
# 3. Cross-validation following Lehman et al. (2016) ----
#_____________________________________________________________________________________________________________

# create blank data.frame
all.scores.1 <- data.frame()

for(x in unique(scaled.data$Animal)) {
  
  elkID <- x
  
  # test data
  test.data <- scaled.data %>% filter(Animal == elkID)
  
  # train data
  train.data <- scaled.data %>% filter(Animal != elkID)
  
  # fit model 1 (global)
  top.model.1struc <- glmmTMB(Case ~ EVI.2pix +
                                     canopy + canopy.Q +
                                     IJI.250 + IJI.250.Q +
                                     ED.500 +
                                     TRI +
                                     TPI.250 +
                                     SRI.500.P +
                                     dRoad.1000 +
                                (1 | Animal),
                              family = poisson,
                              data = train.data,
                              doFit = FALSE)
  
  top.model.1struc$parameters$theta[1] = log(1e3)
  top.model.1struc$mapArg = list(theta = factor(c(NA)))
  
  top.model.1 <- glmmTMB:::fitTMB(top.model.1struc)
  
  # calculate RSF scores on test data
  test.data <- test.data %>% mutate(RSF.score.1 = exp(top.model.1$fit$par[2]*EVI.2pix + 
                                                      top.model.1$fit$par[3]*canopy +
                                                      top.model.1$fit$par[4]*canopy.Q + 
                                                      top.model.1$fit$par[5]*IJI.250 +
                                                      top.model.1$fit$par[6]*IJI.250.Q + 
                                                      top.model.1$fit$par[7]*ED.500 +
                                                      top.model.1$fit$par[8]*TRI + 
                                                      top.model.1$fit$par[9]*TPI.250 +
                                                      top.model.1$fit$par[10]*SRI.500.P + 
                                                      top.model.1$fit$par[11]*dRoad.1000))
  
  # fit model 2 (FORCOV + TOPO + HUMAN)
  top.model.2struc <- glmmTMB(Case ~ EVI.2pix +
                                     canopy + canopy.Q +
                                     TRI +
                                     TPI.250 +
                                     SRI.500.P +
                                     dRoad.1000 +
                                (1 | Animal),
                              family = poisson,
                              data = train.data,
                              doFit = FALSE)
  
  top.model.2struc$parameters$theta[1] = log(1e3)
  top.model.2struc$mapArg = list(theta = factor(c(NA)))
  
  top.model.2 <- glmmTMB:::fitTMB(top.model.2struc)
  
  # calculate RSF scores on test data
  test.data <- test.data %>% mutate(RSF.score.2 = exp(top.model.2$fit$par[2]*EVI.2pix + 
                                                      top.model.2$fit$par[3]*canopy +
                                                      top.model.2$fit$par[4]*canopy.Q + 
                                                      top.model.2$fit$par[5]*TRI +
                                                      top.model.2$fit$par[6]*TPI.250 +
                                                      top.model.2$fit$par[7]*SRI.500.P +
                                                      top.model.2$fit$par[8]*dRoad.1000))
  
  # add weighted RSF scores
  test.data <- test.data %>% mutate(RSF.score.1.w = RSF.score.1*0.79,
                                    RSF.score.2.w = RSF.score.1*0.18)
  
  # create scores data.frame and calculate model-averaged RSF score
  test.scores <- test.data %>% dplyr::select(Animal, Case, RSF.score.1, RSF.score.2, RSF.score.1.w, RSF.score.2.w)
  
  test.scores$RSF.score.avg <- rowSums(test.data[ ,c("RSF.score.1.w", "RSF.score.2.w")])
  
  # create 10 quantiles of RSF scores
  cv.quant.1 <- quantile(test.scores$RSF.score.1, probs = seq(0, 1, 0.1))
  cv.quant.2 <- quantile(test.scores$RSF.score.2, probs = seq(0, 1, 0.1))
  cv.quant.avg <- quantile(test.scores$RSF.score.avg, probs = seq(0, 1, 0.1))
  
  # create a bin variable
  test.scores$bin.1 <- rep(NA, length(test.scores$RSF.score.1))
  test.scores$bin.2 <- rep(NA, length(test.scores$RSF.score.2))
  test.scores$bin.avg <- rep(NA, length(test.scores$RSF.score.avg))
  
  # for loops to classify RSF scores in 1 of 10 bins
  # model 1 scores
  for (j in 1:10){
    
    test.scores$bin.1[test.scores$RSF.score.1 >= cv.quant.1[j] & test.scores$RSF.score.1 < cv.quant.1[j+1]] = j
    
  }
  
  # model 2 scores
  for (j in 1:10){
    
    test.scores$bin.2[test.scores$RSF.score.2 >= cv.quant.2[j] & test.scores$RSF.score.2 < cv.quant.2[j+1]] = j
    
  }
  
  # average scores
  for (j in 1:10){
    
    test.scores$bin.avg[test.scores$RSF.score.avg >= cv.quant.avg[j] & test.scores$RSF.score.avg < cv.quant.avg[j+1]] = j
    
  }
  
  # if there's an NA because of a high RSF score, change the bin to 10
  test.scores$bin.1[is.na(test.scores$bin.1)] <- 10
  test.scores$bin.2[is.na(test.scores$bin.2)] <- 10
  test.scores$bin.avg[is.na(test.scores$bin.avg)] <- 10
  
  # area-adjusted frequency from Roberts et al.
  cv.aa.1 <- table(test.scores$Case, test.scores$bin.1)
  cv.aa.1 <- t(cv.aa.1)
  cv.aa.1 <- as.data.frame.matrix(cv.aa.1)
  cv.aa.1$areaadjusted <- rep(NA, length(10))
  cv.aa.1$areaadjusted <- (cv.aa.1[,2] / sum(cv.aa.1[,2])) / (cv.aa.1[,1] / sum(cv.aa.1[,1]))
  cv.aa.1$bins <- seq(1, 10, by = 1)
  
  cv.aa.2 <- table(test.scores$Case, test.scores$bin.2)
  cv.aa.2 <- t(cv.aa.2)
  cv.aa.2 <- as.data.frame.matrix(cv.aa.2)
  cv.aa.2$areaadjusted <- rep(NA, length(10))
  cv.aa.2$areaadjusted <- (cv.aa.2[,2] / sum(cv.aa.2[,2])) / (cv.aa.2[,1] / sum(cv.aa.2[,1]))
  cv.aa.2$bins <- seq(1, 10, by = 1)
  
  cv.aa.avg <- table(test.scores$Case, test.scores$bin.avg)
  cv.aa.avg <- t(cv.aa.avg)
  cv.aa.avg <- as.data.frame.matrix(cv.aa.avg)
  cv.aa.avg$areaadjusted <- rep(NA, length(10))
  cv.aa.avg$areaadjusted <- (cv.aa.avg[,2] / sum(cv.aa.avg[,2])) / (cv.aa.avg[,1] / sum(cv.aa.avg[,1]))
  cv.aa.avg$bins <- seq(1, 10, by = 1)
  
  # create sum-friendly df
  cv.indiv.1 <- as.data.frame(t(cv.aa.1[ ,c(3)]))
  cv.indiv.2 <- as.data.frame(t(cv.aa.2[ ,c(3)]))
  cv.indiv.avg <- as.data.frame(t(cv.aa.avg[ ,c(3)]))
  
  # add column to label which model scores
  cv.indiv.1$model <- 1
  cv.indiv.2$model <- 2
  cv.indiv.avg$model <- "avg"
  
  # bind to master data frame
  all.scores.1 <- rbind(all.scores.1, cv.indiv.1, cv.indiv.2, cv.indiv.avg)
  
}

#_____________________________________________________________________________________________________________
# 3. AAF vs. bin correlation ----
#_____________________________________________________________________________________________________________

# model 1 
all.scores.model.1 <- all.scores.1 %>% filter(model == "1")

# column sums
sums.model.1.matrix <- cbind(as.vector(colSums(all.scores.model.1[ ,c(1:10)])), 1:10)

sums.model.1.df <- as.data.frame(sums.model.1.matrix)

# model 2 
all.scores.model.2 <- all.scores.1 %>% filter(model == "2")

# column sums
sums.model.2.matrix <- cbind(as.vector(colSums(all.scores.model.2[ ,c(1:10)])), 1:10)

sums.model.2.df <- as.data.frame(sums.model.2.matrix)

# model avg 
all.scores.model.avg <- all.scores.1 %>% filter(model == "avg")

# column sums
sums.model.avg.matrix <- cbind(as.vector(colSums(all.scores.model.avg[ ,c(1:10)])), 1:10)

sums.model.avg.df <- as.data.frame(sums.model.avg.matrix)

# correlation between AAf and bins
cor.test(sums.model.1.df$V1, sums.model.1.df$V2, method = "spearman", correct = FALSE)
cor.test(sums.model.2.df$V1, sums.model.2.df$V2, method = "spearman", correct = FALSE)
cor.test(sums.model.avg.df$V1, sums.model.avg.df$V2, method = "spearman", correct = FALSE)

#_____________________________________________________________________________________________________________
# 4. AAF vs. bin plot ----
#_____________________________________________________________________________________________________________

# add "model" column and bind together
sums.model.1.df$Model <- "1"
sums.model.2.df$Model <- "2"
sums.model.avg.df$Model <- "avg"

sums.model.all <- rbind(sums.model.1.df, sums.model.2.df, sums.model.avg.df)

sums.model.all$Model <- factor(sums.model.all$Model,
                               levels = c("1", "2", "avg"),
                               labels = c("Global", "FORCOV + TOPO + HUMAN", "Average"))

# plot
ggplot(data = sums.model.all, 
       aes(x = V2, 
           y = V1, 
           group = Model, 
           color = Model,
           linetype = Model)) +
  theme_bw() +
  geom_line(size = 1.25) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.35, 0.8)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  xlab("RSF score bin") +
  ylab("Sum of area-adjusted frequencies of used locations") +
  scale_color_manual(values = c("#eb1e2c", "#f9a729", "black")) +
  scale_linetype_manual(values = c("twodash", "twodash", "solid"))

#_____________________________________________________________________________________________________________
# 5. Save image ----
#_____________________________________________________________________________________________________________

save.image("partsite_cv.RData")

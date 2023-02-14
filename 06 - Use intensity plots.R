# Title: Parturition site selection
# Subtitle: 06 Use intensity plots
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began:  24 Aug 2021
# Date completed: 22 Sep 2021
# Date modified: 11 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)     # modeling
library(cowplot)     # arranging multiple plots

#_____________________________________________________________________________________________________________
# 2. Read in .Rdata ----
#_____________________________________________________________________________________________________________

load("partsite_models.RData")

top.model.1 <- final.models[[15]]
top.model.2 <- final.models[[13]]

#_____________________________________________________________________________________________________________
# 3. Calculate RSF scores and paste to unscaled dataset ----
#_____________________________________________________________________________________________________________

all.data <- all.data %>% mutate(RSF.score.1 = exp(top.model.1$fit$par[2]*scaled.data$EVI.2pix +
                                                   top.model.1$fit$par[3]*scaled.data$canopy +
                                                   top.model.1$fit$par[4]*scaled.data$canopy.Q +
                                                   top.model.1$fit$par[5]*scaled.data$IJI.250 +
                                                   top.model.1$fit$par[6]*scaled.data$IJI.250.Q +
                                                   top.model.1$fit$par[7]*scaled.data$ED.500 +
                                                   top.model.1$fit$par[8]*scaled.data$TRI +
                                                   top.model.1$fit$par[9]*scaled.data$TPI.250 +
                                                   top.model.1$fit$par[10]*scaled.data$SRI.500.P +
                                                   top.model.1$fit$par[11]*scaled.data$dRoad.1000),
                                 RSF.score.2 = exp(top.model.2$fit$par[2]*scaled.data$EVI.2pix +
                                                   top.model.2$fit$par[3]*scaled.data$canopy +
                                                   top.model.2$fit$par[4]*scaled.data$canopy.Q +
                                                   top.model.2$fit$par[5]*scaled.data$TRI +
                                                   top.model.2$fit$par[6]*scaled.data$TPI.250 +
                                                   top.model.2$fit$par[7]*scaled.data$SRI +
                                                   top.model.2$fit$par[8]*scaled.data$dRoad.1000)) %>%
                          mutate(RSF.score.avg = (RSF.score.1*0.79 + RSF.score.2*0.18))

# subset for only available points
avail.data <- all.data %>% filter(Case == 0)

# used data for rugs
used.data <- all.data %>% filter(Case == 1)

#_____________________________________________________________________________________________________________
# 4. Create individual plots ----
#_____________________________________________________________________________________________________________
# 4a. EVI ----
#_____________________________________________________________________________________________________________

plot.EVI <- ggplot(data = avail.data, aes(x = EVI.2pix, y = RSF.score.avg)) +
                      geom_vline(xintercept = mean(avail.data$EVI.1pix),
                                 linetype = "dashed") +
                      theme_bw() +
                      geom_smooth(method = "gam",
                                  color = "#66CC33",
                                  fill = "#66CC33",
                                  size = 1.5,
                                  alpha = 0.25) +
                      xlab(expression(EVI[1000])) +
                      theme(panel.grid = element_blank(), 
                            axis.title.y = element_blank()) +
                      scale_x_continuous(breaks = seq(3000, 9000, 2000),
                                         labels = seq(3000, 9000, 2000)/10000) +
                      scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.EVI

#_____________________________________________________________________________________________________________
# 4b. canopy ----
#_____________________________________________________________________________________________________________

plot.canopy <- ggplot(data = avail.data, aes(x = canopy, y = RSF.score.avg)) +
                   geom_vline(xintercept = mean(avail.data$canopy),
                              linetype = "dashed") +
                   theme_bw() +
                   geom_smooth(method = "gam",
                               color = "#66CC33",
                               fill = "#66CC33",
                               size = 1.5,
                               alpha = 0.25) +
                   xlab(expression(Canopy)) +
                   theme(panel.grid = element_blank(), 
                         axis.title.y = element_blank(),
                         axis.text.y = element_blank()) +
                   scale_x_continuous(breaks = seq(0, 100, 25),
                                      labels = seq(0, 100, 25)) +
                   scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.canopy

#_____________________________________________________________________________________________________________
# 4c. IJI.250 ----
#_____________________________________________________________________________________________________________

plot.IJI <- ggplot(data = avail.data, aes(x = IJI.250, y = RSF.score.avg)) +
                      geom_vline(xintercept = mean(avail.data$IJI.250),
                                 linetype = "dashed") +
                      theme_bw() +
                      geom_smooth(method = "gam",
                                  color = "#3399FF",
                                  fill = "#3399FF",
                                  size = 1.5,
                                  alpha = 0.25) +
                      xlab(expression(IJI[250])) +
                      theme(panel.grid = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank()) +
                      scale_x_continuous(breaks = seq(0, 100, 25)) +
                      scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.IJI

#_____________________________________________________________________________________________________________
# 4d. ED.500 ----
#_____________________________________________________________________________________________________________

plot.ED <- ggplot(data = avail.data, aes(x = ED.500, y = RSF.score.avg)) +
                   geom_vline(xintercept = mean(avail.data$ED.500),
                              linetype = "dashed") +
                   theme_bw() +
                   geom_smooth(method = "gam",
                               color = "#3399FF",
                               fill = "#3399FF",
                               size = 1.5,
                               alpha = 0.25) +
                   xlab(expression(ED[500])) +
                   theme(panel.grid = element_blank(),
                         axis.title.y = element_blank(),
                         axis.text.y = element_blank()) +
                   scale_x_continuous(breaks = seq(0, 10, 2)) +
                   scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.ED

#_____________________________________________________________________________________________________________
# 4e. TRI ----
#_____________________________________________________________________________________________________________

plot.TRI <- ggplot(data = avail.data, aes(x = TRI, y = RSF.score.avg)) +
                      geom_vline(xintercept = mean(avail.data$TRI),
                                 linetype = "dashed") +
                      theme_bw() +
                      geom_smooth(method = "gam",
                                  color = "#006633",
                                  fill = "#006633",
                                  size = 1.5,
                                  alpha = 0.25) +
                      xlab(expression(TRI)) +
                      theme(panel.grid = element_blank(),
                            axis.title.y = element_blank()) +
                      scale_x_continuous(breaks = seq(0, 15, 3)) +
                      scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.TRI

#_____________________________________________________________________________________________________________
# 4f. TPI.250 ----
#_____________________________________________________________________________________________________________

plot.TPI <- ggplot(data = avail.data, aes(x = TPI.250, y = RSF.score.avg)) +
                      geom_vline(xintercept = mean(avail.data$TPI.250),
                                 linetype = "dashed") +
                      theme_bw() +
                      geom_smooth(method = "gam",
                                  color = "#006633",
                                  fill = "#006633",
                                  size = 1.5,
                                  alpha = 0.25) +
                      xlab(expression(TPI[250])) +
                      theme(panel.grid = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank()) +
                      scale_x_continuous(breaks = seq(-4, 4, 2)) +
                      scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.TPI

#_____________________________________________________________________________________________________________
# 4g. SRI ----
#_____________________________________________________________________________________________________________

plot.SRI <- ggplot(data = avail.data, aes(x = SRI.500, y = RSF.score.avg)) +
                      geom_vline(xintercept = mean(avail.data$SRI),
                                 linetype = "dashed") +
                      theme_bw() +
                      geom_smooth(method = "gam",
                                  color = "#006633",
                                  fill = "#006633",
                                  size = 1.5,
                                  alpha = 0.25) +
                      xlab(expression(SRI[500])) +
                      theme(panel.grid = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y = element_blank()) +
                      scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.SRI

#_____________________________________________________________________________________________________________
# 4h. dRoad.1000 ----
#_____________________________________________________________________________________________________________

plot.dRoad <- ggplot(data = avail.data, aes(x = dRoad.1000, y = RSF.score.avg)) +
                     geom_vline(xintercept = mean(avail.data$dRoad.1000),
                                linetype = "dashed") +
                     theme_bw() +
                     geom_smooth(method = "gam",
                                 color = "#FF3300",
                                 fill = "#FF3300",
                                 size = 1.5,
                                 alpha = 0.25) +
                     xlab(expression(dRoad[1000])) +
                     theme(panel.grid = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.y = element_blank()) +
                     scale_x_continuous(breaks = seq(0, 8000, 2500)) +
                     scale_y_continuous(breaks = seq(0, 9, 3)) +
                     coord_cartesian(ylim = c(0, 10))

plot.dRoad

#_____________________________________________________________________________________________________________
# 5. Combine plots ----
#_____________________________________________________________________________________________________________

plot_grid(plot.EVI, plot.canopy, plot.IJI, plot.ED, plot.TRI, plot.TPI, plot.SRI, plot.dRoad,
          nrow = 2, ncol = 4,
          rel_widths = c(1.1, 1, 1, 1))


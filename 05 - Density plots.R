# Title: Parturition site selection
# Subtitle: 05 - Density plots of final variables
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 1 Dec 2021
# Date completed: 1 Dec 2021
# Date modified: 11 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)     # modeling
library(cowplot)     # arranging multiple plots
library(shades)      # darker shades

#_____________________________________________________________________________________________________________
# 2. Read in .Rdata ----
#_____________________________________________________________________________________________________________

load("partsite_models.RData")

#_____________________________________________________________________________________________________________
# 3. Create individual plots ----
#_____________________________________________________________________________________________________________
# 3a. EVI ----
#_____________________________________________________________________________________________________________

plot.EVI <- ggplot(data = all.data, aes(x = EVI.2pix)) +
                   theme_bw() +
                   geom_density(aes(linetype = as.factor(Case),
                                    fill = as.factor(Case)),
                                alpha = 0.5) +
                   scale_fill_manual(values = c("darkgray", "#66CC33")) +
                   scale_linetype_manual(values = c("dashed", "solid")) +
                   xlab(expression(EVI[1000])) +
                   theme(panel.grid = element_blank(), 
                         axis.title.y = element_blank(),
                         legend.position = "none") +
                   scale_x_continuous(breaks = seq(3000, 9000, 2000),
                                      labels = seq(3000, 9000, 2000)/10000)

plot.EVI

#_____________________________________________________________________________________________________________
# 3b. canopy ----
#_____________________________________________________________________________________________________________

plot.canopy <- ggplot(data = all.data, aes(x = canopy)) +
                      theme_bw() +
                      geom_density(aes(linetype = as.factor(Case),
                                       fill = as.factor(Case)),
                                   alpha = 0.5) +
                      scale_fill_manual(values = c("darkgray", "#66CC33")) +
                      scale_linetype_manual(values = c("dashed", "solid")) +
                      xlab(expression(Canopy)) +
                      theme(panel.grid = element_blank(), 
                            axis.title.y = element_blank(),
                            legend.position = "none")

plot.canopy

#_____________________________________________________________________________________________________________
# 3c. IJI.250 ----
#_____________________________________________________________________________________________________________

plot.IJI <- ggplot(data = all.data, aes(x = IJI.250)) +
                      theme_bw() +
                      geom_density(aes(linetype = as.factor(Case),
                                       fill = as.factor(Case)),
                                   alpha = 0.5) +
                      scale_fill_manual(values = c("darkgray", "#3399FF")) +
                      scale_linetype_manual(values = c("dashed", "solid")) +
                      xlab(expression(IJI[250])) +
                      theme(panel.grid = element_blank(), 
                            axis.title.y = element_blank(),
                            legend.position = "none")

plot.IJI

#_____________________________________________________________________________________________________________
# 3d. ED.500 ----
#_____________________________________________________________________________________________________________

plot.ED <- ggplot(data = all.data, aes(x = ED.500)) +
                  theme_bw() +
                  geom_density(aes(linetype = as.factor(Case),
                                   fill = as.factor(Case)),
                               alpha = 0.5) +
                  scale_fill_manual(values = c("darkgray", "#3399FF")) +
                  scale_linetype_manual(values = c("dashed", "solid")) +
                  xlab(expression(ED[500])) +
                  theme(panel.grid = element_blank(), 
                        axis.title.y = element_blank(),
                        legend.position = "none")

plot.ED

#_____________________________________________________________________________________________________________
# 3e. TRI ----
#_____________________________________________________________________________________________________________

plot.TRI <- ggplot(data = all.data, aes(x = TRI)) +
                   theme_bw() +
                   geom_density(aes(linetype = as.factor(Case),
                                    fill = as.factor(Case)),
                                alpha = 0.5) +
                   scale_fill_manual(values = c("darkgray", "#006633")) +
                   scale_linetype_manual(values = c("dashed", "solid")) +
                   xlab(expression(TRI)) +
                   theme(panel.grid = element_blank(), 
                         axis.title.y = element_blank(),
                         legend.position = "none")

plot.TRI

#_____________________________________________________________________________________________________________
# 3f. TPI.250 ----
#_____________________________________________________________________________________________________________

plot.TPI <- ggplot(data = all.data, aes(x = TPI.250)) +
                   theme_bw() +
                   geom_density(aes(linetype = as.factor(Case),
                                    fill = as.factor(Case)),
                                alpha = 0.5) +
                   scale_fill_manual(values = c("darkgray", "#006633")) +
                   scale_linetype_manual(values = c("dashed", "solid")) +
                   xlab(expression(TPI[250])) +
                   theme(panel.grid = element_blank(), 
                         axis.title.y = element_blank(),
                         legend.position = "none") +
                   scale_x_continuous(breaks = seq(-4, 4, 2))
                 
plot.TPI

#_____________________________________________________________________________________________________________
# 3g. SRI.500 ----
#_____________________________________________________________________________________________________________

plot.SRI <- ggplot(data = all.data, aes(x = SRI.500)) +
                   theme_bw() +
                   geom_density(aes(linetype = as.factor(Case),
                                    fill = as.factor(Case)),
                                alpha = 0.5) +
                   scale_fill_manual(values = c("darkgray", "#006633")) +
                   scale_linetype_manual(values = c("dashed", "solid")) +
                   xlab(expression(SRI[500])) +
                   theme(panel.grid = element_blank(), 
                         axis.title.y = element_blank(),
                         legend.position = "none")

plot.SRI

#_____________________________________________________________________________________________________________
# 3h. dRoad.1000 ----
#_____________________________________________________________________________________________________________

plot.dRoad <- ggplot(data = all.data, aes(x = dRoad.1000)) +
                     theme_bw() +
                     geom_density(aes(linetype = as.factor(Case),
                                      fill = as.factor(Case)),
                                  alpha = 0.5) +
                     scale_fill_manual(values = c("darkgray", "#FF3300")) +
                     scale_linetype_manual(values = c("dashed", "solid")) +
                     xlab(expression(dRoad[1000])) +
                     theme(panel.grid = element_blank(), 
                           axis.title.y = element_blank(),
                           legend.position = "none") +
                    scale_x_continuous(breaks = seq(0, 8000, 2500))

plot.dRoad

#_____________________________________________________________________________________________________________
# 5. Combine plots ----
#_____________________________________________________________________________________________________________

plot_grid(plot.EVI, plot.canopy, plot.IJI, plot.ED, plot.TRI, plot.TPI, plot.SRI, plot.dRoad,
          nrow = 2, ncol = 4)

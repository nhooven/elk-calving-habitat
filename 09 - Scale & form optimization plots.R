# Title: Parturition site selection
# Subtitle: 09 - Scale & form optimization plots
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 22 Sep 2021
# Date completed: 11 Oct 2021
# Date modified: 2 Dec 2021
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(cowplot)

#_____________________________________________________________________________________________________________
# 2. Read in dataset ----
#_____________________________________________________________________________________________________________

optim.aic <- read.csv("scale_form_optim.csv")

#_____________________________________________________________________________________________________________
# 3. Reorder factor levels ----
#_____________________________________________________________________________________________________________

optim.aic$Variable <- factor(optim.aic$Variable,
                             levels = c("AI.MF", "Canopy", "Open", "MatureForest", "ED",
                                        "dRoad", "IJI", 
                                        "TPI", "TRI", "SRI", "EVI"))

optim.aic$Scale <- factor(optim.aic$Scale,
                          levels = c("point", "250", "500", "1000"))

optim.aic$Form <- factor(optim.aic$Form,
                          levels = c("linear", "quadratic", "pseudothreshold"))

#_____________________________________________________________________________________________________________
# 4. Create facetted plot (delta AICc) ----
#_____________________________________________________________________________________________________________

ggplot(data = optim.aic, aes(x = factor(Scale), 
                             y = Delta.AICc, 
                             group = Form,
                             color = Form,
                             fill = Form,
                             shape = Form)) +
       theme_bw() +
       facet_wrap(~Variable, scales = "free_y") +
       geom_line(size = 1.25) +
       geom_point(size = 2.5,
                  color = "black") +
       theme(axis.title.x = element_blank(),
             legend.background = element_rect(fill = alpha("white", 0)),
             legend.position = c(0.87, 0.13)) +
       ylab(paste("\U0394", "AICc")) +
       scale_fill_brewer(palette = "Set1") +
       scale_color_brewer(palette = "Set1") +
       scale_y_reverse() +
       scale_shape_manual(values = c(21, 22, 24))

#_____________________________________________________________________________________________________________
# 5. Create facetted plot (AICc) ----
#_____________________________________________________________________________________________________________

ggplot(data = optim.aic, aes(x = factor(Scale), 
                             y = AICc, 
                             group = Form,
                             color = Form,
                             fill = Form,
                             shape = Form)) +
        theme_bw() +
        facet_wrap(~Variable) +
        geom_hline(yintercept = 1254.59) +
        geom_line(size = 1.25) +
        geom_point(size = 2.5,
                   color = "black") +
        theme(axis.title.x = element_blank(),
              legend.background = element_rect(fill = alpha("white", 0)),
              legend.position = c(0.87, 0.13)) +
        ylab(paste("AICc")) +
        scale_fill_brewer(palette = "Set1") +
        scale_color_brewer(palette = "Set1") +
        scale_y_reverse() +
        scale_shape_manual(values = c(21, 22, 24))

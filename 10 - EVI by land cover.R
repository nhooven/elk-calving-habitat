# Title: Parturition site selection
# Subtitle: 10 - EVI by land cover class
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 26 Oct 2021
# Date completed: 26 Oct 2021
# Date modified: 13 Feb 2023
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(mefa4)

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

zonal <- read.csv("Zonal tables.csv")

zonal$MEAN.prop <- zonal$MEAN / 10000

zonal.mean <- zonal %>% filter(Date %in% c("2020_mean", "2021_mean"))

zonal.weeks <- zonal %>% filter(Date %notin% c("2020_mean", "2021_mean"))

#_____________________________________________________________________________________________________________
# 3. Plot ----
#_____________________________________________________________________________________________________________

ggplot() +
  theme_bw() +
  facet_wrap(~Year, nrow = 2) +
  geom_line(data = zonal.weeks, aes(x = Week, 
                                    y = MEAN.prop, 
                                    color = CLASS),
            size = 1.5) +
  geom_point(data = zonal.weeks, aes(x = Week, 
                                     y = MEAN.prop, 
                                     fill = CLASS),
             size = 3,
             shape = 21) +
  scale_color_manual(values = c("gray", "red", "darkgreen", "orange", "blue", "green")) +
  scale_fill_manual(values = c("gray", "red", "darkgreen", "orange", "blue", "green")) +
  scale_x_continuous(breaks = c(1:6),
                     labels = c("5-8", "5-24", "6-9", "6-25", "7-11", "7-27")) +
  ylab("Mean EVI by class")

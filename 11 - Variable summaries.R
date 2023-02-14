# Title: Parturition site selection
# Subtitle: 11 - Variable summaries
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 3 Dec 2021
# Date completed: 3 Dec 2021
# Date modified: 15 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(plotrix)

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

all.data20.21 <- read.csv("part_site_data.csv")
all.data22 <- read.csv("part_site_data_2022.csv")

all.data20.21.1 <- all.data20.21 %>% dplyr::select(x, y, t, Animal, Year, Case, 
                                                   dRoad, dOpen, dMatureForest,
                                                   TPI, TRI, SRI, canopy, 
                                                   dRoad.250, open.250, mf.250,
                                                   TPI.250, TRI.250, SRI.250, canopy.250,
                                                   IJI.250, ED.250, 
                                                   dRoad.500, open.500, mf.500,
                                                   TPI.500, TRI.500, SRI.500, canopy.500,
                                                   IJI.500, ED.500, 
                                                   dRoad.1000, open.1000, mf.1000,
                                                   TPI.1000, TRI.1000, SRI.1000, canopy.1000,
                                                   IJI.1000, ED.1000, 
                                                   EVI, EVI.1pix, EVI.2pix,
                                                   ai.mf.250, ai.mf.500, ai.mf.1000)

all.data22.1 <- all.data22 %>% dplyr::select(x, y, t, Animal, Year, Case, 
                                             dRoad, dOpen, dMatureForest,
                                             TPI, TRI, SRI, canopy, 
                                             dRoad.250, open.250, mf.250,
                                             TPI.250, TRI.250, SRI.250, canopy.250,
                                             IJI.250, ED.250, 
                                             dRoad.500, open.500, mf.500,
                                             TPI.500, TRI.500, SRI.500, canopy.500,
                                             IJI.500, ED.500, 
                                             dRoad.1000, open.1000, mf.1000,
                                             TPI.1000, TRI.1000, SRI.1000, canopy.1000,
                                             IJI.1000, ED.1000, 
                                             EVI, EVI.1pix, EVI.2pix,
                                             ai.mf.250, ai.mf.500, ai.mf.1000)

# bind together
all.data <- rbind(all.data20.21.1, all.data22.1)

# scale EVI variables correctly
all.data$EVI <- all.data$EVI / 10000
all.data$EVI.1pix <- all.data$EVI.1pix / 10000
all.data$EVI.2pix <- all.data$EVI.2pix / 10000

#_____________________________________________________________________________________________________________
# 3. Group by Case and summarize ----
#_____________________________________________________________________________________________________________

mean.data <- all.data %>% group_by(Case) %>%
                             summarize(across(.cols = dRoad:ai.mf.1000,
                                              mean,
                                              na.rm = TRUE)) %>%
                          mutate("stat" = "mean")

sd.data <- all.data %>% group_by(Case) %>%
                          summarize(across(.cols = dRoad:ai.mf.1000,
                                           sd,
                                           na.rm = TRUE)) %>%
                         mutate("stat" = "sd")

summ.data <- rbind(mean.data, sd.data)   

summ.data.t <- t(summ.data)

write.table(summ.data.t, "clipboard", sep = "\t")

#_____________________________________________________________________________________________________________
# 4. Range used for analysis ----
#_____________________________________________________________________________________________________________

range.data <- all.data %>% summarize(across(.cols = dRoad:ai.mf.1000,
                                            range,
                                            na.rm = TRUE))

range.data.t <- t(range.data)

write.table(range.data.t, "clipboard", sep = "\t")

#_____________________________________________________________________________________________________________
# 5. Number of available locations ----
#_____________________________________________________________________________________________________________

all.data.avail <- all.data %>% filter(Case == 0)

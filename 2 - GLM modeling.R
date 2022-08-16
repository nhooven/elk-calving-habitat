# Title: Parturition site selection
# Subtitle: 2 - GLM modeling (MCP)
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 14 Apr 2021
# Date completed: 15 Apr 2021
# Date modified: 9 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(glmmTMB)     # modeling
library(AICcmodavg)  # AIC table
library(performance) # VIF

#_____________________________________________________________________________________________________________
# 2. Read in and combine data  ----
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

#_____________________________________________________________________________________________________________
# 4. Transformations and standardized variables (create new data.frame)  ----
#_____________________________________________________________________________________________________________

# pseudothreshold transformations
all.data <- all.data %>% mutate(dRoad.P = log(dRoad),
                                dOpen.P = log(dOpen),
                                dMatureForest.P = log(dMatureForest),
                                TPI.P = log(TPI - min(TPI)),
                                SRI.P = log(SRI - min(SRI)),
                                TRI.P = log(TRI),
                                canopy.P = log(canopy),
                                dRoad.250.P = log(dRoad.250),
                                open.250.P = log(open.250),
                                mf.250.P = log(mf.250),
                                TPI.250.P = log(TPI.250 - min(TPI.250)),
                                SRI.250.P = log(SRI.250 - min(SRI.250)),
                                TRI.250.P = log(TRI.250),
                                canopy.250.P = log(canopy.250),
                                IJI.250.P = log(IJI.250),
                                ED.250.P = log(ED.250),
                                dRoad.500.P = log(dRoad.500),
                                open.500.P = log(open.500),
                                mf.500.P = log(mf.500),
                                TPI.500.P = log(TPI.500 - min(TPI.500)),
                                SRI.500.P = log(SRI.500 - min(SRI.500)),
                                TRI.500.P = log(TRI.500),
                                canopy.500.P = log(canopy.500),
                                IJI.500.P = log(IJI.500),
                                ED.500.P = log(ED.500),
                                dRoad.1000.P = log(dRoad.1000),
                                open.1000.P = log(open.1000),
                                mf.1000.P = log(mf.1000),
                                TPI.1000.P = log(TPI.1000 - min(TPI.1000)),
                                SRI.1000.P = log(SRI.1000 - min(SRI.1000)),
                                TRI.1000.P = log(TRI.1000),
                                canopy.1000.P = log(canopy.1000),
                                IJI.1000.P = log(IJI.1000),
                                ED.1000.P = log(ED.1000),
                                EVI.P = log(EVI),
                                EVI.1pix.P = log(EVI.1pix),
                                EVI.2pix.P = log(EVI.2pix),
                                ai.mf.250.P = log(ai.mf.250), 
                                ai.mf.500.P = log(ai.mf.500), 
                                ai.mf.1000.P = log(ai.mf.1000))

# replace -Inf with 0s
all.data[mapply(is.infinite, all.data)] <- 0

# quadratic transformations
all.data <- all.data %>% mutate(dRoad.Q = (dRoad)^2,
                                 dOpen.Q = (dOpen)^2,
                                 dMatureForest.Q = (dMatureForest)^2,
                                 TPI.Q = (TPI)^2,
                                 SRI.Q = (SRI)^2,
                                 TRI.Q = (TRI)^2,
                                 canopy.Q = (canopy)^2,
                                 dRoad.250.Q = (dRoad.250)^2,
                                 open.250.Q = (open.250)^2,
                                 mf.250.Q = (mf.250)^2,
                                 TPI.250.Q = (TPI.250)^2,
                                 SRI.250.Q = (SRI.250)^2,
                                 TRI.250.Q = (TRI.250)^2,
                                 canopy.250.Q = (canopy.250)^2,
                                 IJI.250.Q = (IJI.250)^2,
                                 ED.250.Q = (ED.250)^2,
                                 dRoad.500.Q = (dRoad.500)^2,
                                 open.500.Q = (open.500)^2,
                                 mf.500.Q = (mf.500)^2,
                                 TPI.500.Q = (TPI.500)^2,
                                 SRI.500.Q = (SRI.500)^2,
                                 TRI.500.Q = (TRI.500)^2,
                                 canopy.500.Q = (canopy.500)^2,
                                 IJI.500.Q = (IJI.500)^2,
                                 ED.500.Q = (ED.500)^2,
                                 dRoad.1000.Q = (dRoad.1000)^2,
                                 open.1000.Q = (open.1000)^2,
                                 mf.1000.Q = (mf.1000)^2,
                                 TPI.1000.Q = (TPI.1000)^2,
                                 SRI.1000.Q = (SRI.1000)^2,
                                 TRI.1000.Q = (TRI.1000)^2,
                                 canopy.1000.Q = (canopy.1000)^2,
                                 IJI.1000.Q = (IJI.1000)^2,
                                 ED.1000.Q = (ED.1000)^2,
                                 EVI.Q = (EVI)^2,
                                 EVI.1pix.Q = (EVI.1pix)^2,
                                 EVI.2pix.Q = (EVI.2pix)^2,
                                 ai.mf.250.Q = (ai.mf.250)^2, 
                                 ai.mf.500.Q = (ai.mf.500)^2, 
                                 ai.mf.1000.Q = (ai.mf.1000)^2)

scaled.data <- all.data %>% mutate_at(c(7:126), scale)

#_____________________________________________________________________________________________________________
# 5. Scale and functional form optimization ----

modnames.1 <- c("point.linear", "point.log", "point.quad", "250.linear", "250.log", "250.quad",
                "500.linear", "500.log", "500.quad", "1000.linear", "1000.log", "1000.quad",
                "null")

modnames.2 <- c("250.linear", "250.log", "250.quad", "500.linear", "500.log", "500.quad", 
                "1000.linear", "1000.log", "1000.quad", "null")
#_____________________________________________________________________________________________________________
# 5b. dRoad ----
#_____________________________________________________________________________________________________________

mod.dRoad <- list()

# point - linear
mod.dRoad1.struc <- glmmTMB(Case ~ dRoad + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad1.struc$parameters$theta[1] = log(1e3)
mod.dRoad1.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[1]] <- glmmTMB:::fitTMB(mod.dRoad1.struc)

# point - log
mod.dRoad2.struc <- glmmTMB(Case ~ dRoad.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad2.struc$parameters$theta[1] = log(1e3)
mod.dRoad2.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[2]] <- glmmTMB:::fitTMB(mod.dRoad2.struc)

# point - quad
mod.dRoad3.struc <- glmmTMB(Case ~ dRoad + dRoad.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad3.struc$parameters$theta[1] = log(1e3)
mod.dRoad3.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[3]] <- glmmTMB:::fitTMB(mod.dRoad3.struc)

# 250 - linear
mod.dRoad4.struc <- glmmTMB(Case ~ dRoad.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad4.struc$parameters$theta[1] = log(1e3)
mod.dRoad4.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[4]] <- glmmTMB:::fitTMB(mod.dRoad4.struc)

# 250 - log
mod.dRoad5.struc <- glmmTMB(Case ~ dRoad.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad5.struc$parameters$theta[1] = log(1e3)
mod.dRoad5.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[5]] <- glmmTMB:::fitTMB(mod.dRoad5.struc)

# 250 - quad
mod.dRoad6.struc <- glmmTMB(Case ~ dRoad.250 + dRoad.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad6.struc$parameters$theta[1] = log(1e3)
mod.dRoad6.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[6]] <- glmmTMB:::fitTMB(mod.dRoad6.struc)

# 500 - linear
mod.dRoad7.struc <- glmmTMB(Case ~ dRoad.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad7.struc$parameters$theta[1] = log(1e3)
mod.dRoad7.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[7]] <- glmmTMB:::fitTMB(mod.dRoad7.struc)

# 500 - log
mod.dRoad8.struc <- glmmTMB(Case ~ dRoad.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad8.struc$parameters$theta[1] = log(1e3)
mod.dRoad8.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[8]] <- glmmTMB:::fitTMB(mod.dRoad8.struc)

# 500 - quad
mod.dRoad9.struc <- glmmTMB(Case ~ dRoad.500 + dRoad.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad9.struc$parameters$theta[1] = log(1e3)
mod.dRoad9.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[9]] <- glmmTMB:::fitTMB(mod.dRoad9.struc)

# 1000 - linear
mod.dRoad10.struc <- glmmTMB(Case ~ dRoad.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad10.struc$parameters$theta[1] = log(1e3)
mod.dRoad10.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[10]] <- glmmTMB:::fitTMB(mod.dRoad10.struc)

# 1000 - log
mod.dRoad11.struc <- glmmTMB(Case ~ dRoad.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad11.struc$parameters$theta[1] = log(1e3)
mod.dRoad11.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[11]] <- glmmTMB:::fitTMB(mod.dRoad11.struc)

# 1000 - quad
mod.dRoad12.struc <- glmmTMB(Case ~ dRoad.1000 + dRoad.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad12.struc$parameters$theta[1] = log(1e3)
mod.dRoad12.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[12]] <- glmmTMB:::fitTMB(mod.dRoad12.struc)

# null
mod.dRoad13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dRoad13.struc$parameters$theta[1] = log(1e3)
mod.dRoad13.struc$mapArg = list(theta = factor(c(NA)))

mod.dRoad[[13]] <- glmmTMB:::fitTMB(mod.dRoad13.struc)

aictab(mod.dRoad, modnames = modnames.1)

write.table(aictab(mod.dRoad, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.dRoad[[10]])

#_____________________________________________________________________________________________________________
# 5c. dOpen ----
#_____________________________________________________________________________________________________________

mod.dOpen <- list()

# point - linear
mod.dOpen1.struc <- glmmTMB(Case ~ dOpen + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen1.struc$parameters$theta[1] = log(1e3)
mod.dOpen1.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[1]] <- glmmTMB:::fitTMB(mod.dOpen1.struc)

# point - log
mod.dOpen2.struc <- glmmTMB(Case ~ dOpen.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen2.struc$parameters$theta[1] = log(1e3)
mod.dOpen2.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[2]] <- glmmTMB:::fitTMB(mod.dOpen2.struc)

# point - quad
mod.dOpen3.struc <- glmmTMB(Case ~ dOpen + dOpen.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen3.struc$parameters$theta[1] = log(1e3)
mod.dOpen3.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[3]] <- glmmTMB:::fitTMB(mod.dOpen3.struc)

# 250 - linear
mod.dOpen4.struc <- glmmTMB(Case ~ open.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen4.struc$parameters$theta[1] = log(1e3)
mod.dOpen4.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[4]] <- glmmTMB:::fitTMB(mod.dOpen4.struc)

# 250 - log
mod.dOpen5.struc <- glmmTMB(Case ~ open.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen5.struc$parameters$theta[1] = log(1e3)
mod.dOpen5.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[5]] <- glmmTMB:::fitTMB(mod.dOpen5.struc)

# 250 - quad
mod.dOpen6.struc <- glmmTMB(Case ~ open.250 + open.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen6.struc$parameters$theta[1] = log(1e3)
mod.dOpen6.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[6]] <- glmmTMB:::fitTMB(mod.dOpen6.struc)

# 500 - linear
mod.dOpen7.struc <- glmmTMB(Case ~ open.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen7.struc$parameters$theta[1] = log(1e3)
mod.dOpen7.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[7]] <- glmmTMB:::fitTMB(mod.dOpen7.struc)

# 500 - log
mod.dOpen8.struc <- glmmTMB(Case ~ open.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen8.struc$parameters$theta[1] = log(1e3)
mod.dOpen8.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[8]] <- glmmTMB:::fitTMB(mod.dOpen8.struc)

# 500 - quad
mod.dOpen9.struc <- glmmTMB(Case ~ open.500 + open.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen9.struc$parameters$theta[1] = log(1e3)
mod.dOpen9.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[9]] <- glmmTMB:::fitTMB(mod.dOpen9.struc)

# 1000 - linear
mod.dOpen10.struc <- glmmTMB(Case ~ open.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen10.struc$parameters$theta[1] = log(1e3)
mod.dOpen10.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[10]] <- glmmTMB:::fitTMB(mod.dOpen10.struc)

# 1000 - log
mod.dOpen11.struc <- glmmTMB(Case ~ open.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen11.struc$parameters$theta[1] = log(1e3)
mod.dOpen11.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[11]] <- glmmTMB:::fitTMB(mod.dOpen11.struc)

# 1000 - quad
mod.dOpen12.struc <- glmmTMB(Case ~ open.1000 + open.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen12.struc$parameters$theta[1] = log(1e3)
mod.dOpen12.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[12]] <- glmmTMB:::fitTMB(mod.dOpen12.struc)

# null
mod.dOpen13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dOpen13.struc$parameters$theta[1] = log(1e3)
mod.dOpen13.struc$mapArg = list(theta = factor(c(NA)))

mod.dOpen[[13]] <- glmmTMB:::fitTMB(mod.dOpen13.struc)

aictab(mod.dOpen, modnames = modnames.1)

write.table(aictab(mod.dOpen, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.dOpen[[1]])

#_____________________________________________________________________________________________________________
# 5d. dMatureForest ----
#_____________________________________________________________________________________________________________

mod.dMatureForest <- list()

# point - linear
mod.dMatureForest1.struc <- glmmTMB(Case ~ dMatureForest + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest1.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest1.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[1]] <- glmmTMB:::fitTMB(mod.dMatureForest1.struc)

# point - log
mod.dMatureForest2.struc <- glmmTMB(Case ~ dMatureForest.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest2.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest2.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[2]] <- glmmTMB:::fitTMB(mod.dMatureForest2.struc)

# point - quad
mod.dMatureForest3.struc <- glmmTMB(Case ~ dMatureForest + dMatureForest.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest3.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest3.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[3]] <- glmmTMB:::fitTMB(mod.dMatureForest3.struc)

# 250 - linear
mod.dMatureForest4.struc <- glmmTMB(Case ~ mf.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest4.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest4.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[4]] <- glmmTMB:::fitTMB(mod.dMatureForest4.struc)

# 250 - log
mod.dMatureForest5.struc <- glmmTMB(Case ~ mf.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest5.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest5.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[5]] <- glmmTMB:::fitTMB(mod.dMatureForest5.struc)

# 250 - quad
mod.dMatureForest6.struc <- glmmTMB(Case ~ mf.250 + mf.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest6.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest6.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[6]] <- glmmTMB:::fitTMB(mod.dMatureForest6.struc)

# 500 - linear
mod.dMatureForest7.struc <- glmmTMB(Case ~ mf.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest7.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest7.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[7]] <- glmmTMB:::fitTMB(mod.dMatureForest7.struc)

# 500 - log
mod.dMatureForest8.struc <- glmmTMB(Case ~ mf.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest8.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest8.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[8]] <- glmmTMB:::fitTMB(mod.dMatureForest8.struc)

# 500 - quad
mod.dMatureForest9.struc <- glmmTMB(Case ~ mf.500 + mf.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest9.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest9.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[9]] <- glmmTMB:::fitTMB(mod.dMatureForest9.struc)

# 1000 - linear
mod.dMatureForest10.struc <- glmmTMB(Case ~ mf.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest10.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest10.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[10]] <- glmmTMB:::fitTMB(mod.dMatureForest10.struc)

# 1000 - log
mod.dMatureForest11.struc <- glmmTMB(Case ~ mf.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest11.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest11.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[11]] <- glmmTMB:::fitTMB(mod.dMatureForest11.struc)

# 1000 - quad
mod.dMatureForest12.struc <- glmmTMB(Case ~ mf.1000 + mf.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest12.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest12.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[12]] <- glmmTMB:::fitTMB(mod.dMatureForest12.struc)

# null
mod.dMatureForest13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.dMatureForest13.struc$parameters$theta[1] = log(1e3)
mod.dMatureForest13.struc$mapArg = list(theta = factor(c(NA)))

mod.dMatureForest[[13]] <- glmmTMB:::fitTMB(mod.dMatureForest13.struc)

aictab(mod.dMatureForest, modnames = modnames.1)

write.table(aictab(mod.dMatureForest, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.dMatureForest[[9]])
summary(mod.dMatureForest[[7]])

#_____________________________________________________________________________________________________________
# 5f. ED ----
#_____________________________________________________________________________________________________________

mod.ED <- list()

# 250 - linear
mod.ED4.struc <- glmmTMB(Case ~ ED.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED4.struc$parameters$theta[1] = log(1e3)
mod.ED4.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[1]] <- glmmTMB:::fitTMB(mod.ED4.struc)

# 250 - log
mod.ED5.struc <- glmmTMB(Case ~ ED.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED5.struc$parameters$theta[1] = log(1e3)
mod.ED5.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[2]] <- glmmTMB:::fitTMB(mod.ED5.struc)

# 250 - quad
mod.ED6.struc <- glmmTMB(Case ~ ED.250 + ED.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED6.struc$parameters$theta[1] = log(1e3)
mod.ED6.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[3]] <- glmmTMB:::fitTMB(mod.ED6.struc)

# 500 - linear
mod.ED7.struc <- glmmTMB(Case ~ ED.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED7.struc$parameters$theta[1] = log(1e3)
mod.ED7.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[4]] <- glmmTMB:::fitTMB(mod.ED7.struc)

# 500 - log
mod.ED8.struc <- glmmTMB(Case ~ ED.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED8.struc$parameters$theta[1] = log(1e3)
mod.ED8.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[5]] <- glmmTMB:::fitTMB(mod.ED8.struc)

# 500 - quad
mod.ED9.struc <- glmmTMB(Case ~ ED.500 + ED.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED9.struc$parameters$theta[1] = log(1e3)
mod.ED9.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[6]] <- glmmTMB:::fitTMB(mod.ED9.struc)

# 1000 - linear
mod.ED10.struc <- glmmTMB(Case ~ ED.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED10.struc$parameters$theta[1] = log(1e3)
mod.ED10.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[7]] <- glmmTMB:::fitTMB(mod.ED10.struc)

# 1000 - log
mod.ED11.struc <- glmmTMB(Case ~ ED.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED11.struc$parameters$theta[1] = log(1e3)
mod.ED11.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[8]] <- glmmTMB:::fitTMB(mod.ED11.struc)

# 1000 - quad
mod.ED12.struc <- glmmTMB(Case ~ ED.1000 + ED.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED12.struc$parameters$theta[1] = log(1e3)
mod.ED12.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[9]] <- glmmTMB:::fitTMB(mod.ED12.struc)

# null
mod.ED13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.ED13.struc$parameters$theta[1] = log(1e3)
mod.ED13.struc$mapArg = list(theta = factor(c(NA)))

mod.ED[[10]] <- glmmTMB:::fitTMB(mod.ED13.struc)

aictab(mod.ED, modnames = modnames.2)

write.table(aictab(mod.ED, modnames = modnames.2), "clipboard", sep = "\t")

summary(mod.ED[[4]])

#_____________________________________________________________________________________________________________
# 5g. TPI ----
#_____________________________________________________________________________________________________________

mod.TPI <- list()

# point - linear
mod.TPI1.struc <- glmmTMB(Case ~ TPI + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI1.struc$parameters$theta[1] = log(1e3)
mod.TPI1.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[1]] <- glmmTMB:::fitTMB(mod.TPI1.struc)

# point - log
mod.TPI2.struc <- glmmTMB(Case ~ TPI.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI2.struc$parameters$theta[1] = log(1e3)
mod.TPI2.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[2]] <- glmmTMB:::fitTMB(mod.TPI2.struc)

# point - quad
mod.TPI3.struc <- glmmTMB(Case ~ TPI + TPI.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI3.struc$parameters$theta[1] = log(1e3)
mod.TPI3.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[3]] <- glmmTMB:::fitTMB(mod.TPI3.struc)

# 250 - linear
mod.TPI4.struc <- glmmTMB(Case ~ TPI.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI4.struc$parameters$theta[1] = log(1e3)
mod.TPI4.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[4]] <- glmmTMB:::fitTMB(mod.TPI4.struc)

# 250 - log
mod.TPI5.struc <- glmmTMB(Case ~ TPI.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI5.struc$parameters$theta[1] = log(1e3)
mod.TPI5.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[5]] <- glmmTMB:::fitTMB(mod.TPI5.struc)

# 250 - quad
mod.TPI6.struc <- glmmTMB(Case ~ TPI.250 + TPI.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI6.struc$parameters$theta[1] = log(1e3)
mod.TPI6.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[6]] <- glmmTMB:::fitTMB(mod.TPI6.struc)

# 500 - linear
mod.TPI7.struc <- glmmTMB(Case ~ TPI.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI7.struc$parameters$theta[1] = log(1e3)
mod.TPI7.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[7]] <- glmmTMB:::fitTMB(mod.TPI7.struc)

# 500 - log
mod.TPI8.struc <- glmmTMB(Case ~ TPI.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI8.struc$parameters$theta[1] = log(1e3)
mod.TPI8.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[8]] <- glmmTMB:::fitTMB(mod.TPI8.struc)

# 500 - quad
mod.TPI9.struc <- glmmTMB(Case ~ TPI.500 + TPI.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI9.struc$parameters$theta[1] = log(1e3)
mod.TPI9.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[9]] <- glmmTMB:::fitTMB(mod.TPI9.struc)

# 1000 - linear
mod.TPI10.struc <- glmmTMB(Case ~ TPI.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI10.struc$parameters$theta[1] = log(1e3)
mod.TPI10.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[10]] <- glmmTMB:::fitTMB(mod.TPI10.struc)

# 1000 - log
mod.TPI11.struc <- glmmTMB(Case ~ TPI.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI11.struc$parameters$theta[1] = log(1e3)
mod.TPI11.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[11]] <- glmmTMB:::fitTMB(mod.TPI11.struc)

# 1000 - quad
mod.TPI12.struc <- glmmTMB(Case ~ TPI.1000 + TPI.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI12.struc$parameters$theta[1] = log(1e3)
mod.TPI12.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[12]] <- glmmTMB:::fitTMB(mod.TPI12.struc)

# null
mod.TPI13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TPI13.struc$parameters$theta[1] = log(1e3)
mod.TPI13.struc$mapArg = list(theta = factor(c(NA)))

mod.TPI[[13]] <- glmmTMB:::fitTMB(mod.TPI13.struc)

aictab(mod.TPI, modnames = modnames.1)

write.table(aictab(mod.TPI, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.TPI[[4]])

#_____________________________________________________________________________________________________________
# 5h. TRI ----
#_____________________________________________________________________________________________________________

mod.TRI <- list()

# point - linear
mod.TRI1.struc <- glmmTMB(Case ~ TRI + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI1.struc$parameters$theta[1] = log(1e3)
mod.TRI1.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[1]] <- glmmTMB:::fitTMB(mod.TRI1.struc)

# point - log
mod.TRI2.struc <- glmmTMB(Case ~ TRI.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI2.struc$parameters$theta[1] = log(1e3)
mod.TRI2.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[2]] <- glmmTMB:::fitTMB(mod.TRI2.struc)

# point - quad
mod.TRI3.struc <- glmmTMB(Case ~ TRI + TRI.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI3.struc$parameters$theta[1] = log(1e3)
mod.TRI3.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[3]] <- glmmTMB:::fitTMB(mod.TRI3.struc)

# 250 - linear
mod.TRI4.struc <- glmmTMB(Case ~ TRI.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI4.struc$parameters$theta[1] = log(1e3)
mod.TRI4.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[4]] <- glmmTMB:::fitTMB(mod.TRI4.struc)

# 250 - log
mod.TRI5.struc <- glmmTMB(Case ~ TRI.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI5.struc$parameters$theta[1] = log(1e3)
mod.TRI5.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[5]] <- glmmTMB:::fitTMB(mod.TRI5.struc)

# 250 - quad
mod.TRI6.struc <- glmmTMB(Case ~ TRI.250 + TRI.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI6.struc$parameters$theta[1] = log(1e3)
mod.TRI6.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[6]] <- glmmTMB:::fitTMB(mod.TRI6.struc)

# 500 - linear
mod.TRI7.struc <- glmmTMB(Case ~ TRI.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI7.struc$parameters$theta[1] = log(1e3)
mod.TRI7.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[7]] <- glmmTMB:::fitTMB(mod.TRI7.struc)

# 500 - log
mod.TRI8.struc <- glmmTMB(Case ~ TRI.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI8.struc$parameters$theta[1] = log(1e3)
mod.TRI8.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[8]] <- glmmTMB:::fitTMB(mod.TRI8.struc)

# 500 - quad
mod.TRI9.struc <- glmmTMB(Case ~ TRI.500 + TRI.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI9.struc$parameters$theta[1] = log(1e3)
mod.TRI9.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[9]] <- glmmTMB:::fitTMB(mod.TRI9.struc)

# 1000 - linear
mod.TRI10.struc <- glmmTMB(Case ~ TRI.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI10.struc$parameters$theta[1] = log(1e3)
mod.TRI10.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[10]] <- glmmTMB:::fitTMB(mod.TRI10.struc)

# 1000 - log
mod.TRI11.struc <- glmmTMB(Case ~ TRI.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI11.struc$parameters$theta[1] = log(1e3)
mod.TRI11.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[11]] <- glmmTMB:::fitTMB(mod.TRI11.struc)

# 1000 - quad
mod.TRI12.struc <- glmmTMB(Case ~ TRI.1000 + TRI.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI12.struc$parameters$theta[1] = log(1e3)
mod.TRI12.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[12]] <- glmmTMB:::fitTMB(mod.TRI12.struc)

# null
mod.TRI13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.TRI13.struc$parameters$theta[1] = log(1e3)
mod.TRI13.struc$mapArg = list(theta = factor(c(NA)))

mod.TRI[[13]] <- glmmTMB:::fitTMB(mod.TRI13.struc)

aictab(mod.TRI, modnames = modnames.1)

write.table(aictab(mod.TRI, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.TRI[[1]])

#_____________________________________________________________________________________________________________
# 5i. SRI ----
#_____________________________________________________________________________________________________________

mod.SRI <- list()

# point - linear
mod.SRI1.struc <- glmmTMB(Case ~ SRI + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI1.struc$parameters$theta[1] = log(1e3)
mod.SRI1.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[1]] <- glmmTMB:::fitTMB(mod.SRI1.struc)

# point - log
mod.SRI2.struc <- glmmTMB(Case ~ SRI.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI2.struc$parameters$theta[1] = log(1e3)
mod.SRI2.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[2]] <- glmmTMB:::fitTMB(mod.SRI2.struc)

# point - quad
mod.SRI3.struc <- glmmTMB(Case ~ SRI + SRI.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI3.struc$parameters$theta[1] = log(1e3)
mod.SRI3.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[3]] <- glmmTMB:::fitTMB(mod.SRI3.struc)

# 250 - linear
mod.SRI4.struc <- glmmTMB(Case ~ SRI.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI4.struc$parameters$theta[1] = log(1e3)
mod.SRI4.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[4]] <- glmmTMB:::fitTMB(mod.SRI4.struc)

# 250 - log
mod.SRI5.struc <- glmmTMB(Case ~ SRI.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI5.struc$parameters$theta[1] = log(1e3)
mod.SRI5.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[5]] <- glmmTMB:::fitTMB(mod.SRI5.struc)

# 250 - quad
mod.SRI6.struc <- glmmTMB(Case ~ SRI.250 + SRI.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI6.struc$parameters$theta[1] = log(1e3)
mod.SRI6.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[6]] <- glmmTMB:::fitTMB(mod.SRI6.struc)

# 500 - linear
mod.SRI7.struc <- glmmTMB(Case ~ SRI.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI7.struc$parameters$theta[1] = log(1e3)
mod.SRI7.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[7]] <- glmmTMB:::fitTMB(mod.SRI7.struc)

# 500 - log
mod.SRI8.struc <- glmmTMB(Case ~ SRI.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI8.struc$parameters$theta[1] = log(1e3)
mod.SRI8.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[8]] <- glmmTMB:::fitTMB(mod.SRI8.struc)

# 500 - quad
mod.SRI9.struc <- glmmTMB(Case ~ SRI.500 + SRI.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI9.struc$parameters$theta[1] = log(1e3)
mod.SRI9.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[9]] <- glmmTMB:::fitTMB(mod.SRI9.struc)

# 1000 - linear
mod.SRI10.struc <- glmmTMB(Case ~ SRI.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI10.struc$parameters$theta[1] = log(1e3)
mod.SRI10.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[10]] <- glmmTMB:::fitTMB(mod.SRI10.struc)

# 1000 - log
mod.SRI11.struc <- glmmTMB(Case ~ SRI.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI11.struc$parameters$theta[1] = log(1e3)
mod.SRI11.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[11]] <- glmmTMB:::fitTMB(mod.SRI11.struc)

# 1000 - quad
mod.SRI12.struc <- glmmTMB(Case ~ SRI.1000 + SRI.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI12.struc$parameters$theta[1] = log(1e3)
mod.SRI12.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[12]] <- glmmTMB:::fitTMB(mod.SRI12.struc)

# null
mod.SRI13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.SRI13.struc$parameters$theta[1] = log(1e3)
mod.SRI13.struc$mapArg = list(theta = factor(c(NA)))

mod.SRI[[13]] <- glmmTMB:::fitTMB(mod.SRI13.struc)

aictab(mod.SRI, modnames = modnames.1)

write.table(aictab(mod.SRI, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.SRI[[8]])

#_____________________________________________________________________________________________________________
# 5j. canopy ----
#_____________________________________________________________________________________________________________

mod.canopy <- list()

# point - linear
mod.canopy1.struc <- glmmTMB(Case ~ canopy + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy1.struc$parameters$theta[1] = log(1e3)
mod.canopy1.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[1]] <- glmmTMB:::fitTMB(mod.canopy1.struc)

# point - log
mod.canopy2.struc <- glmmTMB(Case ~ canopy.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy2.struc$parameters$theta[1] = log(1e3)
mod.canopy2.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[2]] <- glmmTMB:::fitTMB(mod.canopy2.struc)

# point - quad
mod.canopy3.struc <- glmmTMB(Case ~ canopy + canopy.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy3.struc$parameters$theta[1] = log(1e3)
mod.canopy3.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[3]] <- glmmTMB:::fitTMB(mod.canopy3.struc)

# 250 - linear
mod.canopy4.struc <- glmmTMB(Case ~ canopy.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy4.struc$parameters$theta[1] = log(1e3)
mod.canopy4.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[4]] <- glmmTMB:::fitTMB(mod.canopy4.struc)

# 250 - log
mod.canopy5.struc <- glmmTMB(Case ~ canopy.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy5.struc$parameters$theta[1] = log(1e3)
mod.canopy5.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[5]] <- glmmTMB:::fitTMB(mod.canopy5.struc)

# 250 - quad
mod.canopy6.struc <- glmmTMB(Case ~ canopy.250 + canopy.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy6.struc$parameters$theta[1] = log(1e3)
mod.canopy6.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[6]] <- glmmTMB:::fitTMB(mod.canopy6.struc)

# 500 - linear
mod.canopy7.struc <- glmmTMB(Case ~ canopy.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy7.struc$parameters$theta[1] = log(1e3)
mod.canopy7.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[7]] <- glmmTMB:::fitTMB(mod.canopy7.struc)

# 500 - log
mod.canopy8.struc <- glmmTMB(Case ~ canopy.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy8.struc$parameters$theta[1] = log(1e3)
mod.canopy8.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[8]] <- glmmTMB:::fitTMB(mod.canopy8.struc)

# 500 - quad
mod.canopy9.struc <- glmmTMB(Case ~ canopy.500 + canopy.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy9.struc$parameters$theta[1] = log(1e3)
mod.canopy9.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[9]] <- glmmTMB:::fitTMB(mod.canopy9.struc)

# 1000 - linear
mod.canopy10.struc <- glmmTMB(Case ~ canopy.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy10.struc$parameters$theta[1] = log(1e3)
mod.canopy10.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[10]] <- glmmTMB:::fitTMB(mod.canopy10.struc)

# 1000 - log
mod.canopy11.struc <- glmmTMB(Case ~ canopy.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy11.struc$parameters$theta[1] = log(1e3)
mod.canopy11.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[11]] <- glmmTMB:::fitTMB(mod.canopy11.struc)

# 1000 - quad
mod.canopy12.struc <- glmmTMB(Case ~ canopy.1000 + canopy.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy12.struc$parameters$theta[1] = log(1e3)
mod.canopy12.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[12]] <- glmmTMB:::fitTMB(mod.canopy12.struc)

# null
mod.canopy13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.canopy13.struc$parameters$theta[1] = log(1e3)
mod.canopy13.struc$mapArg = list(theta = factor(c(NA)))

mod.canopy[[13]] <- glmmTMB:::fitTMB(mod.canopy13.struc)

aictab(mod.canopy, modnames = modnames.1)

write.table(aictab(mod.canopy, modnames = modnames.1), "clipboard", sep = "\t")

summary(mod.canopy[[3]])

#_____________________________________________________________________________________________________________
# 5k. IJI ----
#_____________________________________________________________________________________________________________

mod.IJI <- list()

# 250 - linear
mod.IJI4.struc <- glmmTMB(Case ~ IJI.250 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI4.struc$parameters$theta[1] = log(1e3)
mod.IJI4.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[1]] <- glmmTMB:::fitTMB(mod.IJI4.struc)

# 250 - log
mod.IJI5.struc <- glmmTMB(Case ~ IJI.250.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI5.struc$parameters$theta[1] = log(1e3)
mod.IJI5.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[2]] <- glmmTMB:::fitTMB(mod.IJI5.struc)

# 250 - quad
mod.IJI6.struc <- glmmTMB(Case ~ IJI.250 + IJI.250.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI6.struc$parameters$theta[1] = log(1e3)
mod.IJI6.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[3]] <- glmmTMB:::fitTMB(mod.IJI6.struc)

# 500 - linear
mod.IJI7.struc <- glmmTMB(Case ~ IJI.500 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI7.struc$parameters$theta[1] = log(1e3)
mod.IJI7.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[4]] <- glmmTMB:::fitTMB(mod.IJI7.struc)

# 500 - log
mod.IJI8.struc <- glmmTMB(Case ~ IJI.500.P +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI8.struc$parameters$theta[1] = log(1e3)
mod.IJI8.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[5]] <- glmmTMB:::fitTMB(mod.IJI8.struc)

# 500 - quad
mod.IJI9.struc <- glmmTMB(Case ~ IJI.500 + IJI.500.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI9.struc$parameters$theta[1] = log(1e3)
mod.IJI9.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[6]] <- glmmTMB:::fitTMB(mod.IJI9.struc)

# 1000 - linear
mod.IJI10.struc <- glmmTMB(Case ~ IJI.1000 + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI10.struc$parameters$theta[1] = log(1e3)
mod.IJI10.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[7]] <- glmmTMB:::fitTMB(mod.IJI10.struc)

# 1000 - log
mod.IJI11.struc <- glmmTMB(Case ~ IJI.1000.P + 
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI11.struc$parameters$theta[1] = log(1e3)
mod.IJI11.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[8]] <- glmmTMB:::fitTMB(mod.IJI11.struc)

# 1000 - quad
mod.IJI12.struc <- glmmTMB(Case ~ IJI.1000 + IJI.1000.Q +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI12.struc$parameters$theta[1] = log(1e3)
mod.IJI12.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[9]] <- glmmTMB:::fitTMB(mod.IJI12.struc)

# null
mod.IJI13.struc <- glmmTMB(Case ~ 1 +
                                 (1 | Animal),
                        family = poisson,
                        data = scaled.data,
                        doFit = FALSE)

mod.IJI13.struc$parameters$theta[1] = log(1e3)
mod.IJI13.struc$mapArg = list(theta = factor(c(NA)))

mod.IJI[[10]] <- glmmTMB:::fitTMB(mod.IJI13.struc)

aictab(mod.IJI, modnames = modnames.2)

write.table(aictab(mod.IJI, modnames = modnames.2), "clipboard", sep = "\t")

summary(mod.IJI[[3]])

check_collinearity(mod.IJI[[3]])

#_____________________________________________________________________________________________________________
# 5m. EVI ----
#_____________________________________________________________________________________________________________

mod.EVI <- list()

# 250 - linear
mod.EVI4.struc <- glmmTMB(Case ~ EVI + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI4.struc$parameters$theta[1] = log(1e3)
mod.EVI4.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[1]] <- glmmTMB:::fitTMB(mod.EVI4.struc)

# 250 - log
mod.EVI5.struc <- glmmTMB(Case ~ EVI.P + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI5.struc$parameters$theta[1] = log(1e3)
mod.EVI5.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[2]] <- glmmTMB:::fitTMB(mod.EVI5.struc)

# 250 - quad
mod.EVI6.struc <- glmmTMB(Case ~ EVI + EVI.Q +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI6.struc$parameters$theta[1] = log(1e3)
mod.EVI6.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[3]] <- glmmTMB:::fitTMB(mod.EVI6.struc)

# 500 - linear
mod.EVI7.struc <- glmmTMB(Case ~ EVI.1pix +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI7.struc$parameters$theta[1] = log(1e3)
mod.EVI7.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[4]] <- glmmTMB:::fitTMB(mod.EVI7.struc)

# 500 - log
mod.EVI8.struc <- glmmTMB(Case ~ EVI.1pix.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI8.struc$parameters$theta[1] = log(1e3)
mod.EVI8.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[5]] <- glmmTMB:::fitTMB(mod.EVI8.struc)

# 500 - quad
mod.EVI9.struc <- glmmTMB(Case ~ EVI.1pix + EVI.1pix.Q +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI9.struc$parameters$theta[1] = log(1e3)
mod.EVI9.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[6]] <- glmmTMB:::fitTMB(mod.EVI9.struc)

# 1000 - linear
mod.EVI10.struc <- glmmTMB(Case ~ EVI.2pix + 
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.EVI10.struc$parameters$theta[1] = log(1e3)
mod.EVI10.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[7]] <- glmmTMB:::fitTMB(mod.EVI10.struc)

# 1000 - log
mod.EVI11.struc <- glmmTMB(Case ~ EVI.2pix.P + 
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.EVI11.struc$parameters$theta[1] = log(1e3)
mod.EVI11.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[8]] <- glmmTMB:::fitTMB(mod.EVI11.struc)

# 1000 - quad
mod.EVI12.struc <- glmmTMB(Case ~ EVI.2pix + EVI.2pix.Q +
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.EVI12.struc$parameters$theta[1] = log(1e3)
mod.EVI12.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[9]] <- glmmTMB:::fitTMB(mod.EVI12.struc)

# null
mod.EVI13.struc <- glmmTMB(Case ~1 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.EVI13.struc$parameters$theta[1] = log(1e3)
mod.EVI13.struc$mapArg = list(theta = factor(c(NA)))

mod.EVI[[10]] <- glmmTMB:::fitTMB(mod.EVI13.struc)

aictab(mod.EVI, modnames = modnames.2)

write.table(aictab(mod.EVI, modnames = modnames.2), "clipboard", sep = "\t")

summary(mod.EVI[[7]])

#_____________________________________________________________________________________________________________
# 5o. ai.mf ----
#_____________________________________________________________________________________________________________

mod.ai.mf <- list()

# 250 - linear
mod.ai.mf4.struc <- glmmTMB(Case ~ ai.mf.250 + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf4.struc$parameters$theta[1] = log(1e3)
mod.ai.mf4.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[1]] <- glmmTMB:::fitTMB(mod.ai.mf4.struc)

# 250 - log
mod.ai.mf5.struc <- glmmTMB(Case ~ ai.mf.250.P + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf5.struc$parameters$theta[1] = log(1e3)
mod.ai.mf5.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[2]] <- glmmTMB:::fitTMB(mod.ai.mf5.struc)

# 250 - quad
mod.ai.mf6.struc <- glmmTMB(Case ~ ai.mf.250 + ai.mf.250.Q +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf6.struc$parameters$theta[1] = log(1e3)
mod.ai.mf6.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[3]] <- glmmTMB:::fitTMB(mod.ai.mf6.struc)

# 500 - linear
mod.ai.mf7.struc <- glmmTMB(Case ~ ai.mf.500 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf7.struc$parameters$theta[1] = log(1e3)
mod.ai.mf7.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[4]] <- glmmTMB:::fitTMB(mod.ai.mf7.struc)

# 500 - log
mod.ai.mf8.struc <- glmmTMB(Case ~ ai.mf.500.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf8.struc$parameters$theta[1] = log(1e3)
mod.ai.mf8.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[5]] <- glmmTMB:::fitTMB(mod.ai.mf8.struc)

# 500 - quad
mod.ai.mf9.struc <- glmmTMB(Case ~ ai.mf.500 + ai.mf.500.Q +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.ai.mf9.struc$parameters$theta[1] = log(1e3)
mod.ai.mf9.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[6]] <- glmmTMB:::fitTMB(mod.ai.mf9.struc)

# 1000 - linear
mod.ai.mf10.struc <- glmmTMB(Case ~ ai.mf.1000 + 
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.ai.mf10.struc$parameters$theta[1] = log(1e3)
mod.ai.mf10.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[7]] <- glmmTMB:::fitTMB(mod.ai.mf10.struc)

# 1000 - log
mod.ai.mf11.struc <- glmmTMB(Case ~ ai.mf.1000.P + 
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.ai.mf11.struc$parameters$theta[1] = log(1e3)
mod.ai.mf11.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[8]] <- glmmTMB:::fitTMB(mod.ai.mf11.struc)

# 1000 - quad
mod.ai.mf12.struc <- glmmTMB(Case ~ ai.mf.1000 + ai.mf.1000.Q +
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.ai.mf12.struc$parameters$theta[1] = log(1e3)
mod.ai.mf12.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[9]] <- glmmTMB:::fitTMB(mod.ai.mf12.struc)

# null
mod.ai.mf13.struc <- glmmTMB(Case ~ 1 +
                              (1 | Animal),
                            family = poisson,
                            data = scaled.data,
                            doFit = FALSE)

mod.ai.mf13.struc$parameters$theta[1] = log(1e3)
mod.ai.mf13.struc$mapArg = list(theta = factor(c(NA)))

mod.ai.mf[[10]] <- glmmTMB:::fitTMB(mod.ai.mf13.struc)

aictab(mod.ai.mf, modnames = modnames.2)

write.table(aictab(mod.ai.mf, modnames = modnames.2), "clipboard", sep = "\t")

summary(mod.ai.mf[[4]])

#_____________________________________________________________________________________________________________
# 6. Density plots ----
#_____________________________________________________________________________________________________________

# dRoad
ggplot(data = all.data, aes(x = dRoad.1000)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# Open
ggplot(data = all.data, aes(x = dOpen)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5) +
  coord_cartesian(xlim = c(0, 1500))

# mf.500
ggplot(data = all.data, aes(x = mf.500)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# ED.500
ggplot(data = all.data, aes(x = ED.500)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# TPI
ggplot(data = all.data, aes(x = TPI.250)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# TRI
ggplot(data = all.data, aes(x = TRI)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# SRI.500.P
ggplot(data = all.data, aes(x = SRI.500.P)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# Canopy
ggplot(data = all.data, aes(x = canopy)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# IJI.250
ggplot(data = all.data, aes(x = IJI.250)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# EVI.1000
ggplot(data = all.data, aes(x = EVI.2pix)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

# ai.mf.500
ggplot(data = all.data, aes(x = ai.mf.500)) +
       theme_bw() +
       geom_density(aes(fill = as.factor(Case), color = as.factor(Case)),
                    alpha = 0.5)

#_____________________________________________________________________________________________________________
# 7. Correlation ----
#_____________________________________________________________________________________________________________

# Pearson's
cor.table.2 <- cor(x = scaled.data[, c("dRoad.1000", "dOpen", "mf.500",
                                       "ED.500", "TPI.250", "TRI", "SRI.500.P", "canopy", "canopy.Q", 
                                       "IJI.250", "IJI.250.Q", 
                                       "EVI.2pix", "ai.mf.500")],
                  method = "pearson")

write.table(cor.table.2, "clipboard", sep = "\t")

# Spearman's
cor.table.3 <- cor(x = scaled.data[, c("dRoad.1000", "dOpen", "mf.500",
                                       "ED.500", "TPI.250", "TRI", "SRI.500.P", "canopy", "canopy.Q", 
                                       "IJI.250", "IJI.250.Q", 
                                       "EVI.2pix", "ai.mf.500")],
                   method = "spearman")

write.table(cor.table.3, "clipboard", sep = "\t")

#_____________________________________________________________________________________________________________
# 8. Which correlated variables to keep? ----
#_____________________________________________________________________________________________________________

mod.corr <- list()

# mf.500
mod.corr1.struc <- glmmTMB(Case ~ mf.500 + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.corr1.struc$parameters$theta[1] = log(1e3)
mod.corr1.struc$mapArg = list(theta = factor(c(NA)))

mod.corr[[1]] <- glmmTMB:::fitTMB(mod.corr1.struc)

# ED.500
mod.corr2.struc <- glmmTMB(Case ~ ED.500 + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.corr2.struc$parameters$theta[1] = log(1e3)
mod.corr2.struc$mapArg = list(theta = factor(c(NA)))

mod.corr[[2]] <- glmmTMB:::fitTMB(mod.corr2.struc)

# ai.mf.500
mod.corr3.struc <- glmmTMB(Case ~ ai.mf.500 + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.corr3.struc$parameters$theta[1] = log(1e3)
mod.corr3.struc$mapArg = list(theta = factor(c(NA)))

mod.corr[[3]] <- glmmTMB:::fitTMB(mod.corr3.struc)

# dOpen
mod.corr4.struc <- glmmTMB(Case ~ dOpen + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.corr4.struc$parameters$theta[1] = log(1e3)
mod.corr4.struc$mapArg = list(theta = factor(c(NA)))

mod.corr[[4]] <- glmmTMB:::fitTMB(mod.corr4.struc)

# canopy
mod.corr5.struc <- glmmTMB(Case ~ canopy + canopy.Q + 
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

mod.corr5.struc$parameters$theta[1] = log(1e3)
mod.corr5.struc$mapArg = list(theta = factor(c(NA)))

mod.corr[[5]] <- glmmTMB:::fitTMB(mod.corr5.struc)

aictab(mod.corr, modnames = c("mf", "ed", "ai", "open", "canopy"))

summary(mod.corr[[3]])

#_____________________________________________________________________________________________________________
# 9. Final model selection - basic groupings ----
#_____________________________________________________________________________________________________________

final.models <- list()

# 1. FORCOV
final.models1.struc <- glmmTMB(Case ~ EVI.2pix +
                                      canopy + canopy.Q +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models1.struc$parameters$theta[1] = log(1e3)
final.models1.struc$mapArg = list(theta = factor(c(NA)))

final.models[[1]] <- glmmTMB:::fitTMB(final.models1.struc)

# 2. LAND
final.models2.struc <- glmmTMB(Case ~ IJI.250 + IJI.250.Q +
                                      ED.500 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models2.struc$parameters$theta[1] = log(1e3)
final.models2.struc$mapArg = list(theta = factor(c(NA)))

final.models[[2]] <- glmmTMB:::fitTMB(final.models2.struc)

# 3. TOPO
final.models3.struc <- glmmTMB(Case ~ TRI +
                                      TPI.250 +
                                      SRI.500.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models3.struc$parameters$theta[1] = log(1e3)
final.models3.struc$mapArg = list(theta = factor(c(NA)))

final.models[[3]] <- glmmTMB:::fitTMB(final.models3.struc)

# 4. HUMAN
final.models4.struc <- glmmTMB(Case ~ dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models4.struc$parameters$theta[1] = log(1e3)
final.models4.struc$mapArg = list(theta = factor(c(NA)))

final.models[[4]] <- glmmTMB:::fitTMB(final.models4.struc)

# 5. FORCOV + LAND
final.models5.struc <- glmmTMB(Case ~ EVI.2pix +
                                      canopy + canopy.Q +
                                      IJI.250 + IJI.250.Q +
                                      ED.500 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models5.struc$parameters$theta[1] = log(1e3)
final.models5.struc$mapArg = list(theta = factor(c(NA)))

final.models[[5]] <- glmmTMB:::fitTMB(final.models5.struc)

# 6. FORCOV + TOPO
final.models6.struc <- glmmTMB(Case ~ EVI.2pix +
                                      canopy + canopy.Q +
                                      TRI +
                                      TPI.250 +
                                      SRI.500.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models6.struc$parameters$theta[1] = log(1e3)
final.models6.struc$mapArg = list(theta = factor(c(NA)))

final.models[[6]] <- glmmTMB:::fitTMB(final.models6.struc)

# 7. FORCOV + HUMAN
final.models7.struc <- glmmTMB(Case ~ EVI.2pix +
                                      canopy + canopy.Q +
                                      dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models7.struc$parameters$theta[1] = log(1e3)
final.models7.struc$mapArg = list(theta = factor(c(NA)))

final.models[[7]] <- glmmTMB:::fitTMB(final.models7.struc)

# 8. LAND + TOPO
final.models8.struc <- glmmTMB(Case ~ IJI.250 + IJI.250.Q +
                                      ED.500 +
                                      TRI +
                                      TPI.250 +
                                      SRI.500.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models8.struc$parameters$theta[1] = log(1e3)
final.models8.struc$mapArg = list(theta = factor(c(NA)))

final.models[[8]] <- glmmTMB:::fitTMB(final.models8.struc)

# 9. LAND + HUMAN
final.models9.struc <- glmmTMB(Case ~ IJI.250 + IJI.250.Q +
                                      ED.500 +
                                      dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models9.struc$parameters$theta[1] = log(1e3)
final.models9.struc$mapArg = list(theta = factor(c(NA)))

final.models[[9]] <- glmmTMB:::fitTMB(final.models9.struc)

# 10. TOPO + HUMAN
final.models10.struc <- glmmTMB(Case ~ TRI +
                                       TPI.250 +
                                       SRI.500.P +
                                       dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models10.struc$parameters$theta[1] = log(1e3)
final.models10.struc$mapArg = list(theta = factor(c(NA)))

final.models[[10]] <- glmmTMB:::fitTMB(final.models10.struc)

# 11. FORCOV + LAND + TOPO
final.models11.struc <- glmmTMB(Case ~ EVI.2pix +
                                       canopy + canopy.Q +
                                       IJI.250 + IJI.250.Q +
                                       ED.500 +
                                       TRI +
                                       TPI.250 +
                                       SRI.500.P +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models11.struc$parameters$theta[1] = log(1e3)
final.models11.struc$mapArg = list(theta = factor(c(NA)))

final.models[[11]] <- glmmTMB:::fitTMB(final.models11.struc)

# 12. FORCOV + LAND + HUMAN
final.models12.struc <- glmmTMB(Case ~ EVI.2pix +
                                       canopy + canopy.Q +
                                       IJI.250 + IJI.250.Q +
                                       ED.500 +
                                       dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models12.struc$parameters$theta[1] = log(1e3)
final.models12.struc$mapArg = list(theta = factor(c(NA)))

final.models[[12]] <- glmmTMB:::fitTMB(final.models12.struc)

# 13. FORCOV + TOPO + HUMAN
final.models13.struc <- glmmTMB(Case ~ EVI.2pix +
                                       canopy + canopy.Q +
                                       TRI +
                                       TPI.250 +
                                       SRI.500.P +
                                       dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models13.struc$parameters$theta[1] = log(1e3)
final.models13.struc$mapArg = list(theta = factor(c(NA)))

final.models[[13]] <- glmmTMB:::fitTMB(final.models13.struc)

# 14. LAND + TOPO + HUMAN
final.models14.struc <- glmmTMB(Case ~ IJI.250 + IJI.250.Q +
                                       ED.500 +
                                       TRI +
                                       TPI.250 +
                                       SRI.500.P +
                                       dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models14.struc$parameters$theta[1] = log(1e3)
final.models14.struc$mapArg = list(theta = factor(c(NA)))

final.models[[14]] <- glmmTMB:::fitTMB(final.models14.struc)

# 15. Global
final.models15.struc <- glmmTMB(Case ~ EVI.2pix +
                                       canopy + canopy.Q +
                                       IJI.250 + IJI.250.Q +
                                       ED.500 +
                                       TRI +
                                       TPI.250 +
                                       SRI.500.P +
                                       dRoad.1000 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models15.struc$parameters$theta[1] = log(1e3)
final.models15.struc$mapArg = list(theta = factor(c(NA)))

final.models[[15]] <- glmmTMB:::fitTMB(final.models15.struc)

# 16. Null
final.models16.struc <- glmmTMB(Case ~ 1 +
                             (1 | Animal),
                           family = poisson,
                           data = scaled.data,
                           doFit = FALSE)

final.models16.struc$parameters$theta[1] = log(1e3)
final.models16.struc$mapArg = list(theta = factor(c(NA)))

final.models[[16]] <- glmmTMB:::fitTMB(final.models16.struc)

aictab(final.models, modnames = c("FORCOV", "LAND", "TOPO", "HUMAN",
                                  "FORCOV + LAND", "FORCOV + TOPO",
                                  "FORCOV + HUMAN", "LAND + TOPO", 
                                  "LAND + HUMAN", "TOPO + HUMAN", 
                                  "FORCOV + LAND + TOPO", "FORCOV + LAND + HUMAN",
                                  "FORCOV + TOPO + HUMAN", "LAND + TOPO + HUMAN",
                                  "Global", "Null"))

summary(final.models[[15]])
summary(final.models[[13]])

check_collinearity(final.models[[15]])
check_collinearity(final.models[[13]])

write.table(aictab(final.models,
                   modnames = c("FORCOV", "LAND", "TOPO", "HUMAN",
                                "FORCOV + LAND", "FORCOV + TOPO",
                                "FORCOV + HUMAN", "LAND + TOPO", 
                                "LAND + HUMAN", "TOPO + HUMAN", 
                                "FORCOV + LAND + TOPO", "FORCOV + LAND + HUMAN",
                                "FORCOV + TOPO + HUMAN", "LAND + TOPO + HUMAN",
                                "Global", "Null")), 
            "clipboard", sep = "\t")

#_____________________________________________________________________________________________________________
# 10. Save image ----
#_____________________________________________________________________________________________________________

save.image(file = "partsite_models.RData")

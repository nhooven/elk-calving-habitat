# Title: Parturition site selection
# Subtitle: 8 - Sample size and coefficient convergence (continous covariates only)
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 12 Oct 2021
# Date completed: 12 Oct 2021
# Date modified: 3 Aug 2022
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(lubridate)
library(sp)               # work with spatial objects
library(raster)           # work with rasters
library(adehabitatHR)     # fit MCPs
library(rgeos)            # buffer
library(landscapemetrics) # sample landscape metrics
library(mefa4)            # notin function
library(glmmTMB)

#_____________________________________________________________________________________________________________
# 2. Read in data and work with CRSs ----
#_____________________________________________________________________________________________________________

# import "sites" datasets
sites <- read.csv("sites.csv")
sites.2022 <- read.csv("sites_2022.csv")

# import relocation data
elk.data <- read.csv("D:/Elk project/Data analysis/Raw data processing/Relocations_vectronic_1.csv")
elk.data.2022 <- read.csv("Relocations_2022.csv")

# get dates into correct formats
elk.data$t <- as.POSIXct(mdy_hm(elk.data$t, tz = "America/New_York"))
elk.data.2022$t <- as.POSIXct(elk.data.2022$t, tz = "America/New_York")

# Raster directory
raster.dir <- "D:/Elk project/Elk Zone rasters (7-20-21)"

# Read in each raster
dDeveloped    <- raster(paste0(raster.dir, "/", "dDeveloped.tif"))

TRI           <- raster(paste0(raster.dir, "/", "TRI_10.tif"))

TPI           <- raster(paste0(raster.dir, "/", "TPI_10.tif"))

SRI           <- raster(paste0(raster.dir, "/", "SRI_10.tif"))

dRoad         <- raster(paste0(raster.dir, "/", "dRoad.tif"))

dEdge         <- raster(paste0(raster.dir, "/", "dEdge.tif"))

dYoungForest  <- raster(paste0(raster.dir, "/", "dYoungForest.tif"))

dMatureForest <- raster(paste0(raster.dir, "/", "dMatureForest.tif"))

dOpen         <- raster(paste0(raster.dir, "/", "dOpen.tif"))

canopy        <- raster(paste0(raster.dir, "/", "canopy.tif"))

# define raster projection
raster.proj <- dDeveloped@crs

# define lat long CRS
latlong.projection <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# define UTM CRS
UTM.projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")

#_____________________________________________________________________________________________________________
# 3. Convert sites to UTM ----
#_____________________________________________________________________________________________________________
# 3a. 2020 and 2021 sites ----
#_____________________________________________________________________________________________________________

# convert to SpatialPoints
sites.sp <- SpatialPoints(coords = sites[ ,c("Longitude", "Latitude")],
                          proj4string = latlong.projection)

# transform to UTM
sites.sp.UTM <- spTransform(sites.sp, UTM.projection)

# rename and bind to original dataframe and rename
UTM.sites <- as.data.frame(sites.sp.UTM@coords)
colnames(UTM.sites) <- c("x", "y")

sites <- cbind(sites, UTM.sites)

# keep only the columns we need, rename to match elk.data, add case = 1
sites.1 <- sites %>% dplyr::select(x, y, t = BirthDate, Animal = CollarID, Year) %>%
                     mutate(Case = 1)

#_____________________________________________________________________________________________________________
# 3b. 2022 sites ----
#_____________________________________________________________________________________________________________

# convert to SpatialPoints
sites.2022.sp <- SpatialPoints(coords = sites.2022[ ,c("Longitude", "Latitude")],
                               proj4string = latlong.projection)

# transform to UTM
sites.2022.sp.UTM <- spTransform(sites.2022.sp, UTM.projection)

# rename and bind to original dataframe and rename
UTM.sites.2022 <- as.data.frame(sites.2022.sp.UTM@coords)
colnames(UTM.sites.2022) <- c("x", "y")

sites.2022 <- cbind(sites.2022, UTM.sites.2022)

# keep only the columns we need, rename to match elk.data, add case = 1
sites.2022.1 <- sites.2022 %>% dplyr::select(x, y, t = BirthDate, Animal = CollarID, Year) %>%
                               mutate(Case = 1)

# bind together
sites.all <- rbind(sites.1, sites.2022.1)

#_____________________________________________________________________________________________________________
# 4. Sample 5000 available locations from each home range ----
#_____________________________________________________________________________________________________________

# bind data together
elk.data.2022 <- elk.data.2022 %>% dplyr::select(-X)

elk.data.all <- rbind(elk.data, elk.data.2022)

available.points <- data.frame()

# for loop to run through each animal's data
for (i in unique(sites.all$Animal)) {
  
  elkID <- i
  
  # subset elk.data
  indiv.data <- elk.data.all %>% filter(Animal == elkID)
  
  # define period prior to parturition
  part.date <- as.POSIXct(sites.all$t[sites.all$Animal == elkID],
                          tryFormats = "%m/%d/%Y",
                          tz = "America/New_York")
  
  indiv.data.period <- indiv.data %>% filter(t < part.date)
  
  # create a SpatialPoints object
  indiv.period.sp <- SpatialPoints(coords = indiv.data.period[ ,c("x", "y")],
                                   proj4string = UTM.projection)
  
  # fit an MCP
  indiv.mcp <- mcp(indiv.period.sp, 
                   percent = 100,
                   unin = "m",
                   unout = "km2")
  
  indiv.avail.sp <- spsample(indiv.mcp, n = 5000, type = "regular")
  
  # create data.frame
  indiv.avail <- data.frame(x = indiv.avail.sp@coords[,1],
                            y = indiv.avail.sp@coords[,2],
                            t = NA,
                            Animal = elkID,
                            Year = sites.all$Year[sites.all$Animal == elkID],
                            Case = 0)
  
  # bind to master data.frame
  available.points <- rbind(available.points, indiv.avail)
  
}

#_____________________________________________________________________________________________________________
# 5. Sample from continuous rasters ----
#_____________________________________________________________________________________________________________
# 5a. Used locations ----
#_____________________________________________________________________________________________________________

# used points 
sites.sp <- SpatialPoints(coords = sites.all[ ,c("x", "y")],
                          proj4string = UTM.projection)

# reproject to raster CRS
sites.sp.proj <- spTransform(sites.sp, raster.proj)
  
# sample from each raster
# distance rasters
sites.all$dDeveloped        <- raster::extract(dDeveloped, sites.sp.proj, method = "simple")
      
sites.all$dRoad             <- raster::extract(dRoad, sites.sp.proj, method = "simple")
      
sites.all$dOpen             <- raster::extract(dOpen, sites.sp.proj, method = "simple")
  
sites.all$dYoungForest      <- raster::extract(dYoungForest, sites.sp.proj, method = "simple")
  
sites.all$dMatureForest     <- raster::extract(dMatureForest, sites.sp.proj, method = "simple")
  
sites.all$dEdge             <- raster::extract(dEdge, sites.sp.proj, method = "simple")
  
# topography
sites.all$TPI               <- raster::extract(TPI, sites.sp.proj, method = "simple")
  
sites.all$TRI               <- raster::extract(TRI, sites.sp.proj, method = "simple")
  
sites.all$SRI               <- raster::extract(SRI, sites.sp.proj, method = "simple")
  
# continuous vegetation
sites.all$canopy            <- raster::extract(canopy, sites.sp.proj, method = "simple")

#_____________________________________________________________________________________________________________
# 5b. Available locations ----
#_____________________________________________________________________________________________________________

# used points 
avail.sp <- SpatialPoints(coords = available.points[ ,c("x", "y")],
                          proj4string = UTM.projection)

# reproject to raster CRS
avail.sp.proj <- spTransform(avail.sp, raster.proj)
  
# sample from each raster
# distance rasters
available.points$dDeveloped        <- raster::extract(dDeveloped, avail.sp.proj, method = "simple")
      
available.points$dRoad             <- raster::extract(dRoad, avail.sp.proj, method = "simple")
      
available.points$dOpen             <- raster::extract(dOpen, avail.sp.proj, method = "simple")
  
available.points$dYoungForest      <- raster::extract(dYoungForest, avail.sp.proj, method = "simple")
  
available.points$dMatureForest     <- raster::extract(dMatureForest, avail.sp.proj, method = "simple")
  
available.points$dEdge             <- raster::extract(dEdge, avail.sp.proj, method = "simple")
  
# topography
available.points$TPI               <- raster::extract(TPI, avail.sp.proj, method = "simple")
  
available.points$TRI               <- raster::extract(TRI, avail.sp.proj, method = "simple")
  
available.points$SRI               <- raster::extract(SRI, avail.sp.proj, method = "simple")
  
# continuous vegetation
available.points$canopy            <- raster::extract(canopy, avail.sp.proj, method = "simple")

#_____________________________________________________________________________________________________________
# 6. Subsample data ----
#_____________________________________________________________________________________________________________

sub.points <- data.frame()

avail.seq <- seq(25, 500, 25)

for (i in 1:10) {
  
  run <- i
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    sampled.points <- available.points %>% group_by(Animal) %>%
                                           slice_sample(n = n.avail) %>%
                                           ungroup()
    
    all.points <- rbind(sites.all, sampled.points)
    
    # add sampled n and run variables
    all.points <- all.points %>% mutate("run" = run, "samp.size" = n.avail)
    
    # bind to master df
    sub.points <- rbind(sub.points, all.points)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7. Model and extract coefficients ----
#_____________________________________________________________________________________________________________
# 7a. dDeveloped ----
#_____________________________________________________________________________________________________________ 

dDeveloped.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    # dDeveloped
    model1.struc <- glmmTMB(Case ~ dDeveloped +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dDeveloped",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dDeveloped.coef <- rbind(dDeveloped.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7b. dRoad ----
#_____________________________________________________________________________________________________________ 

dRoad.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ dRoad +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dRoad",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dRoad.coef <- rbind(dRoad.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7c. dOpen ----
#_____________________________________________________________________________________________________________ 

dOpen.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ dOpen +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dOpen",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dOpen.coef <- rbind(dOpen.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7d. dYoungForest ----
#_____________________________________________________________________________________________________________ 

dYoungForest.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ dYoungForest +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dYoungForest",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dYoungForest.coef <- rbind(dYoungForest.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7e. dMatureForest ----
#_____________________________________________________________________________________________________________ 

dMatureForest.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ dMatureForest +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dMatureForest",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dMatureForest.coef <- rbind(dMatureForest.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7f. dEdge ----
#_____________________________________________________________________________________________________________ 

dEdge.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ dEdge +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "dEdge",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    dEdge.coef <- rbind(dEdge.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7g. TRI ----
#_____________________________________________________________________________________________________________ 

TRI.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ TRI +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "TRI",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    TRI.coef <- rbind(TRI.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7h. TPI ----
#_____________________________________________________________________________________________________________ 

TPI.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ TPI +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "TPI",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    TPI.coef <- rbind(TPI.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7i. SRI ----
#_____________________________________________________________________________________________________________ 

SRI.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ SRI +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "SRI",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    SRI.coef <- rbind(SRI.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 7j. Canopy ----
#_____________________________________________________________________________________________________________ 

canopy.coef <- data.frame() 

for (i in 1:10) {
  
  run.1 <- i 
  
  for (j in avail.seq) {
    
    n.avail <- j
    
    x1.data <- sub.points %>% filter(run == run.1, samp.size == n.avail)
    
    x1.data.scaled <- x1.data %>% mutate_at(c(7:16), scale)
    
    # fit model
    model1.struc <- glmmTMB(Case ~ canopy +
                               (1 | Animal),
                             family = poisson,
                             data = x1.data.scaled,
                             doFit = FALSE)
  
    model1.struc$parameters$theta[1] = log(1e3)
    model1.struc$mapArg = list(theta = factor(c(NA)))
    
    model1 <- glmmTMB:::fitTMB(model1.struc)
    
    model.coefs <- data.frame("run" = run.1,
                              "sample" = n.avail,
                              "var" = "canopy",
                              "value" = model1$fit$par[2])
    
    # bind to master frame
    canopy.coef <- rbind(canopy.coef, model.coefs)
    
  }
  
}

#_____________________________________________________________________________________________________________
# 8. Bind all dfs together ----
#_____________________________________________________________________________________________________________ 

all.coef <- rbind(dRoad.coef, dOpen.coef, dYoungForest.coef, dMatureForest.coef,
                  dEdge.coef, TRI.coef, TPI.coef, SRI.coef, canopy.coef)

all.coef <- rbind(dRoad.coef, dOpen.coef, dMatureForest.coef,
                  TRI.coef, TPI.coef, SRI.coef, canopy.coef)

#_____________________________________________________________________________________________________________
# 9. Plot ----
#_____________________________________________________________________________________________________________ 

ggplot(data = all.coef, aes(x = sample, y = value, group = sample)) +
  theme_bw() +
  facet_wrap(~var, scales = "free_y") +
  geom_vline(xintercept = 300) +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_boxplot(alpha = 0.5) +
  ylab("Selection coefficient") +
  xlab("Availability sample per parturition site")

#_____________________________________________________________________________________________________________
# 10. Save image ----
#_____________________________________________________________________________________________________________

save.image("sensitivity.RData")

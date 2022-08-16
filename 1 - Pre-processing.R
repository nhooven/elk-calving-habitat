# Title: Parturition site selection
# Subtitle: 1 - Pre-processing (random points from MCP)
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 14 Apr 2021
# Date completed: 14 Apr 2021
# Date modified: 16 Nov 2021
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)               # work with spatial objects
library(raster)           # work with rasters
library(adehabitatHR)     # fit MCPs
library(rgeos)            # buffer
library(landscapemetrics) # sample landscape metrics
library(mefa4)            # notin function

#_____________________________________________________________________________________________________________
# 2. Read in data and work with CRSs ----
#_____________________________________________________________________________________________________________

# import "sites" dataset
sites <- read.csv("sites.csv")

# import relocation data
elk.data <- read.csv("D:/Elk project/Data analysis/Raw data processing/Relocations_vectronic_1.csv")

# get dates into correct formats
elk.data$t <- as.POSIXct(elk.data$t, tz = "America/New_York")

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

landfire.1    <- raster(paste0(raster.dir, "/", "landfire_1.tif"))

landfire.2    <- raster(paste0(raster.dir, "/", "landfire_2.tif"))

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
# 4. Fit MCPs to all previous relocations and sample available points ----
#_____________________________________________________________________________________________________________

# define number of available points to sample
n.avail <- 300

# blank data.frame to hold all points
available.points <- data.frame()

# for loop to run through each animal's data
for (i in unique(sites.1$Animal)) {
  
  elkID <- i
  
  # subset elk.data
  indiv.data <- elk.data %>% filter(Animal == elkID)
  
  # define period prior to parturition
  part.date <- as.POSIXct(sites.1$t[sites.1$Animal == elkID],
                          tryFormats = "%m/%d/%Y",
                          tz = "America/New_York")
  
  indiv.data.period <- indiv.data %>% dplyr::filter(t < part.date)
  
  # create a SpatialPoints object
  indiv.period.sp <- SpatialPoints(coords = indiv.data.period[ ,c("x", "y")],
                                   proj4string = UTM.projection)
  
  # fit an MCP
  indiv.mcp <- mcp(indiv.period.sp, 
                   percent = 100,
                   unin = "m",
                   unout = "km2")
  
  # sample 300 points
  indiv.avail.sp <- spsample(indiv.mcp, n = n.avail, type = "regular")
  
  # create data.frame
  indiv.avail <- data.frame(x = indiv.avail.sp@coords[,1],
                            y = indiv.avail.sp@coords[,2],
                            t = NA,
                            Animal = elkID,
                            Year = sites.1$Year[sites.1$Animal == elkID],
                            Case = 0)
  
  # bind to master data.frame
  available.points <- rbind(available.points, indiv.avail)
  
}

#_____________________________________________________________________________________________________________
# 5. Extract covariates/sample landscape metrics for parturition sites ----
#_____________________________________________________________________________________________________________

# create a SpatialPointsDataFrame
used.spdf <- SpatialPointsDataFrame(coords = sites.1[ ,c("x", "y")],
                                    proj4string = UTM.projection,
                                    data = sites.1)

# reproject to raster CRS
used.spdf.proj <- spTransform(used.spdf, raster.proj)

#_____________________________________________________________________________________________________________
# 5a. Sample from continuous rasters ----
#_____________________________________________________________________________________________________________

# distance rasters
used.spdf.proj@data$dDeveloped        <- raster::extract(dDeveloped, used.spdf.proj, method = "simple")
      
used.spdf.proj@data$dRoad             <- raster::extract(dRoad, used.spdf.proj, method = "simple")
      
used.spdf.proj@data$dOpen             <- raster::extract(dOpen, used.spdf.proj, method = "simple")

used.spdf.proj@data$dYoungForest      <- raster::extract(dYoungForest, used.spdf.proj, method = "simple")

used.spdf.proj@data$dMatureForest     <- raster::extract(dMatureForest, used.spdf.proj, method = "simple")

used.spdf.proj@data$dEdge             <- raster::extract(dEdge, used.spdf.proj, method = "simple")

# topography
used.spdf.proj@data$TPI               <- raster::extract(TPI, used.spdf.proj, method = "simple")

used.spdf.proj@data$TRI               <- raster::extract(TRI, used.spdf.proj, method = "simple")

used.spdf.proj@data$SRI               <- raster::extract(SRI, used.spdf.proj, method = "simple")

# continuous vegetation
used.spdf.proj@data$canopy            <- raster::extract(canopy, used.spdf.proj, method = "simple")

#_____________________________________________________________________________________________________________
# 5b. Sample landscape metrics - landfire.1 ----
#_____________________________________________________________________________________________________________

# create vector of buffer sizes
buffer.sizes <- c(250, 500, 1000)

# 6-class LANDFIRE
# sample_lsm (% young forest, % mature forest, % open, IJI, np)
used.lsm <- buffer.sizes %>% 
            set_names() %>% 
            map_dfr(~sample_lsm(landfire.1, 
                                used.spdf.proj, 
                                shape = "circle",
                                what = c("lsm_c_pland",
                                         "lsm_l_iji",
                                         "lsm_l_np"), 
                                size = .), 
                    .id = "buffer")

# create 4, 5, and 8 for those strata that don't have them
# If class 4 (successional forest) does not exist, make one
nonexistent4.250 <- function(x) {
  
  if (4 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 250])
    c(250, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
}

nonexistent4.500 <- function(x) {
  
  if (4 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 500])
    c(500, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
}

nonexistent4.1000 <- function(x) {
  
  if (4 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 1000])
    c(1000, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
}

# Apply functions over all plot IDs
used.lsm4.250LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent4.250)
used.lsm4.500LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent4.500)
used.lsm4.1000LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent4.1000)

# Bind all vectors together from list
used.lsm4.250 <- plyr::ldply(used.lsm4.250LIST, rbind)
used.lsm4.500 <- plyr::ldply(used.lsm4.500LIST, rbind)
used.lsm4.1000 <- plyr::ldply(used.lsm4.1000LIST, rbind)

# bind together dataframes for all plot_id:buffer combinations that lack 4s
used.lsm4 <- rbind(used.lsm4.250, used.lsm4.500, used.lsm4.1000)

# If class 5 (mature forest) does not exist, make one
nonexistent5.250 <- function(x) {
  
  if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 250])
    c(250, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
}

nonexistent5.500 <- function(x) {
  
  if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 500])
    c(500, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
}

nonexistent5.1000 <- function(x) {
  
  if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 1000])
    c(1000, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
}

# Apply functions over all plot IDs
used.lsm5.250LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.250)
used.lsm5.500LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.500)
used.lsm5.1000LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.1000)

# Bind all vectors together from list
used.lsm5.250 <- plyr::ldply(used.lsm5.250LIST, rbind)
used.lsm5.500 <- plyr::ldply(used.lsm5.500LIST, rbind)
used.lsm5.1000 <- plyr::ldply(used.lsm5.1000LIST, rbind)

# bind together dataframes for all plot_id:buffer combinations that lack 5s
used.lsm5 <- rbind(used.lsm5.250, used.lsm5.500, used.lsm5.1000)

# If class 6 (open) does not exist, make one
nonexistent6.250 <- function(x) {
  
  if (6 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 250])
    c(250, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
}

nonexistent6.500 <- function(x) {
  
  if (6 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 500])
    c(500, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
}

nonexistent6.1000 <- function(x) {
  
  if (6 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 1000])
    c(1000, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
}

# Apply functions over all plot IDs
used.lsm6.250LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent6.250)
used.lsm6.500LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent6.500)
used.lsm6.1000LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent6.1000)

# Bind all vectors together from list
used.lsm6.250 <- plyr::ldply(used.lsm6.250LIST, rbind)
used.lsm6.500 <- plyr::ldply(used.lsm6.500LIST, rbind)
used.lsm6.1000 <- plyr::ldply(used.lsm6.1000LIST, rbind)

# bind together dataframes for all plot_id:buffer combinations that lack 6s
used.lsm6 <- rbind(used.lsm6.250, used.lsm6.500, used.lsm6.1000)

# bind all new data frames together
used.lsm.extra <- rbind(used.lsm4, used.lsm5, used.lsm6)

# change names and bind to original data.frame
colnames(used.lsm.extra) <- colnames(used.lsm)

used.lsm <- rbind(used.lsm, used.lsm.extra)

# make sure the "value" column is numeric
used.lsm$value <- as.numeric(used.lsm$value)

# plot_id is the row index from used.spdf.proj
# get into a useable format
used.lsm.df <- data.frame(index = 1:length(unique(used.lsm$plot_id)),
                          yf.250 = NA,
                          mf.250 = NA,
                          open.250 = NA,
                          IJI.250 = NA,
                          np.250 = NA,
                          yf.500 = NA,
                          mf.500 = NA,
                          open.500 = NA,
                          IJI.500 = NA,
                          np.500 = NA,
                          yf.1000 = NA,
                          mf.1000 = NA,
                          open.1000 = NA,
                          IJI.1000 = NA,
                          np.1000 = NA)

# fill data.frame
for (y in 1:length(unique(used.lsm$plot_id))) {
  
  row.index <- y
  
  # buffer - 250
  used.lsm.df[row.index, "yf.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                     used.lsm$class == 4 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "mf.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                     used.lsm$class == 5 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "open.250"] <- (used.lsm$value[used.lsm$buffer == 250 &
                                                        used.lsm$class == 6 &
                                                        used.lsm$metric == "pland" &
                                                        used.lsm$plot_id == row.index])
  
  used.lsm.df[row.index, "IJI.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                      used.lsm$metric == "iji" &
                                                      used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "np.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                      used.lsm$metric == "np" &
                                                      used.lsm$plot_id == row.index]
  
  # buffer - 500
  used.lsm.df[row.index, "yf.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                     used.lsm$class == 4 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "mf.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                     used.lsm$class == 5 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "open.500"] <- (used.lsm$value[used.lsm$buffer == 500 &
                                                        used.lsm$class == 6 &
                                                        used.lsm$metric == "pland" &
                                                        used.lsm$plot_id == row.index])
  
  used.lsm.df[row.index, "IJI.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                      used.lsm$metric == "iji" &
                                                      used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "np.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                      used.lsm$metric == "np" &
                                                      used.lsm$plot_id == row.index]
  
  # buffer - 1000
  used.lsm.df[row.index, "yf.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                     used.lsm$class == 4 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "mf.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                     used.lsm$class == 5 &
                                                     used.lsm$metric == "pland" &
                                                     used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "open.1000"] <- (used.lsm$value[used.lsm$buffer == 1000 &
                                                        used.lsm$class == 6 &
                                                        used.lsm$metric == "pland" &
                                                        used.lsm$plot_id == row.index])
  
  used.lsm.df[row.index, "IJI.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                      used.lsm$metric == "iji" &
                                                      used.lsm$plot_id == row.index]
  
  used.lsm.df[row.index, "np.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                      used.lsm$metric == "np" &
                                                      used.lsm$plot_id == row.index]
  
}

# replace NAs in IJI with 100
used.lsm.df$IJI.250[is.na(used.lsm.df$IJI.250)] <- 100
used.lsm.df$IJI.500[is.na(used.lsm.df$IJI.500)] <- 100
used.lsm.df$IJI.1000[is.na(used.lsm.df$IJI.1000)] <- 100

#_____________________________________________________________________________________________________________
# 5d. Bind continuous covariates and landscape metrics together ----
#_____________________________________________________________________________________________________________

used.lsm.df.1 <- used.lsm.df[ ,c(2:16)]

used.data.all <- cbind(used.spdf.proj@data, used.lsm.df.1)

#_____________________________________________________________________________________________________________
# 6. Extract covariates/sample landscape metrics for random points ----
#_____________________________________________________________________________________________________________

# create a SpatialPointsDataFrame
avail.spdf <- SpatialPointsDataFrame(coords = available.points[ ,c("x", "y")],
                                     proj4string = UTM.projection,
                                     data = available.points)

# reproject to raster CRS
avail.spdf.proj <- spTransform(avail.spdf, raster.proj)

#_____________________________________________________________________________________________________________
# 6a. Sample from continuous rasters ----
#_____________________________________________________________________________________________________________

# distance rasters
avail.spdf.proj@data$dDeveloped        <- raster::extract(dDeveloped, avail.spdf.proj, method = "simple")
      
avail.spdf.proj@data$dRoad             <- raster::extract(dRoad, avail.spdf.proj, method = "simple")
      
avail.spdf.proj@data$dOpen             <- raster::extract(dOpen, avail.spdf.proj, method = "simple")

avail.spdf.proj@data$dYoungForest      <- raster::extract(dYoungForest, avail.spdf.proj, method = "simple")

avail.spdf.proj@data$dMatureForest     <- raster::extract(dMatureForest, avail.spdf.proj, method = "simple")

avail.spdf.proj@data$dEdge             <- raster::extract(dEdge, avail.spdf.proj, method = "simple")

# topography
avail.spdf.proj@data$TPI               <- raster::extract(TPI, avail.spdf.proj, method = "simple")

avail.spdf.proj@data$TRI               <- raster::extract(TRI, avail.spdf.proj, method = "simple")

avail.spdf.proj@data$SRI               <- raster::extract(SRI, avail.spdf.proj, method = "simple")

# continuous vegetation
avail.spdf.proj@data$canopy            <- raster::extract(canopy, avail.spdf.proj, method = "simple")

#_____________________________________________________________________________________________________________
# Calculate landscape metrics by individual ----
#_____________________________________________________________________________________________________________

avail.data.all.master <- data.frame()

for (z in unique(avail.spdf.proj$Animal)) {
  
  # subset based upon individual
  avail.spdf.proj.indiv <- avail.spdf.proj[avail.spdf.proj$Animal == z, ]
  
  #_____________________________________________________________________________________________________________
  # 6b. Sample landscape metrics - landfire.1 ----
  #_____________________________________________________________________________________________________________
  
  # 6-class LANDFIRE
  # sample_lsm (% young forest, % mature forest, % open, IJI)
  avail.lsm <- buffer.sizes %>% 
               set_names() %>% 
               map_dfr(~sample_lsm(landfire.1, 
                                   avail.spdf.proj.indiv, 
                                   shape = "circle",
                                   what = c("lsm_c_pland",
                                            "lsm_l_iji",
                                            "lsm_l_np"), 
                                   size = .), 
                       .id = "buffer")
  
  # create 4, 5, and 8 for those strata that don't have them
  # If class 4 (successional forest) does not exist, make one
  nonexistent4.250 <- function(x) {
    
    if (4 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 250])
      c(250, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent4.500 <- function(x) {
    
    if (4 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 500])
      c(500, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent4.1000 <- function(x) {
    
    if (4 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 1000])
      c(1000, 1, "class", 4, NA, "pland", 0, paste0(x), NA)
  }
  
  # Apply functions over all plot IDs
  avail.lsm4.250LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent4.250)
  avail.lsm4.500LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent4.500)
  avail.lsm4.1000LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent4.1000)
  
  # Bind all vectors together from list
  avail.lsm4.250 <- plyr::ldply(avail.lsm4.250LIST, rbind)
  avail.lsm4.500 <- plyr::ldply(avail.lsm4.500LIST, rbind)
  avail.lsm4.1000 <- plyr::ldply(avail.lsm4.1000LIST, rbind)
  
  # bind together dataframes for all plot_id:buffer combinations that lack 4s
  avail.lsm4 <- rbind(avail.lsm4.250, avail.lsm4.500, avail.lsm4.1000)
  
  # If class 5 (mature forest) does not exist, make one
  nonexistent5.250 <- function(x) {
    
    if (5 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 250])
      c(250, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent5.500 <- function(x) {
    
    if (5 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 500])
      c(500, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent5.1000 <- function(x) {
    
    if (5 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 1000])
      c(1000, 1, "class", 5, NA, "pland", 0, paste0(x), NA)
  }
  
  # Apply functions over all plot IDs
  avail.lsm5.250LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent5.250)
  avail.lsm5.500LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent5.500)
  avail.lsm5.1000LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent5.1000)
  
  # Bind all vectors together from list
  avail.lsm5.250 <- plyr::ldply(avail.lsm5.250LIST, rbind)
  avail.lsm5.500 <- plyr::ldply(avail.lsm5.500LIST, rbind)
  avail.lsm5.1000 <- plyr::ldply(avail.lsm5.1000LIST, rbind)
  
  # bind together dataframes for all plot_id:buffer combinations that lack 5s
  avail.lsm5 <- rbind(avail.lsm5.250, avail.lsm5.500, avail.lsm5.1000)
  
  # If class 6 (open) does not exist, make one
  nonexistent6.250 <- function(x) {
    
    if (6 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 250])
      c(250, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent6.500 <- function(x) {
    
    if (6 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 500])
      c(500, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
  }
  
  nonexistent6.1000 <- function(x) {
    
    if (6 %notin% avail.lsm$class[avail.lsm$plot_id == x & avail.lsm$buffer == 1000])
      c(1000, 1, "class", 6, NA, "pland", 0, paste0(x), NA)
  }
  
  # Apply functions over all plot IDs
  avail.lsm6.250LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent6.250)
  avail.lsm6.500LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent6.500)
  avail.lsm6.1000LIST <- lapply(min(avail.lsm$plot_id):max(avail.lsm$plot_id), nonexistent6.1000)
  
  # Bind all vectors together from list
  avail.lsm6.250 <- plyr::ldply(avail.lsm6.250LIST, rbind)
  avail.lsm6.500 <- plyr::ldply(avail.lsm6.500LIST, rbind)
  avail.lsm6.1000 <- plyr::ldply(avail.lsm6.1000LIST, rbind)
  
  # bind together dataframes for all plot_id:buffer combinations that lack 6s
  avail.lsm6 <- rbind(avail.lsm6.250, avail.lsm6.500, avail.lsm6.1000)
  
  # bind all new data frames together
  avail.lsm.extra <- rbind(avail.lsm4, avail.lsm5, avail.lsm6)
  
  # change names and bind to original data.frame
  colnames(avail.lsm.extra) <- colnames(avail.lsm)
  
  avail.lsm <- rbind(avail.lsm, avail.lsm.extra)
  
  # make sure the "value" column is numeric
  avail.lsm$value <- as.numeric(avail.lsm$value)
  
  # plot_id is the row index from avail.spdf.proj
  # get into a useable format
  avail.lsm.df <- data.frame(index = 1:length(unique(avail.lsm$plot_id)),
                             yf.250 = NA,
                             mf.250 = NA,
                             open.250 = NA,
                             IJI.250 = NA,
                             np.250 = NA,
                             yf.500 = NA,
                             mf.500 = NA,
                             open.500 = NA,
                             IJI.500 = NA,
                             np.500 = NA,
                             yf.1000 = NA,
                             mf.1000 = NA,
                             open.1000 = NA,
                             IJI.1000 = NA,
                             np.1000 = NA)
  
  # fill data.frame
  for (y in 1:length(unique(avail.lsm$plot_id))) {
    
    row.index <- y
    
    # buffer - 250
    avail.lsm.df[row.index, "yf.250"] <- avail.lsm$value[avail.lsm$buffer == 250 &
                                                       avail.lsm$class == 4 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "mf.250"] <- avail.lsm$value[avail.lsm$buffer == 250 &
                                                       avail.lsm$class == 5 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "open.250"] <- (avail.lsm$value[avail.lsm$buffer == 250 &
                                                          avail.lsm$class == 6 &
                                                          avail.lsm$metric == "pland" &
                                                          avail.lsm$plot_id == row.index])
    
    avail.lsm.df[row.index, "IJI.250"] <- avail.lsm$value[avail.lsm$buffer == 250 &
                                                        avail.lsm$metric == "iji" &
                                                        avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "np.250"] <- avail.lsm$value[avail.lsm$buffer == 250 &
                                                        avail.lsm$metric == "np" &
                                                        avail.lsm$plot_id == row.index]
    
    # buffer - 500
    avail.lsm.df[row.index, "yf.500"] <- avail.lsm$value[avail.lsm$buffer == 500 &
                                                       avail.lsm$class == 4 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "mf.500"] <- avail.lsm$value[avail.lsm$buffer == 500 &
                                                       avail.lsm$class == 5 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "open.500"] <- (avail.lsm$value[avail.lsm$buffer == 500 &
                                                          avail.lsm$class == 6 &
                                                          avail.lsm$metric == "pland" &
                                                          avail.lsm$plot_id == row.index])
    
    avail.lsm.df[row.index, "IJI.500"] <- avail.lsm$value[avail.lsm$buffer == 500 &
                                                        avail.lsm$metric == "iji" &
                                                        avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "np.500"] <- avail.lsm$value[avail.lsm$buffer == 500 &
                                                        avail.lsm$metric == "np" &
                                                        avail.lsm$plot_id == row.index]
    
    # buffer - 1000
    avail.lsm.df[row.index, "yf.1000"] <- avail.lsm$value[avail.lsm$buffer == 1000 &
                                                       avail.lsm$class == 4 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "mf.1000"] <- avail.lsm$value[avail.lsm$buffer == 1000 &
                                                       avail.lsm$class == 5 &
                                                       avail.lsm$metric == "pland" &
                                                       avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "open.1000"] <- (avail.lsm$value[avail.lsm$buffer == 1000 &
                                                          avail.lsm$class == 6 &
                                                          avail.lsm$metric == "pland" &
                                                          avail.lsm$plot_id == row.index])
    
    avail.lsm.df[row.index, "IJI.1000"] <- avail.lsm$value[avail.lsm$buffer == 1000 &
                                                        avail.lsm$metric == "iji" &
                                                        avail.lsm$plot_id == row.index]
    
    avail.lsm.df[row.index, "np.1000"] <- avail.lsm$value[avail.lsm$buffer == 1000 &
                                                        avail.lsm$metric == "np" &
                                                        avail.lsm$plot_id == row.index]
  }
  
  # replace NAs in IJI with 100
  avail.lsm.df$IJI.250[is.na(avail.lsm.df$IJI.250)] <- 100
  avail.lsm.df$IJI.500[is.na(avail.lsm.df$IJI.500)] <- 100
  avail.lsm.df$IJI.1000[is.na(avail.lsm.df$IJI.1000)] <- 100

    #_____________________________________________________________________________________________________________
    # 6d. Bind continuous covariates and landscape metrics together ----
    #_____________________________________________________________________________________________________________
    
    avail.data.all <- cbind(avail.spdf.proj.indiv@data, avail.lsm.df)
    
    # bind together
    avail.data.all.master <- rbind(avail.data.all.master, avail.data.all)
  
    
  }

#_____________________________________________________________________________________________________________
# 7. Multi-scale means ----
#_____________________________________________________________________________________________________________

# raster directory 2
raster.dir.2 <- "C:/Users/User/Desktop/GIS/Raster processing/Continuous buffers"

# Read in each raster
dDeveloped.250 <- raster(paste0(raster.dir.2, "/", "dDeveloped_250.tif"))
dDeveloped.500 <- raster(paste0(raster.dir.2, "/", "dDeveloped_500.tif"))
dDeveloped.1000 <- raster(paste0(raster.dir.2, "/", "dDeveloped_1000.tif"))

dRoad.250 <- raster(paste0(raster.dir.2, "/", "dRoad_250.tif"))
dRoad.500 <- raster(paste0(raster.dir.2, "/", "dRoad_500.tif"))
dRoad.1000 <- raster(paste0(raster.dir.2, "/", "dRoad_1000.tif"))

TRI.250 <- raster(paste0(raster.dir.2, "/", "TRI_250.tif"))
TRI.500 <- raster(paste0(raster.dir.2, "/", "TRI_500.tif"))
TRI.1000 <- raster(paste0(raster.dir.2, "/", "TRI_1000.tif"))

TPI.250 <- raster(paste0(raster.dir.2, "/", "TPI_250.tif"))
TPI.500 <- raster(paste0(raster.dir.2, "/", "TPI_500.tif"))
TPI.1000 <- raster(paste0(raster.dir.2, "/", "TPI_1000.tif"))

SRI.250 <- raster(paste0(raster.dir.2, "/", "SRI_250.tif"))
SRI.500 <- raster(paste0(raster.dir.2, "/", "SRI_500.tif"))
SRI.1000 <- raster(paste0(raster.dir.2, "/", "SRI_1000.tif"))

canopy.250 <- raster(paste0(raster.dir.2, "/", "canopy_250.tif"))
canopy.500 <- raster(paste0(raster.dir.2, "/", "canopy_500.tif"))
canopy.1000 <- raster(paste0(raster.dir.2, "/", "canopy_1000.tif"))

ED.250 <- raster(paste0(raster.dir.2, "/", "ED_250_clip.tif"))
ED.500 <- raster(paste0(raster.dir.2, "/", "ED_500_clip.tif"))
ED.1000 <- raster(paste0(raster.dir.2, "/", "ED_1000_clip.tif"))

#_____________________________________________________________________________________________________________
# 7a. Used points ----
#_____________________________________________________________________________________________________________

# create a SpatialPointsDataFrame
used.spdf.2 <- SpatialPointsDataFrame(coords = used.data.all[ ,c("x", "y")],
                                      proj4string = UTM.projection,
                                      data = used.data.all)

# reproject to raster CRS
used.spdf.2.proj <- spTransform(used.spdf.2, raster.proj)

# distance rasters
# dDeveloped
used.spdf.2.proj@data$dDeveloped.250 <- raster::extract(dDeveloped.250, 
                                                        used.spdf.2.proj, 
                                                        method = "simple")
used.spdf.2.proj@data$dDeveloped.500 <- raster::extract(dDeveloped.500, 
                                                        used.spdf.2.proj, 
                                                        method = "simple")
used.spdf.2.proj@data$dDeveloped.1000 <- raster::extract(dDeveloped.1000, 
                                                        used.spdf.2.proj, 
                                                        method = "simple")

# dRoad
used.spdf.2.proj@data$dRoad.250 <- raster::extract(dRoad.250, 
                                                   used.spdf.2.proj, 
                                                   method = "simple")
used.spdf.2.proj@data$dRoad.500 <- raster::extract(dRoad.500, 
                                                   used.spdf.2.proj, 
                                                   method = "simple")
used.spdf.2.proj@data$dRoad.1000 <- raster::extract(dRoad.1000, 
                                                    used.spdf.2.proj, 
                                                    method = "simple")
# topography
# TPI
used.spdf.2.proj@data$TPI.250 <- raster::extract(TPI.250, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$TPI.500 <- raster::extract(TPI.500, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$TPI.1000 <- raster::extract(TPI.1000, 
                                                  used.spdf.2.proj, 
                                                  method = "simple")

# TRI
used.spdf.2.proj@data$TRI.250 <- raster::extract(TRI.250, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$TRI.500 <- raster::extract(TRI.500, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$TRI.1000 <- raster::extract(TRI.1000, 
                                                  used.spdf.2.proj, 
                                                  method = "simple")

# SRI
used.spdf.2.proj@data$SRI.250 <- raster::extract(SRI.250, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$SRI.500 <- raster::extract(SRI.500, 
                                                 used.spdf.2.proj, 
                                                 method = "simple")
used.spdf.2.proj@data$SRI.1000 <- raster::extract(SRI.1000, 
                                                  used.spdf.2.proj, 
                                                  method = "simple")

# continuous vegetation
# canopy
used.spdf.2.proj@data$canopy.250 <- raster::extract(canopy.250, 
                                                    used.spdf.2.proj, 
                                                    method = "simple")
used.spdf.2.proj@data$canopy.500 <- raster::extract(canopy.500, 
                                                    used.spdf.2.proj, 
                                                    method = "simple")
used.spdf.2.proj@data$canopy.1000 <- raster::extract(canopy.1000, 
                                                     used.spdf.2.proj, 
                                                     method = "simple")

# ED
used.spdf.2.proj@data$ED.250 <- raster::extract(ED.250, 
                                                    used.spdf.2.proj, 
                                                    method = "simple")
used.spdf.2.proj@data$ED.500 <- raster::extract(ED.500, 
                                                    used.spdf.2.proj, 
                                                    method = "simple")
used.spdf.2.proj@data$ED.1000 <- raster::extract(ED.1000, 
                                                     used.spdf.2.proj, 
                                                     method = "simple")

# attribute spdf data to new df
used.data.all.2 <- used.spdf.2.proj@data

# add to simple samples 
used.data.all.3 <- cbind(used.spdf.proj@data, used.data.all.2[ ,c(7:42)])

#_____________________________________________________________________________________________________________
# 7b. Available points ----
#_____________________________________________________________________________________________________________

# create a SpatialPointsDataFrame
avail.spdf.2 <- SpatialPointsDataFrame(coords = avail.data.all.master[ ,c("x", "y")],
                                       proj4string = UTM.projection,
                                       data = avail.data.all.master)

# reproject to raster CRS
avail.spdf.2.proj <- spTransform(avail.spdf.2, raster.proj)

# distance rasters
# dDeveloped
avail.spdf.2.proj@data$dDeveloped.250 <- raster::extract(dDeveloped.250, 
                                                         avail.spdf.2.proj, 
                                                         method = "simple")
avail.spdf.2.proj@data$dDeveloped.500 <- raster::extract(dDeveloped.500, 
                                                         avail.spdf.2.proj, 
                                                         method = "simple")
avail.spdf.2.proj@data$dDeveloped.1000 <- raster::extract(dDeveloped.1000, 
                                                          avail.spdf.2.proj, 
                                                          method = "simple")

# dRoad
avail.spdf.2.proj@data$dRoad.250 <- raster::extract(dRoad.250, 
                                                    avail.spdf.2.proj, 
                                                    method = "simple")
avail.spdf.2.proj@data$dRoad.500 <- raster::extract(dRoad.500, 
                                                    avail.spdf.2.proj, 
                                                    method = "simple")
avail.spdf.2.proj@data$dRoad.1000 <- raster::extract(dRoad.1000, 
                                                     avail.spdf.2.proj, 
                                                     method = "simple")
# topography
# TPI
avail.spdf.2.proj@data$TPI.250 <- raster::extract(TPI.250, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$TPI.500 <- raster::extract(TPI.500, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$TPI.1000 <- raster::extract(TPI.1000, 
                                                   avail.spdf.2.proj, 
                                                   method = "simple")

# TRI
avail.spdf.2.proj@data$TRI.250 <- raster::extract(TRI.250, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$TRI.500 <- raster::extract(TRI.500, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$TRI.1000 <- raster::extract(TRI.1000, 
                                                   avail.spdf.2.proj, 
                                                   method = "simple")

# SRI
avail.spdf.2.proj@data$SRI.250 <- raster::extract(SRI.250, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$SRI.500 <- raster::extract(SRI.500, 
                                                  avail.spdf.2.proj, 
                                                  method = "simple")
avail.spdf.2.proj@data$SRI.1000 <- raster::extract(SRI.1000, 
                                                   avail.spdf.2.proj, 
                                                   method = "simple")

# continuous vegetation
# canopy
avail.spdf.2.proj@data$canopy.250 <- raster::extract(canopy.250, 
                                                     avail.spdf.2.proj, 
                                                     method = "simple")
avail.spdf.2.proj@data$canopy.500 <- raster::extract(canopy.500, 
                                                     avail.spdf.2.proj, 
                                                     method = "simple")
avail.spdf.2.proj@data$canopy.1000 <- raster::extract(canopy.1000, 
                                                      avail.spdf.2.proj, 
                                                      method = "simple")

# ED
avail.spdf.2.proj@data$ED.250 <- raster::extract(ED.250, 
                                                     avail.spdf.2.proj, 
                                                     method = "simple")
avail.spdf.2.proj@data$ED.500 <- raster::extract(ED.500, 
                                                     avail.spdf.2.proj, 
                                                     method = "simple")
avail.spdf.2.proj@data$ED.1000 <- raster::extract(ED.1000, 
                                                      avail.spdf.2.proj, 
                                                      method = "simple")

# attribute spdf data to new df
avail.data.all.2 <- avail.spdf.2.proj@data

# bind used and available together
part.data.1 <- rbind(used.data.all.3, avail.data.all.2)
  
#_____________________________________________________________________________________________________________
# 8. EVI means and SDs ----
#_____________________________________________________________________________________________________________

# raster data
# raster directory
EVI.dir <- "G:/Elk project/EVI rasters"

evi.508.20 <- raster(paste0(EVI.dir, "/", "5_08_20.tif"))
evi.509.21 <- raster(paste0(EVI.dir, "/", "5_09_21.tif"))

evi.524.20 <- raster(paste0(EVI.dir, "/", "5_24_20.tif"))
evi.525.21 <- raster(paste0(EVI.dir, "/", "5_25_21.tif"))

evi.609.20 <- raster(paste0(EVI.dir, "/", "6_09_20.tif"))
evi.610.21 <- raster(paste0(EVI.dir, "/", "6_10_21.tif"))

evi.625.20 <- raster(paste0(EVI.dir, "/", "6_25_20.tif"))
evi.626.21 <- raster(paste0(EVI.dir, "/", "6_26_21.tif"))

evi.711.20 <- raster(paste0(EVI.dir, "/", "7_11_20.tif"))
evi.712.21 <- raster(paste0(EVI.dir, "/", "7_12_21.tif"))

evi.727.20 <- raster(paste0(EVI.dir, "/", "7_27_20.tif"))
evi.728.21 <- raster(paste0(EVI.dir, "/", "7_28_21.tif"))

# raster proj
EVI.proj <- evi.508.20@crs

#_____________________________________________________________________________________________________________
# 8a. Loop through each animal ----
#_____________________________________________________________________________________________________________

part.data.2 <- data.frame()

# convert part.data.1$t to date
part.data.1$t <- as.POSIXct(part.data.1$t, tryFormats = "%m/%d/%Y")

for (x in unique(part.data.1$Animal)) {
  
  CollarID <- x
  
  # subset individual data
  indiv.data <- part.data.1 %>% filter(Animal == CollarID)
  
  # convert to SpatialPoints
  indiv.sp <- SpatialPoints(coords = indiv.data[c("x", "y")],
                            proj4string = UTM.projection)
  
  # project sp
  indiv.sp.proj <- spTransform(indiv.sp, EVI.proj)
  
  # define which EVI raster to use
  indiv.date <- indiv.data$t[1]
  
  if (indiv.date > as.Date("2020-05-07") & indiv.date < as.Date("2020-05-24")) {
    
    use.raster <- evi.508.20
    
  } else {
    
    if (indiv.date > as.Date("2020-05-23") & indiv.date < as.Date("2020-06-09")) {
      
      use.raster <- evi.524.20
      
    } else {
      
      if (indiv.date > as.Date("2020-06-08") & indiv.date < as.Date("2020-06-25")) {
        
        use.raster <- evi.609.20
        
      } else {
        
        if (indiv.date > as.Date("2020-06-24") & indiv.date < as.Date("2020-07-11")) {
          
          use.raster <- evi.625.20
          
        } else {
          
          if (indiv.date > as.Date("2020-07-10") & indiv.date < as.Date("2020-07-28")) {
            
            use.raster <- evi.711.20
            
          } else {
            
            if (indiv.date > as.Date("2020-07-27") & indiv.date < as.Date("2020-08-15")) {
              
              use.raster <- evi.727.20
              
            } else {
              
              if (indiv.date > as.Date("2021-05-08") & indiv.date < as.Date("2021-05-25")) {
                
                use.raster <- evi.509.21
                
              } else {
                
                if (indiv.date > as.Date("2021-05-24") & indiv.date < as.Date("2021-06-10")) {
                  
                  use.raster <- evi.525.21
                  
                } else {
                  
                  if (indiv.date > as.Date("2021-06-09") & indiv.date < as.Date("2021-06-26")) {
                    
                    use.raster <- evi.610.21
                    
                  } else {
                    
                    if (indiv.date > as.Date("2021-06-25") & indiv.date < as.Date("2021-07-12")) {
                      
                      use.raster <- evi.626.21
                      
                    } else {
                      
                      if (indiv.date > as.Date("2021-07-11") & indiv.date < as.Date("2021-07-29")) {
                        
                        use.raster <- evi.712.21 
                        
                      } else {
                        
                        use.raster <- evi.728.21
                        
                      }
                      
                    }
                    
                  }
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      }
    }
  }
  
  # extract EVI value for each location
  indiv.data$EVI <- raster::extract(use.raster, indiv.sp.proj, method = "simple")
  
  # mean of 1-pixel and 2-pixel neighborhoods
  # fit MCP to all locations
  indiv.mcp <- mcp(indiv.sp.proj, percent = 100)
  
  # buffer mcp by 2000 m
  indiv.mcp.buffer <- gBuffer(indiv.mcp, width = 2000)
  
  # crop and mask use.raster
  use.raster.crop <- crop(use.raster, indiv.mcp.buffer)
  use.raster.mask <- mask(use.raster.crop, indiv.mcp.buffer)
  
  # 1-pixel focal window
  focal.1.pix <- matrix(1, nrow = 3, ncol = 3)
  
  focal.1.raster <- focal(use.raster.mask, w = focal.1.pix, fun = mean)
  focal.1.sd <- focal(use.raster.mask, w = focal.1.pix, fun = sd)
  
  indiv.data$EVI.1pix <- raster::extract(focal.1.raster, indiv.sp.proj, method = "simple")
  indiv.data$EVI.1pix.sd <- raster::extract(focal.1.sd, indiv.sp.proj, method = "simple")
  
  # 2-pixel focal window
  focal.2.pix <- matrix(1, nrow = 5, ncol = 5)
  
  focal.2.raster <- focal(use.raster.mask, w = focal.2.pix, fun = mean)
  focal.2.sd <- focal(use.raster.mask, w = focal.2.pix, fun = sd)
  
  indiv.data$EVI.2pix <- raster::extract(focal.2.raster, indiv.sp.proj, method = "simple")
  indiv.data$EVI.2pix.sd <- raster::extract(focal.2.sd, indiv.sp.proj, method = "simple")
  
  # bind to master data frame
  part.data.2 <- rbind(part.data.2, indiv.data)
  
}

#_____________________________________________________________________________________________________________
# 9. Write all to .csv ----
#_____________________________________________________________________________________________________________

write.csv(part.data.2, "part_site_data.csv")

# Title: Parturition site selection
# Subtitle: 1b - Extra landscape metrics
# Author: Nathan D. Hooven
# Email: nathan.hooven@uky.edu
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 16 Nov 2021
# Date completed: 17 Nov 2021
# Date modified: 30 Nov 2021
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

part.data.4 <- read.csv("part_site_data.csv")

# Raster directory
raster.dir <- "G:/Elk project/Elk Zone rasters (7-20-21)"

landfire.1    <- raster(paste0(raster.dir, "/", "landfire_1.tif"))

# define raster projection
raster.proj <- landfire.1@crs

# define UTM CRS
UTM.projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")

#_____________________________________________________________________________________________________________
# 3. Calculate extra metrics per individual ----
#_____________________________________________________________________________________________________________

# create vector of buffer sizes
buffer.sizes <- c(250, 500, 1000)

extra.metrics <- data.frame()

# lists of names that are yet to go
animals <- unique(part.data.4$Animal)

animals.1 <- animals[2:46]

for (i in animals) {
  
  elkID <- i
  
  # subset data frame
  indiv.data <- part.data.4 %>% filter(Animal == elkID)
  
  # create a SpatialPoints object
  indiv.sp <- SpatialPoints(coords = indiv.data[ ,c("x", "y")],
                            proj4string = UTM.projection)
  
  # project correctly
  indiv.sp.proj <- spTransform(indiv.sp, raster.proj)
  
  # sample_lsm (CLUMPY by class, LPI by class, shape_mn by class)
  used.lsm <- buffer.sizes %>% 
              set_names() %>% 
              map_dfr(~sample_lsm(landfire.1, 
                                  indiv.sp.proj, 
                                  shape = "circle",
                                  what = c("lsm_c_ai",
                                           "lsm_l_pr"), 
                                  size = .), 
                      .id = "buffer")
  
  # change plot_id to numeric
  used.lsm$plot_id <- as.numeric(used.lsm$plot_id)
  
  # If class 5 (mf) ai does not exist, make one
  nonexistent5.250 <- function(x) {
    
    if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 250 & used.lsm$metric == "ai"])
      c(250, 1, "class", 5, NA, "ai", 0, paste0(x), NA)
  }
  
  nonexistent5.500 <- function(x) {
    
    if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 500 & used.lsm$metric == "ai"])
      c(500, 1, "class", 5, NA, "ai", 0, paste0(x), NA)
  }
  
  nonexistent5.1000 <- function(x) {
    
    if (5 %notin% used.lsm$class[used.lsm$plot_id == x & used.lsm$buffer == 1000 & used.lsm$metric == "ai"])
      c(1000, 1, "class", 5, NA, "ai", 0, paste0(x), NA)
  }
  
  # Apply functions over all plot IDs
  used.lsm5.250LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.250)
  used.lsm5.500LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.500)
  used.lsm5.1000LIST <- lapply(min(used.lsm$plot_id):max(used.lsm$plot_id), nonexistent5.1000)
  
  # Bind all vectors together from list
  used.lsm5.250 <- plyr::ldply(used.lsm5.250LIST, rbind)
  used.lsm5.500 <- plyr::ldply(used.lsm5.500LIST, rbind)
  used.lsm5.1000 <- plyr::ldply(used.lsm5.1000LIST, rbind)
  
  # bind together dataframes for all plot_id:buffer combinations that lack 5
  used.lsm5 <- rbind(used.lsm5.250, used.lsm5.500, used.lsm5.1000)
  
  # change names and bind to original data.frame
  if (is_empty(used.lsm5) == FALSE){
    
    colnames(used.lsm5) <- colnames(used.lsm)
    
  }  else {
    
    NULL
    
  }
  
  used.lsm <- rbind(used.lsm, used.lsm5)
  
  # get into a useable format
  used.lsm.df <- data.frame(index = 1:length(unique(used.lsm$plot_id)),
                             ai.mf.250 = NA,
                             ai.mf.500 = NA,
                             ai.mf.1000 = NA,
                             pr.250 = NA,
                             pr.500 = NA,
                             pr.1000 = NA)
  
  # fill data.frame
  for (y in 1:length(unique(used.lsm$plot_id))) {
    
    row.index <- y
    
    # buffer - 250
    # clumpy.mf
    used.lsm.df[row.index, "ai.mf.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                              used.lsm$class == 5 &
                                                              used.lsm$metric == "ai" &
                                                              used.lsm$plot_id == row.index]
    # pr
    used.lsm.df[row.index, "pr.250"] <- used.lsm$value[used.lsm$buffer == 250 &
                                                           used.lsm$metric == "pr" &
                                                           used.lsm$plot_id == row.index]
    
    # buffer - 500
    # clumpy.mf
    used.lsm.df[row.index, "ai.mf.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                            used.lsm$class == 5 &
                                                            used.lsm$metric == "ai" &
                                                            used.lsm$plot_id == row.index]
    # pr
    used.lsm.df[row.index, "pr.500"] <- used.lsm$value[used.lsm$buffer == 500 &
                                                         used.lsm$metric == "pr" &
                                                         used.lsm$plot_id == row.index]
    
    # buffer - 1000
    # clumpy.mf
    used.lsm.df[row.index, "ai.mf.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                            used.lsm$class == 5 &
                                                            used.lsm$metric == "ai" &
                                                            used.lsm$plot_id == row.index]
    # pr
    used.lsm.df[row.index, "pr.1000"] <- used.lsm$value[used.lsm$buffer == 1000 &
                                                         used.lsm$metric == "pr" &
                                                         used.lsm$plot_id == row.index]
  }
  
  # make sure each column is numeric
  used.lsm.df <- used.lsm.df[ ,c(2:7)] %>% mutate_all(function(x) as.numeric(x)) %>%
                                            mutate(Animal.2 = elkID)
  
  # bind together
  extra.metrics <- rbind(extra.metrics, used.lsm.df)
  
}

# replace NAs with 100
extra.metrics <- extra.metrics %>% replace(is.na(.), 0)

# bind together
part.data.5 <- cbind(part.data.4, extra.metrics[ ,c(1:6)])
  

#_____________________________________________________________________________________________________________
# 5. SD of canopy cover ----
#_____________________________________________________________________________________________________________

# raster directory 2
raster.dir.2 <- "C:/Users/User/Desktop/GIS/Raster processing/Continuous buffers"

# Read in each raster
canopy.sd.250 <- raster(paste0(raster.dir.2, "/", "canopy_sd_250.tif"))
canopy.sd.500 <- raster(paste0(raster.dir.2, "/", "canopy_sd_500.tif"))
canopy.sd.1000 <- raster(paste0(raster.dir.2, "/", "canopy_sd_1000.tif"))

# sample from each raster
all.sp <- SpatialPoints(coords = part.data.5[ ,c("x", "y")],
                        proj4string = UTM.projection)

all.sp.proj <- spTransform(all.sp, raster.proj)

part.data.5$canopy.sd.250 <- raster::extract(canopy.sd.250, 
                                             all.sp.proj, 
                                             method = "simple")

part.data.5$canopy.sd.500 <- raster::extract(canopy.sd.500, 
                                             all.sp.proj, 
                                             method = "simple")

part.data.5$canopy.sd.1000 <- raster::extract(canopy.sd.1000, 
                                              all.sp.proj, 
                                              method = "simple")

#_____________________________________________________________________________________________________________
# 4. Write to csv ----
#_____________________________________________________________________________________________________________

write.csv(part.data.5, "part_site_data.csv")

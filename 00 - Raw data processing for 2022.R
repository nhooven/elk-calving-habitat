# Title: Parturition site selection
# Subtitle: 00 - Raw data processing for 2022
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Mammal Spatial Ecology and Conservation Lab, Washington State University
# Date began: 3 Aug 2022
# Date completed: 3 Aug 2022
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(sp)                # create spatial objects
library(amt)               # create track and resample
library(rgdal)             # write shapefiles
library(lubridate)         # deal with dates
library(adehabitatHR)      # create home ranges
library(adehabitatLT)      # create trajectories
library(stringr)           # string manipulation
library(ctmm)              # apply error models

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

reloc <- read.csv("2022_reloc.csv")

# load error model
load("D:/Elk project/Data analysis/Raw data processing/vectronic_uere.RData")

#_____________________________________________________________________________________________________________
# 3. 2022 - For loop to run through all elk ----
#_____________________________________________________________________________________________________________

# Timeframe
Start.Date <- as.POSIXct("2/2/2022", tz = "America/New_York", tryFormats = "%m/%d/%Y")
End.Date <- as.POSIXct("8/1/2022", tz = "America/New_York", tryFormats = "%m/%d/%Y")

# Create blank data frames to hold all data
elk.relocations.2022 <- data.frame()

for (x in unique(reloc$Collar.ID)){
  
  #_____________________________________________________________________________________________________________
  # 3a. Select data ----
  #_____________________________________________________________________________________________________________
  
  CollarID <- x
  
  Animal  <- reloc %>% filter(Collar.ID == CollarID)
  
  #_____________________________________________________________________________________________________________
  # 3b. Filter relocations by date ----
  #_____________________________________________________________________________________________________________
  
  # coerce local date to POSIX format
  Animal$Timestamp <- as.POSIXct(ymd_hms(Animal$Acq..Time..LMT., tz = "America/New_York"))
  
  # filter dates of interest
  Animal.1 <- Animal %>% filter(Timestamp >= Start.Date & Timestamp <= End.Date)
  
  #_____________________________________________________________________________________________________________
  # 3c. Resample to ~ 13 hours and filter "bad fixes" using collar error model ----
  #_____________________________________________________________________________________________________________
  
  # filter "no fix" relocations out
  Animal.1 <- Animal.1 %>% filter(Fix.Type == "3D Validated")
  
  # define projection
  projection <- CRS("+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs")
  
  # convert to UTM
  Animal.1.sp <- SpatialPoints(coords = Animal.1[ ,c("Longitude.deg.", "Latitude.deg.")],
                               proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  Animal.1.sp.utm <- spTransform(Animal.1.sp, projection)
  
  # add in easting and northing
  Animal.1 <- Animal.1 %>% mutate(Easting = Animal.1.sp.utm@coords[ ,1],
                                  Northing = Animal.1.sp.utm@coords[ ,2])
  
  # create track in amt
  Animal.track <- make_track(Animal.1, 
                             .x = Easting,
                             .y = Northing,
                             .t = Timestamp,
                             crs = projection,
                             DOP = DOP)
  
  # resample
  Animal.track.1 <- track_resample(Animal.track,
                                   rate = hours(13),
                                   tolerance = hours(2),
                                   start = 1)
  
  # create telemetry object
  Animal.telem <- as_telemetry(Animal.track.1,
                               keep = TRUE,
                               timeformat = "auto",
                               timezone = "America/New_York",
                               projection = projection,
                               keep = TRUE)
  
  Animal.telem$HDOP <- Animal.track.1$DOP
  
  # assign error using calibration-data error model
  uere(Animal.telem) <- all.animals.uere
  
  # remove outliers with speed > 0.5 m/s and distance > 10 km
  Animal.telem.out <- outlie(Animal.telem)
  
  out.rules <- Animal.telem.out$distance > 8000 | Animal.telem.out$speed > 0.08
  
  # remove outliers from track
  Animal.track.2 <- Animal.track.1[!out.rules,]
  
  #_____________________________________________________________________________________________________________
  # 3d. Create a trajectory and bind to master df ----
  #_____________________________________________________________________________________________________________
  
  Animal.ltraj <- as.ltraj(xy = Animal.track.2[,c("x_", "y_")],
                           date = Animal.track.2$t_,
                           id = 1,
                           proj4string = projection)
  
  # Name columns "x" and "y" in trajectory
  Animal.ltraj[[1]]$x <- Animal.ltraj[[1]]$x_
  Animal.ltraj[[1]]$y <- Animal.ltraj[[1]]$y_
  
  # Create a SpatialPoints object too
  Animal.sp <- SpatialPoints(coords = Animal.track.2[,c("x_", "y_")],
                             proj4string = projection)
  
  # Create data frame and bind to master relocation data frame
  Animal.relocations <- data.frame(x = Animal.track.2$x_,
                                   y = Animal.track.2$y_,
                                   t = Animal.track.2$t_,
                                   burst = Animal.track.2$burst_,
                                   DOP = Animal.track.2$DOP,
                                   Animal = CollarID,
                                   Year = "2022")
  
  # write to master data frame
  elk.relocations.2022 <- rbind(elk.relocations.2022, Animal.relocations)
  
}

#_____________________________________________________________________________________________________________
# 4. Create season variables ----
#_____________________________________________________________________________________________________________

# create season variables (in a smart way)
elk.relocations.2022 <- elk.relocations.2022 %>% mutate(Season = dplyr::case_when(substr(elk.relocations.2022$t, 6, 7) %in%  c("02", "03", "04") ~ "WS",  
                                                                                  substr(elk.relocations.2022$t, 6, 7) %in%  c("05", "06", "07") ~ "SU", 
                                                                                  substr(elk.relocations.2022$t, 6, 7) %in%  c("08", "09", "10") ~ "UA",
                                                                                  substr(elk.relocations.2022$t, 6, 7) %in%  c("11", "12", "01") ~ "AW"))

#_____________________________________________________________________________________________________________
# 5. Bind and write to final csv ----
#_____________________________________________________________________________________________________________

write.csv(elk.relocations.2022, "Relocations_2022.csv")

# Title: Parturition site selection
# Subtitle: 13 - Map of parturition sites
# Author: Nathan D. Hooven
# Email: nathan.d.hooven@gmail.com
# Affiliation: Department of Forestry and Natural Resources, University of Kentucky
# Date began: 14 Aug 2022
# Date completed: 14 Aug 2022 
# Date modified: 
# R version: 3.6.2

#_____________________________________________________________________________________________________________
# 1. Load in required packages ----
#_____________________________________________________________________________________________________________

library(tidyverse)
library(raster)        # rasters
library(sf)            # simple features
library(ggspatial)     # scale bar and north arrow
library(spData)        # data for inset map
library(cowplot)       # draw inset map

#_____________________________________________________________________________________________________________
# 2. Read in data ----
#_____________________________________________________________________________________________________________

# site data
sites <- read.csv("sites.csv")
sites.2022 <- read.csv("sites_2022.csv")

# bind all together
all.sites <- rbind(sites[ ,c("Year", "Latitude", "Longitude")], sites.2022[ ,c("Year", "Latitude", "Longitude")])

# read in shapefile as sf
elkzone <- st_read("D:/Elk project/Elk Zone rasters (7-20-21)/elkzone_proj2.shp")

# us states from spData
data("us_states", package = "spData")

ky.state <- us_states[us_states$NAME == "Kentucky", ]

plot(ky.state)

#_____________________________________________________________________________________________________________
# 3. Convert relocations to sf and transform everything to UTM ----
#_____________________________________________________________________________________________________________

all.sites.sf <- st_as_sf(x = all.sites,
                        coords = c("Longitude", "Latitude"),
                        crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# transform 
utm.proj <- "+proj=utm +zone=17 +ellps=GRS80 +units=m +no_defs"

all.sites.sf.utm <- st_transform(all.sites.sf, utm.proj)
plot(all.sites.sf.utm)

elkzone.utm <- st_transform(elkzone, utm.proj)
plot(elkzone.utm)

ky.state.utm <- st_transform(ky.state, utm.proj)
plot(ky.state)

#_____________________________________________________________________________________________________________
# 4. Plot in ggplot ----
#_____________________________________________________________________________________________________________

# KY inset map
# extract bounding box of elkzone
elkzone.bb <- st_as_sfc(st_bbox(elkzone.utm))

map.inset <- ggplot() +
             theme_bw() +
             geom_sf(data = ky.state.utm,
                     fill = "lightgray",
                     alpha = 0.15) +
             geom_sf(data = elkzone.utm,
                     fill = "white") +
             geom_sf(data = elkzone.bb,
                     fill = NA,
                     size = 1) +
             theme(panel.grid = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),
                   plot.margin = margin(c(0, 0, 0, 0)))
             
# main map
map.main <- ggplot() +
            theme_bw() +
            geom_sf(data = elkzone.utm,
                    fill = "white") +
            geom_sf(data = all.sites.sf.utm,
                    aes(color = as.factor(Year),
                        shape = as.factor(Year)),
                    size = 2,
                    alpha = 0.5) +
            scale_colour_viridis_d(option = "plasma") +
            coord_sf(crs = "+proj=merc +lat_ts=0 +lon_0=0 +k=1.000000 +x_0=0 +y_0=0 +ellps=GRS80 +datum=WGS84 +units=m +no_defs no_defs") +
            annotation_scale(location = "br", 
                             bar_cols = c("gray", "white"),
                             style = "ticks") +
            annotation_north_arrow(style = north_arrow_fancy_orienteering,
                                   location = "br",
                                   pad_y = unit(0.25, "in")) +
            theme(legend.position = c(0.78, 0.2),
                  legend.title = element_blank(),
                  legend.background = element_rect(colour = 'black', 
                                                   fill = 'white', 
                                                   linetype='solid'))

# draw map together
ggdraw() +
  draw_plot(map.main) +
  draw_plot(map.inset,
            scale = 0.35,
            x = -0.22,
            y = 0.25)

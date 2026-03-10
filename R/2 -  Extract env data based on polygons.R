##### II -  Extract env data based on polygons ####

library(raster)
library(dplyr)
library(reshape2)
library(lubridate)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(exactextractr)
library(tidyverse)
library(ncdf4)
library(tidync)
library(terra)

# open fish data
load("Data/fish_dfclean_complete_230cells.rData")
fish_data2$year <- as.numeric(fish_data2$year)

# open and process temperature and salinity

nc.data <- terra::rast("DataBases/GLOBAL_MULTIYEAR_PHY_001_030_cmems_mod_glo_phy_my_0.083deg_P1M-m.nc")
nc.temp <- grep("bottomT", names(nc.data), value = TRUE)
print(nc.temp)
nc.temp <- subset(nc.data, nc.temp)
dim(nc.temp)

nc.sal <- grep("so", names(nc.data), value = TRUE)
print(nc.sal)
nc.sal <- subset(nc.data, nc.sal)

coords <- fish_data2 %>%
  dplyr::select(cell, longitude, latitude) %>%
  unique()
sppoints <- terra::vect(
  SpatialPoints(coords[,-1],
                proj4string = CRS(crs(nc.sal, proj = TRUE))))
temperature <- terra::extract(nc.temp,
                              sppoints,
                              df = TRUE)

salinity <- terra::extract(nc.sal,
                           sppoints,
                           df = TRUE)

# Combine and process
temp <- cbind(coords, temperature) %>%
  dplyr::select(-ID) %>%
  melt(id.vars = names(coords),
       value.name =  "Temperature",
       variable.name = "rastID") %>%
  # Split rastID into depth and time
  mutate(depth = as.numeric(sub(".*=(.*?)_.*", "\\1", rastID)),
         time_reference =  as.numeric(sub(".*_(\\d+)$", "\\1", rastID))) %>%
  dplyr::select(-rastID)

time_reference <- data.frame(date = as.Date(unique(time(nc.temp))),
                             time_reference = 1:length(unique(time(nc.temp))))

temp <- merge(temp,
             time_reference,
             by = "time_reference") %>%
  mutate(year = year(date),
         month = month(date)) %>%
  dplyr::select(-time_reference)

names(temp)[4] <- "value"

sal <- cbind(coords, salinity) %>%
  dplyr::select(-ID) %>%
  melt(id.vars = names(coords),
       value.name =  "value",
       variable.name = "rastID") %>%
  # Split rastID into depth and time
  mutate(depth = as.numeric(sub(".*=(.*?)_.*", "\\1", rastID)),
         time_reference =  as.numeric(sub(".*_(\\d+)$", "\\1", rastID))) %>%
  dplyr::select(-rastID)

time_reference <- data.frame(date = as.Date(unique(time(nc.sal))),
                             time_reference = 1:length(unique(time(nc.sal))))

sal <- merge(sal,
              time_reference,
              by = "time_reference") %>%
  mutate(year = year(date),
         month = month(date)) %>%
  dplyr::select(-time_reference)

# open and process oxygen

nc_data2 <- terra::rast("DataBases/GLOBAL_MULTIYEAR_BGC_001_029_cmems_mod_glo_bgc_my_0.25deg_P1M-m.nc")
dim(nc_data2)
names(nc_data2)[c(1, 2, 13, 14)]

nc.o2 <- subset(nc_data2, nc.o2)
dim(nc.o2)
names(nc.o2)[c(1, 2, 13, 14)]

coords <- fish_data2 %>%
  dplyr::select(cell, longitude, latitude) %>%
  unique()
names(coords)[2:3] <- c("lon", "lat")

# extract cellid from raster, NOTE that reduce the number of unique coordinates - that's because of extremely close lat and lon 
coords_with_cells <- data.frame(rastID = terra::cellFromXY(nc.o2, matrix(c(coords$lon, coords$lat), ncol = 2)))
coords_with_cells <- cbind(coords, cell_id = coords_with_cells)
cells_rastID <- unique(terra::cellFromXY(nc.o2, matrix(c(coords$lon, coords$lat), ncol = 2)))

# extract o2 values only in cells present in fish data
oxygenextracted <- extract(nc.o2, cells_rastID)

# reshape df and add lat and lon
oxygenextracted$rastID <- cells_rastID
uniq_coords_with_cells <- coords_with_cells %>%
  distinct(rastID, .keep_all = TRUE)
oxygenextracted <- left_join(oxygenextracted, uniq_coords_with_cells,
          by = "rastID")

# combine o2 data with lat, lon, and time
o2 <- oxygenextracted %>%
  melt(id.vars = c("latitude","longitude"),
       value.name =  "o2",
       variable.name = "rastID") %>%
  # Split rastID into depth and time
  mutate(depth = as.numeric(sub(".*=(.*?)_.*", "\\1", rastID)),
         time_reference =  as.numeric(sub(".*_(\\d+)$", "\\1", rastID))) %>%
  dplyr::select(-rastID)

# filter the deepest layer per coord
o2_filtered <- o2 %>%
  dplyr::filter(!is.na(o2) & !is.na(depth)) %>%
  dplyr::group_by(latitude, longitude) %>%
  dplyr::filter(depth == max(depth, na.rm = TRUE)) %>%  
  dplyr::ungroup()
  
# create a time reference df to match year and month instead of ID
time_reference <- data.frame(date = as.Date(unique(time(nc.o2))),
                             time_reference = 1:length(unique(time(nc.o2))))

# merge back to env
o2_filtered <- merge(o2_filtered,
             time_reference,
             by = "time_reference") %>%
  mutate(year = year(date),
         month = month(date)) %>%
  dplyr::select(-time_reference)

# get hexagons id (from fish data)
o2 <- left_join(o2_filtered, coords, by = c("latitude", "longitude"))
names(o2)[3] <- "value"
names(o2)[1:2] <- c("latitude","longitude")

# save environmental variables long format
save(temp, file = "Data/templong_matchcleandf_230cells.rData")
save(sal, file = "Data/sallong_matchcleandf_230cells.rData")
save(o2, file = "Data/o2long_matchcleandf_230cells.rData")

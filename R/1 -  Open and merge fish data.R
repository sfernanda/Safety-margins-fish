##### I -  Open and process fish data ####
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(sf)
library(dggridR)

load(file = "Data/fish_data.rData")

# check
coord <- fish_data %>%
  dplyr::select(longitude, latitude) %>%
  unique()

# Create grid hexagons for fish data 
# Construct a global grid with cells 5x5 degree resolution
dggs <- dgconstruct(aperture = 3, metric = TRUE, 
                    resround = 'down', topology = "HEXAGON", res = 6)

# Get the corresponding grid cells (lat-long pair)
coord$cell <- dgGEO_to_SEQNUM(dggs, coord$longitude, coord$latitude)$seqnum

# Converting SEQNUM to GEO gives the center coordinates of the cells
cellcenters <- dgSEQNUM_to_GEO(dggs, coord$cell)

# Add centroid coordinates
coord$Lon_c <- cellcenters$lon_deg
coord$Lat_c <- cellcenters$lat_deg

# Get the grid cell boundaries for cells
dGrid <- dgcellstogrid(dggs, coord$cell)

# Merge hexagon grid back to the unique coordinates
head(coord)
head(dGrid)
length(unique(coord$cell))
length(unique(dGrid$seqnum))
dHex <- left_join(coord, dGrid, by = c("cell" = "seqnum"))

nrow(dHex) == nrow(coord)

fish_data2 <- fish_data %>%
  left_join(dHex, by = c("longitude", "latitude"))

save(fish_data2, file = "Data/fish_dfclean_complete_230cells.rData")


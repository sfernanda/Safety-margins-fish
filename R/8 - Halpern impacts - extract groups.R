##### VIII - Extracting Halpern impacts manually ##### 

library(raster)
library(tibble)
library(sf)
library(ggplot2)
library(patchwork)
library(dplyr)
library(terra)
library(exactextractr)
library(tidyr)
library(purrr)

##### Relevant info ####
# The impacts are divided in 4 categories: 
# •	Fishing: commercial demersal destructive, commercial demersal nondestructive high bycatch, commercial demersal nondestructive low bycatch, pelagic high bycatch, pelagic low bycatch, artisanal
# •	Climate change: sea surface temperature, ocean acidification, sea level rise
# •	Ocean: shipping
# •	Land-based: nutrient pollution, organic chemical pollution, direct human, light
#
# Averages and slopes are extracted for models:
# Only fishing is extracted individually - average within each grid cell per year (then tot average across years and slope)
# All other groups are cumulative (sum) per yr - then avg across all years within each cell and slope of the trend 


##### 1. List folders to open ####
folders <- list.dirs("Data/Halpern", recursive = FALSE) # download all the data from Halpern et al.2019 repository

grouped_folders <- list(
  fishing = folders[grepl("fishing", folders)],
  climate = folders[grepl("SST|OA|SLR", folders)],
  ocean = folders[grepl("Shipping", folders)],
  landbased = folders[grepl("Nutrient_pollution|Organic_chemical_pollution|Direct_human|Light", folders)])

# Append "/data" to each folder path
grouped_folders <- lapply(grouped_folders, function(x) file.path(x, "data"))

##### 2. Open rasters and sum the impacts ####
# Years to process
years <- 2003:2013

# Example group names - replace with your actual group names
group_names <- names(grouped_folders)

# Initialize list to store summed rasters per group
grouped_yearly_sums <- list()

# Process each group one by one, manually
for (i in seq_along(group_names)) {
  group_name <- group_names[i]
  cat("\n=== Processing group:", group_name, "===\n")
  
  data_folders <- grouped_folders[[group_name]]
  
  # Initialize list for yearly sums in this group
  yearly_sums <- list()
  
  # Process each year manually
  for (j in seq_along(years)) {
    year <- years[j]
    cat("\nProcessing year:", year, "for group:", group_name, "\n")
    
    # Find raster files in all folders of this group for this year
    rasters <- c()  # start empty
    for (folder in data_folders) {
      if (dir.exists(folder)) {
        files <- list.files(path = folder, pattern = paste0(year, "_impact\\.tif$"), full.names = TRUE)
        rasters <- c(rasters, files)
      }
    }
   
    # Read rasters
    import_rast <- list()
    for (k in seq_along(rasters)) {
      cat("  Reading raster:", basename(rasters[k]), "\n")
      import_rast[[k]] <- rast(rasters[k])
    }
    
    # Stack rasters
    rast_stack <- do.call(c, import_rast)
    
    # Sum rasters
    summed_rast <- app(rast_stack, sum, na.rm = TRUE)
    
    # Store summed raster
    yearly_sums[[as.character(year)]] <- summed_rast
  }
  
  # Save results for this group
  grouped_yearly_sums[[group_name]] <- yearly_sums
}

# save individually, as the binded list is too heavy
fishing <- grouped_yearly_sums$fishing
climate <- grouped_yearly_sums$climate
ocean   <- grouped_yearly_sums$ocean
landbased <-grouped_yearly_sums$landbased

##### 3. Get data only for polygons with fish ####
load("Data/fish_dfclean_complete_230cells.rData") # hexagon geometries

# Get unique hexagons
coord_poly <- fish_data2 %>%
  dplyr::select(cell, geometry) %>%
  distinct(cell, .keep_all = TRUE)
coord_poly <- st_as_sf(coord_poly)

first_raster <- fishing[[1]]
coord_poly <- st_transform(coord_poly, crs = st_crs(crs(first_raster))) # reproject polygins to match the raster, not the oposite

# Function to extract mean values per polygon for each year in each category of the list
extract_means_by_category <- function(category_list, polygons_sf) {
  result_list <- list()
  
  for (year_name in names(category_list)) {
    cat("Processing year:", year_name, "\n")
    
    raster_layer <- raster::raster(category_list[[year_name]])
    
    means_vals <- exact_extract(raster_layer, polygons_sf, 'mean')
    
    means_sf <- polygons_sf %>%
      st_drop_geometry() %>%
      mutate(mean = means_vals) %>%
      bind_cols(geometry = st_geometry(polygons_sf)) %>%
      st_as_sf()
    
    result_list[[year_name]] <- means_sf
  }
  
  return(result_list)
}

# Now apply to all categories
fishing_hex <- extract_means_by_category(fishing, coord_poly)
climate_hex <- extract_means_by_category(climate, coord_poly)
ocean_hex   <- extract_means_by_category(ocean, coord_poly)
landbased_hex <- extract_means_by_category(landbased, coord_poly)

# Combine all years into data.frame for all the groups
# for all impacts
fishing_hex_df <- do.call(cbind, fishing_hex) %>% as.data.frame()
colnames(fishing_hex_df) <- paste0("mean_", names(fishing))
fishing_hex_df <- bind_rows(lapply(names(fishing_hex), function(y) {
    fishing_hex[[y]] %>% mutate(year = as.integer(y))
  }))

climate_hex_df <- do.call(cbind, climate_hex) %>% as.data.frame()
colnames(climate_hex_df) <- paste0("mean_", names(climate))
climate_hex_df <- bind_rows(lapply(names(climate_hex), function(y) {
    climate_hex[[y]] %>% mutate(year = as.integer(y))
  }))

ocean_hex_df <- do.call(cbind, ocean_hex) %>% as.data.frame()
colnames(ocean_hex_df) <- paste0("mean_", names(ocean))
ocean_hex_df <- bind_rows(lapply(names(ocean_hex), function(y) {
    ocean_hex[[y]] %>% mutate(year = as.integer(y))
  }))

landbased_hex_df <- do.call(cbind, landbased_hex) %>% as.data.frame()
colnames(landbased_hex_df) <- paste0("mean_", names(landbased))
landbased_hex_df <- bind_rows(lapply(names(landbased_hex), function(y) {
    landbased_hex[[y]] %>% mutate(year = as.integer(y))
  }))

# ggplot(fishing_hex_df %>% filter(year == 2013))+
#   geom_sf(aes(geometry = geometry, fill = mean))


##### 4. Save data extracted per year per grid cell ####
save(fishing_hex_df, file = "Data/fishing_hex_df.rData")
save(climate_hex_df, file = "Data/climate_hex_df.rData")
save(ocean_hex_df, file = "Data/ocean_hex_df.rData")
save(landbased_hex_df, file = "Data/landbased_hex_df.rData")

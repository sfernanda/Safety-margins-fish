##### III - Calculations of max/min environment data #####
library(dplyr)
library(sf)

# raw temperature and salinity data (same resolution)
load(file = "Data/templong_matchcleandf_230cells.rData")
load(file = "Data/sallong_matchcleandf_230cells.rData")
load(file = "Data/o2long_matchcleandf_230cells.rData")

# Function to calculate the warmest and coldest month per year and decade
calculate_extreme_means_per_decade <- function(data) {
  # Calculate warmest and coldest month per year
  annual_extremes <- data %>%
    group_by(year, cell) %>%
    summarise(
    max_value = if (any(!is.na(value))) max(value, na.rm = TRUE) else NA_real_,  # Warmest month
    min_value = if (any(!is.na(value))) min(value, na.rm = TRUE) else NA_real_,  # Coldest month
    longitude = first(longitude),
    latitude = first(latitude),
    .groups = "drop")
  
  # add decade column
  annual_extremes <- annual_extremes %>%
    mutate(decade = paste0((year - min(year)) %/% 10 * 10 + min(year), "_", 
                           (year - min(year)) %/% 10 * 10 + min(year) + 9))
  
  # Calculate the decade averages
  decade_extremes <- annual_extremes %>%
    group_by(decade, cell) %>%
    summarise(
      avg_max_value = mean(max_value, na.rm = TRUE), # Average of warmest months per decade
      avg_min_value = mean(min_value, na.rm = TRUE), # Average of coldest months per decade
      longitude = first(longitude),
      latitude = first(latitude),
      decade = first(decade),
      .groups = "drop")
  
  return(decade_extremes)
}

# Apply the function to env data
extremestemp_dec <- calculate_extreme_means_per_decade(temp)
extremessal_dec <- calculate_extreme_means_per_decade(sal)
extremeso2_dec <- calculate_extreme_means_per_decade(o2)

# merge env data and ecoregion data to fish data
load(file = "Data/fish_dfclean_complete_230cells.rData")

fish_data2$hex_geometry <- fish_data2$geometry
fish_data2_sf <- st_as_sf(fish_data2, coords = c("longitude", "latitude"), crs = 4326)
fish_data2_sf$year <- as.numeric(fish_data2_sf$year)
fish_data2_sf <- fish_data2_sf %>%
  mutate(decade = paste0((year - min(year)) %/% 10 * 10 + min(year), "_", 
                         (year - min(year)) %/% 10 * 10 + min(year) + 9))

fish_env <- left_join(fish_data2_sf, extremestemp_dec, by = c("decade", "cell"))
names(fish_env)[grep("avg_max_value", colnames(fish_env))] <- "max_temp_dec"
names(fish_env)[grep("avg_min_value", colnames(fish_env))] <- "min_temp_dec"

fish_env <- left_join(fish_env, extremessal_dec, by = c("decade", "cell"))
names(fish_env)[grep("avg_max_value", colnames(fish_env))] <- "max_sal_dec"
names(fish_env)[grep("avg_min_value", colnames(fish_env))] <- "min_sal_dec"

fish_env <- left_join(fish_env, extremeso2_dec, by = c("decade", "cell"))
names(fish_env)[grep("avg_max_value", colnames(fish_env))] <- "max_o2_dec"
names(fish_env)[grep("avg_min_value", colnames(fish_env))] <- "min_o2_dec"

fish_env2 <- fish_env[,-c(19,20,23,24)]
save(fish_env2, file = "Data/fish_env_cleandf_230cells.rData")

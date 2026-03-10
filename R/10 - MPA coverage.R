##### X - Test relationship of MPA on ESM  ##### 

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)

# Disable s2 to prevent geometry errors with complex polygons
sf::sf_use_s2(FALSE)

# Load fish hexagon polygons
load("Data/fish_dfclean_complete_230cells.rData")

# Extract unique hexagon geometries
fish_poly <- fish_data2 %>%
  dplyr::select(cell, geometry) %>%
  distinct(cell, .keep_all = TRUE) %>%
  st_as_sf()
fish_poly <- st_wrap_dateline(fish_poly, options = c("WRAPDATELINE=YES",
                                                     "DATELINEOFFSET=180"), quiet = TRUE)
fish_poly_equal <- st_cast(fish_poly, "POLYGON", explode = TRUE)
fish_poly_equal$area_km2 <- as.numeric(st_area(fish_poly_equal)) / 10^6
fish_poly_equal <- fish_poly_equal %>%
  arrange(area_km2) %>%
  slice(-(1:2)) # remove cells in the edge
plot(fish_poly_equal)
summary(fish_poly_equal)

# Get bounding box of your study area
bbox_fish <- st_bbox(fish_poly_equal)

# Read the full MPA shapefile
mpa_all0 <- st_read("Data/MPAs/marine_shp_0/WDPA_WDOECM_May2025_Public_marine_shp-polygons.shp")
mpa_all1 <- st_read("Data/MPAs/marine_shp_1/WDPA_WDOECM_May2025_Public_marine_shp-polygons.shp")
mpa_all2 <- st_read("Data/MPAs/marine_shp_2/WDPA_WDOECM_May2025_Public_marine_shp-polygons.shp")

# Merge all into a single sf object
mpa_all <- rbind(mpa_all0, mpa_all1, mpa_all2)

# Filter: keep only designated marine MPAs, exclude proposed/inscribed
mpa2 <- mpa_all %>%
  filter(MARINE == 2) #!STATUS %in% c("Proposed", "Inscribed"))

# Crop MPA to area of fish polygons
mpa_crop <- st_crop(mpa2, bbox_fish)

# Ensure valid geometries
mpa_crop <- st_make_valid(mpa_crop)
plot(mpa_crop$geometry)
mpa_merged <- st_union(mpa_crop)
plot(mpa_merged)

# Ensure both layers have the same CRS
fish_poly_equal <- st_transform(fish_poly_equal, crs = st_crs(mpa_merged))

# Intersect fish polygons with MPAs
intersections <- st_intersection(fish_poly_equal, mpa_merged)
length(unique(intersections$cell))

# Calculate area of intersection
intersections2 <- intersections
intersections2$area_intersect <- st_area(intersections2)
intersections2$area_intersect <- as.numeric(intersections2$area_intersect) / 1e6

# Summarize intersection area per cell
intersections2$percent_cover <- intersections2$area_intersect/intersections2$area_km2*100
intersections2 <- intersections2 %>% 
  st_drop_geometry()

# Join back to fish polygons
fish_poly2 <- left_join(fish_poly, intersections2, by = "cell")

# Replace NA with 0 for cells with no intersection
fish_poly2$percent_cover[is.na(fish_poly2$percent_cover)] <- 0

# Final result
world_map <- map_data("world")

# save(fish_poly2, file = "Data/MPA_coverage.rData")
load("Data/MPA_coverage.rData")

# Open safety margins metrics after permutations to calculate cumulative safety
load(file = "Data/resultsbiomass_list_cwm_230cells.rData")

# Combine all results into one data frame
final_results <- bind_rows(results_list)
final_results <- final_results %>%
  separate(cell_decade, into = c("cell", "decade"), sep = "_", convert = TRUE)

# mean of all metrics after randomization
final_results2 <- final_results %>%
  dplyr::group_by(cell, hex_geometry) %>%
  dplyr::summarise(
    # preferred tolerances
    usm_temp_mean_pref = mean(usm_temp_pref, na.rm = TRUE),
    lsm_temp_mean_pref = mean(lsm_temp_pref, na.rm = TRUE),
    lsm_sal_mean_pref = mean(lsm_sal_pref, na.rm = TRUE),
    lsm_o2_mean_pref = mean(lsm_o2_pref, na.rm = TRUE),
    .groups = "drop")
str(final_results2)

# rescale safety margin results in a scale 0-1
final_results3 <- final_results2 %>%
  mutate(across(c(usm_temp_mean_pref, lsm_temp_mean_pref,
                  lsm_sal_mean_pref, lsm_o2_mean_pref), ~ scales::rescale(.)))
names(final_results3)[2] <- "geometry"

final_results3 <- final_results3 %>%
  group_by(cell, geometry) %>%
  dplyr::summarise(Safety = usm_temp_mean_pref + lsm_temp_mean_pref + lsm_sal_mean_pref + lsm_o2_mean_pref)

mpa_safety <- merge(final_results3, fish_poly2, by = "cell")
hist(mpa_safety$Safety)
mpa_safety$percent_cover[is.na(mpa_safety$percent_cover)] <- 0

library(mgcv)
mpagam <- gam(data = mpa_safety, Safety ~ s(percent_cover, k = 3))
summary(mpagam)
hist(resid(mpagam))
plot(mpagam)

plotmpa <- ggplot(mpa_safety, aes(y = Safety, x = percent_cover)) +
  geom_point(alpha = .3, size = 2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3), col = "#B84193") +
  labs(y = "Cumulative Safety Margin", x = "Percentage of area covered by MPA") +
  theme_minimal()

##### Plot of fig.2d ####
load(file = "Data/proportion_above0.rData")
load(file = "Data/MPA_coverage.rData")
fish_poly2 <- st_drop_geometry(fish_poly2)
fish_poly2$percent_cover[is.na(fish_poly2$percent_cover)] <- 0

average_proportions <- avg_props_long %>%
  pivot_wider(
    names_from = metric,
    values_from = value)
propsafe_mpa <- merge(fish_poly2, average_proportions, by = "cell")

safety_prop <- merge(propsafe_mpa, final_results3, by = "cell")
safety_prop$Safety <- (safety_prop$Safety - min(safety_prop$Safety, na.rm = TRUE)) /
  (max(safety_prop$Safety, na.rm = TRUE) - min(safety_prop$Safety, na.rm = TRUE)) # standardize 0-1

safety_prop$quadrant <- with(safety_prop, ifelse(prop_safe_avg4 >= 0.5 & Safety >= 0.5, 2,
                                                 ifelse(prop_safe_avg4 < 0.5 & Safety >= 0.5, 1,
                                                        ifelse(prop_safe_avg4 < 0.5 & Safety < 0.5, 3,
                                                               4))))
safety_prop$quadrant <- factor(safety_prop$quadrant,
                               levels = c(1, 2, 3, 4),
                               labels = c("I.", "II.", "III.", "IV."))


# Scatter plot with quadrant colors + MPA % coverage
safety_prop2 <- safety_prop %>%
  dplyr::filter(metric == "all_positive")

ggplot(safety_prop, aes(x = prop_safe_avg4, y = Safety)) +
  annotate("text", x = 0.01, y = 0.99, label = "I", size = 5, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = 0.51, y = 0.99, label = "II", size = 5, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = 0.01, y = 0.49, label = "III", size = 5, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text", x = 0.51, y = 0.49, label = "IV", size = 5, fontface = "bold", hjust = 0, vjust = 1) +
  
  geom_point(aes(color = quadrant, size = percent_cover), alpha = 0.5, show.legend = c(color = FALSE, size = TRUE)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  xlim(0, 1) + ylim(0, 1) +
  scale_color_manual(values = c(
    "III." = "#432818",
    "I." = "#D7AA74",
    "IV." = "#E5DE67",
    "II." = "#bb9457")) +
  labs(x = "Proportion of species above 0",
       y = "Standardized ESMs",
       size = "% area \ncovered \nby MPA") +
  theme_minimal()

##### V - Calculate metrics permutating data using SACs #####
library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr)
library(reshape2)
library(vegan)
library(legendry)
library(scico)
library(patchwork)  # for `inset_element()`

load(file = "Data/filtered_data2_sac_hex_230cells.rData") # filtered after SAC
load(file = "Data/cell_names_hex_230cells.rData")
filtered_data2 <- filtered_data2[, c("haul_id","Species","decade","cell","hex_geometry",
                                     "max_temp_dec","min_temp_dec","max_sal_dec","min_sal_dec","max_o2_dec","min_o2_dec",
                                     "TempPrefMax","TempPrefMin","SalinityPrefMax","SalinityPrefMin","OxyPrefMax","OxyPrefMin")]

ind_metrics_list <- list()

for (j in 1:99){
  print(j)
  
  sub <- subset(filtered_data2, filtered_data2$haul_id %in% cell_names[,j])
  sub <- sub[!(is.na(sub$TempPrefMax) & is.na(sub$TempPrefMin) &
                 is.na(sub$SalinityPrefMax) & is.na(sub$SalinityPrefMin) &
                 is.na(sub$OxyPrefMax) & is.na(sub$OxyPrefMin)), ]
  
  # Calculate metrics
  sub$usm_temp_pref <- sub$TempPrefMax - sub$max_temp_dec
  sub$lsm_temp_pref <- sub$min_temp_dec - sub$TempPrefMin
  sub$lsm_sal_pref <- sub$min_sal_dec - sub$SalinityPrefMin
  sub$lsm_o2_pref   <- sub$min_o2_dec - sub$OxyPrefMin
  
  # Keep only unique species per haul
  sub <- sub[!duplicated(sub[, c("decade", "haul_id", "Species")]), ]
  
  # Store result in the list
  ind_metrics_list[[j]] <- sub
}

# Create a list of iteration numbers corresponding to each element in results_list
iterations <- seq_along(ind_metrics_list)

# Add iteration number to each data frame in the list
ind_metrics_list_with_iterations <- map2(ind_metrics_list, 99,
                                         ~mutate(.x, iteration = .y))

# Calculate proportions for each list element - group by cells or decades + cell
proportions_list <- lapply(ind_metrics_list_with_iterations, function(data) {
  data %>%
    group_by(decade, cell, hex_geometry) %>% 
    summarise(prop_usm_temp_pref = sum(usm_temp_pref > 0, na.rm = TRUE) / n(),
              prop_lsm_temp_pref = sum(lsm_temp_pref > 0, na.rm = TRUE) / n(),
              prop_lsm_sal_pref = sum(lsm_sal_pref > 0, na.rm = TRUE) / n(),
              prop_lsm_o2_pref = sum(lsm_o2_pref > 0, na.rm = TRUE) / n(),
              prop_all_positive = sum(
                usm_temp_pref > 0 & 
                  lsm_temp_pref > 0 & 
                  lsm_sal_pref > 0 & 
                  lsm_o2_pref > 0, na.rm = TRUE) / n(),
              .groups = "drop")
})

# Combine the results from all list elements
combined_proportions <- bind_rows(proportions_list)

# Calculate the average of all proportions across the list
average_proportions <- combined_proportions %>%
  group_by(cell, hex_geometry) %>%
  summarise(
    avg_prop_usm_temp_pref = mean(prop_usm_temp_pref, na.rm = TRUE),
    avg_prop_lsm_temp_pref = mean(prop_lsm_temp_pref, na.rm = TRUE),
    avg_prop_lsm_sal_pref = mean(prop_lsm_sal_pref, na.rm = TRUE),
    avg_prop_lsm_o2_pref = mean(prop_lsm_o2_pref, na.rm = TRUE),
    avg_prop_all_positive = mean(prop_all_positive, na.rm = TRUE,
                                 .groups = "drop"))
average_proportions$prop_safe_avg4 <- rowMeans(st_drop_geometry(average_proportions[, 3:6]), na.rm = TRUE)

avg_props_long <- average_proportions %>%
  pivot_longer(cols = starts_with("avg_prop_"),
               names_to = "metric",
               values_to = "value")
avg_props_long$metric <- gsub("avg_prop_", "", avg_props_long$metric)

save(avg_props_long, file = "Data/proportion_above0.rData")
load(file = "Data/proportion_above0.rData")

# Plot map
world_map <- map_data("world")
average_proportions <- st_as_sf(avg_props_long)
average_proportions <- st_wrap_dateline(average_proportions, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

# Create a custom histogram guide
library(legendry)
library(scico)
library(patchwork)  # for `inset_element()`

# Step 1: Create 6 bins for cumulative impact, remove brackets from labels
breaks <- seq(min(average_proportions$prop_safe_avg4, na.rm = TRUE),
              max(1, na.rm = TRUE),
              length.out = 7)

# Define the custom labels for bins, removing the brackets
labels <- paste0(round(breaks[-length(breaks)] + 0.01, 1), "–", round(breaks[-1] - 0.01, 1))

# Apply the cut with custom labels
average_proportions$impact_bin <- cut(average_proportions$prop_safe_avg4,
                                      breaks = breaks,
                                      include.lowest = TRUE,
                                      right = TRUE,
                                      labels = labels)

# Ensure the factor is clean
average_proportions$impact_bin <- factor(average_proportions$impact_bin)

# Use scico palette with 6 colors
fill_colors <- c("#540b0e", "#9e2a2b", "#e09f3e", "#fff3b0", "#a3b18a", "#3a5a40")
names(fill_colors) <- levels(average_proportions$impact_bin)

# Create map without legend
map_prop <- ggplot(average_proportions) +
  geom_sf(aes(geometry = hex_geometry, fill = prop_safe_avg4 ), color = NA) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray75") +
  coord_sf(xlim = c(-180, 60), ylim = c(20, 85), expand = FALSE) +
  # scico::scale_fill_scico(palette = "buda", name = "Prop of safe", limits = c(0, 1)) +
  scale_fill_gradientn(
    colours = fill_colors,
    name = "Prop of safe",
    limits = c(0, 1)) + 
  theme_minimal() +
  labs(x = "", y = "", title = "(b) Proportion of species with metrics within the safe zone") +
  theme(legend.position = "none")

# Create histogram legend
hist_legend_prop <- ggplot(average_proportions, aes(x = impact_bin)) +
  geom_bar(aes(fill = impact_bin), color = "white", size = 0.01) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  theme_void() +
  scale_x_discrete(labels = labels) +  # Custom labels without brackets
  labs(x = "Proportion of safe species") +
  theme(axis.text.x = element_text(angle = 30, size = 4),
        axis.title.x = element_text(size = 5))

# Combine
map_prop + inset_element(hist_legend_prop, left = 0.02, bottom = 0.01, right = 0.2, top = 0.45)
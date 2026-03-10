##### VI - Metrics weighted by biomass #####

library(ggplot2)
library(sf)
library(dplyr)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr)
library(reshape2)
library(vegan)
library(FD)
library(BAT)
library(picante)
library(stringr)
library(ggdist)

load(file = "Data/filtered_data2_sac_hex_230cells.rData") # filtered after SAC
load(file = "Data/cell_names_hex_230cells.rData")

results_list <- list()

for (j in 1:3){
    print(j)  
  
  sub <- subset(filtered_data2, filtered_data2$haul_id %in% cell_names[,j])
  sub <- st_drop_geometry(sub)
  biommatrix <- sub[,c("Species","wgt_cpua","cell_decade")]
  
  traitmatrix <- sub[,c(10,25,26,29,30,33,34)]
  traitmatrix <- as.data.frame(traitmatrix[!duplicated(traitmatrix),])
  traitmatrix <- na.omit(traitmatrix)
  biommatrix <- biommatrix[biommatrix$Species %in% traitmatrix$Species,]

  cwm_values <- biommatrix %>%
    left_join(traitmatrix, by = "Species") %>%
    group_by(cell_decade)  %>%
    summarise(across(where(is.numeric), ~ sum(.x * wgt_cpua, na.rm = TRUE) / sum(wgt_cpua, na.rm = TRUE), .names = "CWM_{.col}"), 
              .groups = "drop")
  
  sub2 <- as.data.frame(sub[, c("max_temp_dec", "min_temp_dec", "min_sal_dec", "min_o2_dec", 
                                "cell_decade", "hex_geometry", "decade")])
  sub2 <- sub2[!duplicated(sub2[, -which(names(sub2) == "hex_geometry")]), ]
  sub2 <- left_join(cwm_values, sub2, by = "cell_decade")
  
  
  # Calculate metrics
  metrics <- sub2 %>%
    dplyr::group_by(cell_decade) %>%
    dplyr::mutate(# using pref trait tolerance values
                  usm_temp_pref = CWM_TempPrefMax - max_temp_dec,
                  lsm_temp_pref = min_temp_dec - CWM_TempPrefMin,
                  lsm_sal_pref = min_sal_dec - CWM_SalinityPrefMin,
                  lsm_o2_pref = min_o2_dec - CWM_OxyPrefMin)
  
  # Store result in the list
  results_list[[j]] <- metrics
}

# save(results_list, file = "Data/resultsbiomass_list_cwm_230cells.rData")
load(file = "Data/resultsbiomass_list_cwm_230cells.rData")

# Create a list of iteration numbers corresponding to each element in results_list
iterations <- seq_along(results_list)

# Add iteration number to each data frame in the list
results_list_with_iterations <- map2(results_list, 5,
                                     ~mutate(.x, iteration = .y))

# Combine all results into one data frame
final_results <- bind_rows(results_list_with_iterations)

# mean of all metrics after randomization
final_results2 <- final_results %>%
  dplyr::group_by(hex_geometry) %>%
  dplyr::summarise(
    usm_temp_mean_pref = mean(usm_temp_pref, na.rm = TRUE),
    lsm_temp_mean_pref = mean(lsm_temp_pref, na.rm = TRUE),
    lsm_sal_mean_pref = mean(lsm_sal_pref, na.rm = TRUE),
    lsm_o2_mean_pref = mean(lsm_o2_pref, na.rm = TRUE))

# rescale safety margin results in a scale 0-1
final_results2 <- final_results2 %>%
  mutate(across(c(usm_temp_mean_pref, lsm_temp_mean_pref,
                  lsm_sal_mean_pref, lsm_o2_mean_pref), ~ (.-min(., na.rm = TRUE)) / (max(., na.rm = TRUE) - min(., na.rm = TRUE))))


final_results2 <- st_as_sf(final_results2)
final_results2 <- st_wrap_dateline(final_results2, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

final_results2$cumulative_sm <- final_results2$usm_temp_mean_pref + 
  final_results2$lsm_temp_mean_pref + 
  final_results2$lsm_sal_mean_pref + 
  final_results2$lsm_o2_mean_pref
final_results2 <- final_results2 %>%
  mutate(cumulative_sm = (cumulative_sm - min(cumulative_sm, na.rm = TRUE)) /
           (max(cumulative_sm, na.rm = TRUE) - min(cumulative_sm, na.rm = TRUE)))

# Plot maps to visualize results
countries <- ne_countries(scale = "medium", returnclass = "sf")
countries <- map_data("world")

final_results_long <- final_results2 %>%
  select(hex_geometry, cumulative_sm) %>%
  pivot_longer(cols = -hex_geometry, names_to = "Variable", values_to = "Value")

# Create 6 bins for cumulative impact, remove brackets from labels
breaks <- seq(min(final_results_long$Value, na.rm = TRUE),
              max(final_results_long$Value, na.rm = TRUE),
              length.out = 7)

# Define the custom labels for bins, removing the brackets
labels <- paste0(round(breaks[-length(breaks)] + 0.01, 1), "–", round(breaks[-1] - 0.01, 1))

# Apply the cut with custom labels
final_results_long$impact_bin <- cut(final_results_long$Value,
                                     breaks = breaks,
                                     include.lowest = TRUE,
                                     right = TRUE,
                                     labels = labels)

# Ensure the factor is clean
final_results_long$impact_bin <- factor(final_results_long$impact_bin)

# Use scico palette with 6 colors
fill_colors <- c("#540b0e", "#9e2a2b", "#e09f3e", "#fff3b0", "#a3b18a", "#3a5a40")
names(fill_colors) <- levels(final_results_long$impact_bin)

# Create map without legend
map_SM <- ggplot(final_results_long) +
  geom_sf(aes(geometry = hex_geometry, fill = Value ), color = NA) +
  geom_polygon(data = countries, aes(x = long, y = lat, group = group), fill = "gray75") +
  coord_sf(xlim = c(-180, 60), ylim = c(20, 85), expand = FALSE) +
  # scico::scale_fill_scico(palette = "buda", name = "Prop of safe", limits = c(0, 1)) +
  scale_fill_gradientn(
    colours = fill_colors,
    limits = c(0, 1)) +
  theme_minimal() +
  labs(x = "", y = "", title = "(a) Spatial distributions of ESMs") +
  theme(legend.position = "none")

# Create histogram legend
hist_legend_SM <- ggplot(final_results_long, aes(x = impact_bin)) +
  geom_bar(aes(fill = impact_bin), color = "white", size = 0.01) +
  scale_fill_manual(values = fill_colors, guide = "none") +
  theme_void() +
  scale_x_discrete(labels = labels) +  # Custom labels without brackets
  labs(x = "Unsafe   -   Safe") +
  theme(axis.text.x = element_text(angle = 30, size = 4),
        axis.title.x = element_text(size = 5))

# Combine
map_SM + inset_element(hist_legend_SM, left = 0.02, bottom = 0.01, right = 0.2, top = 0.45)


#### Null model ####
all_ranges <- filtered_data2[,c(10,25,26,29,30,33,34)]
all_ranges <- st_drop_geometry(all_ranges)
all_ranges <- all_ranges %>%
  distinct(Species, .keep_all = TRUE)
all_ranges <- na.omit(all_ranges)
all_ranges2 <- as.data.frame(all_ranges)
rownames(all_ranges2) <- all_ranges2[,1]
all_ranges2 <- all_ranges2[,-1]

n_model_list <- list()

# this is an example with 3 permutations
for (j in 1:3){

  sub <- subset(filtered_data2, filtered_data2$haul_id %in% cell_names[,j])
  sub <- st_drop_geometry(sub)
  biommatrix <- sub[,c("Species","wgt_cpua","cell_decade")]
  
  nullmetrics <- data.frame()
  
  for (i in 1:3){
    print(i)
  randm <- randomizeMatrix(all_ranges2, iterations = 1000) #creates randomized occurrence matrix
  traitmatrix <- as.data.frame(randm[rownames(randm) %in% unique(biommatrix$Species),])
  traitmatrix$Species <- rownames(traitmatrix)
  biommatrix <-  biommatrix[biommatrix$Species %in% traitmatrix$Species,]
  
  cwm_values <- biommatrix %>%
    left_join(traitmatrix, by = "Species") %>%
    group_by(cell_decade)  %>%
    summarise(across(where(is.numeric), ~ sum(.x * wgt_cpua, na.rm = TRUE) / sum(wgt_cpua, na.rm = TRUE), .names = "CWM_{.col}"), 
              .groups = "drop")
  
  sub2 <- as.data.frame(sub[, c("max_temp_dec", "min_temp_dec", "min_sal_dec", "min_o2_dec", 
                                "cell_decade", "hex_geometry", "decade")])
  sub2 <- sub2[!duplicated(sub2[, -which(names(sub2) == "hex_geometry")]), ]
  sub2 <- left_join(cwm_values, sub2, by = "cell_decade")
  
  # Calculate metrics
  metrics <- sub2 %>%
    dplyr::group_by(cell_decade) %>%
    dplyr::mutate(# using pref trait tolerance values
      usm_temp_pref = CWM_TempPrefMax - max_temp_dec,
      lsm_temp_pref = min_temp_dec - CWM_TempPrefMin,
      lsm_sal_pref = min_sal_dec - CWM_SalinityPrefMin,
      lsm_o2_pref = min_o2_dec - CWM_OxyPrefMin)
  
  metrics$ite <- i
  nullmetrics <- rbind(nullmetrics, metrics)
  }

  # Store result in the list
  n_model_list[[j]] <- nullmetrics
}

n_model_list_df <- n_model_list

# save(n_model_list_df, file = "Data/nullmodel_sf.rData")
load(file = "Data/nullmodel_sf.rData") # note: this is an example only, need to increase permutation n

for(i in seq_along(n_model_list_df)){
  n_model_list_df[[i]] <- as.data.frame(n_model_list_df[[i]])
  n_model_list_df[[i]]$lap <- i
}

n_model_list_df <- do.call(rbind.data.frame, n_model_list_df)

# mean of all metrics after randomization
n_model_final <- n_model_list_df %>%
  dplyr::group_by(hex_geometry, cell_decade) %>%
  dplyr::summarise(# preferred tolerances
    usm_temp_mean_pref = mean(usm_temp_pref, na.rm = TRUE),
    lsm_temp_mean_pref = mean(lsm_temp_pref, na.rm = TRUE),
    lsm_sal_mean_pref = mean(lsm_sal_pref, na.rm = TRUE),
    lsm_o2_mean_pref = mean(lsm_o2_pref, na.rm = TRUE),
    standev_usm = sd(usm_temp_pref),
    standev_lsm = sd(lsm_temp_pref),
    standev_lsm = sd(lsm_sal_pref),
    standev_lsm = sd(lsm_o2_pref))

n_model_final <- n_model_final %>%
  separate(cell_decade, into = c("cell", "start", "end"), sep = "_") %>%
  mutate(decade = paste0(start, "_", end)) %>%
  select(cell, decade, everything(), -start, -end)

# save(n_model_final, file = "Data/nullmodel_cell_decade_final.rData")

# Calculating SES
# (observed – mean of null distribution/standard deviation (SD)
obs <- final_results2

ses_df <- n_model_final %>%
  left_join(obs, by = "hex_geometry", suffix = c("_null", "_obs")) %>%
  group_by(cell, decade) %>%
  mutate(
    ses_usm_temp = (usm_temp_mean_pref_obs - usm_temp_mean_pref_null) / standev_usm,
    ses_lsm_temp = (lsm_temp_mean_pref_obs - lsm_temp_mean_pref_null) / standev_lsm,
    ses_lsm_sal  = (lsm_sal_mean_pref_obs  - lsm_sal_mean_pref_null)  / standev_lsm,
    ses_lsm_o2   = (lsm_o2_mean_pref_obs   - lsm_o2_mean_pref_null)   / standev_lsm
  ) %>%
  select(hex_geometry, starts_with("ses_"))

# plot SES
# Pivot SES columns to long format
ses_long <- ses_df %>%
  pivot_longer(
    cols = starts_with("ses_"),
    names_to = "variable",
    values_to = "SES")

ses_long <- ses_long %>%
  mutate(
    variable = variable %>%
      str_remove("^ses_") %>%               # Remove "ses_" prefix
      str_replace_all("_", " ") %>%         # Replace underscores with spaces
      str_to_title())
    
# open Spalding ecoregions
ecoregions <- st_read("Data/meow_ecos.shp")
ecoregion_coords <- ecoregions %>%
  mutate(ecoregion_geometry = geometry) %>%
  st_drop_geometry()

# merge to null safety margins
ecoregion_coords <- ecoregion_coords[,c(1,10)]
ses_long <- st_as_sf(ses_long)
polyg <- st_transform(ses_long, crs = st_crs(ecoregions))

# find intersections: polygons that overlay ecoregions
pol_ecoreg <- st_join(polyg, ecoregions, join = st_nearest_feature)
pol_ecoreg <- pol_ecoreg %>%
  left_join(ecoregion_coords, by = "ECO_CODE")

pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Mediterranean Sea", 
                              "Lusitanian", pol_ecoreg$PROVINCE)
pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Tropical Northwestern Atlantic", 
                              "Warm Temperate Northwest Atlantic", pol_ecoreg$PROVINCE)
pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Warm Temperate Northeast Pacific", 
                              "Cold Temperate Northeast Pacific", pol_ecoreg$PROVINCE)

pol_ecoreg <- pol_ecoreg %>%
  mutate(facet_label = case_when(
    PROVINCE == "Lusitanian"                        ~ "Lusitanian",
    PROVINCE == "Arctic"                            ~ "Arctic",
    PROVINCE == "Cold Temperate Northeast Pacific"  ~ "Cold Temp.\nNE Pacific",
    PROVINCE == "Northern European Seas"            ~ "Northern\nEuropean Seas",
    PROVINCE == "Warm Temperate Northwest Atlantic" ~ "Warm Temp.\nNW Atlantic",
    PROVINCE == "Cold Temperate Northwest Atlantic" ~ "Cold Temp.\nNW Atlantic",
    TRUE ~ PROVINCE))

ses_summary <- pol_ecoreg %>%
  group_by(variable, facet_label) %>%
  summarise(
    mean_SES = mean(SES, na.rm = TRUE),
    se = sd(SES, na.rm = TRUE) / sqrt(n()),  # standard error
    lower = mean_SES - 1.96 * se,
    upper = mean_SES + 1.96 * se,
    .groups = "drop") %>%
  mutate(
    signif_type = case_when(
      lower > 0 & upper > 0 ~ "positive",
      lower < 0 & upper < 0 ~ "negative",
      TRUE ~ "ns"))

# "Boxplot" of SES
ggplot(ses_summary, aes(y = variable, x = mean_SES, color = signif_type)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, color = "black") +
  facet_wrap(~facet_label, ncol = 6) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(
    values = c(
      positive = "cyan3",      # significant positive
      negative = "orangered4", # significant negative
      ns = "grey70"            # not significant
    ),
    name = "Significance"
  ) +
  theme_minimal() +
  labs(
    title = "Standardized Effect Size \nper Safety Margin",
    x = "Effect size",
    y = NULL
  )

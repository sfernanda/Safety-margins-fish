##### VII - Spider plot for Spalding Provinces/ecoregions #####

library(patchwork)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyr)
library(purrr)
library(scales)

# open Spalding ecoregions
ecoregions <- st_read("Data/meow_ecos.shp")
ecoregion_coords <- ecoregions %>%
  mutate(ecoregion_geometry = geometry) %>%
  st_drop_geometry()
ecoregion_coords <- ecoregion_coords[,c(1,10)]

# open results (metrics after permutations )
load(file = "Data/resultsbiomass_list_cwm_230cells.rData")

# Create a list of iteration numbers corresponding to each element in results_list
iterations <- seq_along(results_list)

# Add iteration number to each data frame in the list
results_list_with_iterations <- map2(results_list, 99,
                                     ~mutate(.x, iteration = .y))

# Combine all results into one data frame
final_results <- bind_rows(results_list_with_iterations)

# mean of all metrics after randomization
final_results2 <- final_results %>%
  dplyr::group_by(decade, hex_geometry) %>%
  dplyr::summarise(
    usm_temp_mean = mean(usm_temp, na.rm = TRUE),
    lsm_temp_mean = mean(lsm_temp, na.rm = TRUE),
    lsm_sal_mean = mean(lsm_sal, na.rm = TRUE),
    lsm_o2_mean = mean(lsm_o2, na.rm = TRUE),
   
  )


final_results2 <- st_as_sf(final_results2)
polyg <- st_transform(final_results2, crs = st_crs(ecoregions))

# find intersections: polygons that overlay ecoregions
pol_ecoreg <- st_join(polyg, ecoregions, join = st_nearest_feature)
pol_ecoreg <- pol_ecoreg %>%
  left_join(ecoregion_coords, by = "ECO_CODE")

# function for spider plot using ggplot
custom_colors <- rev(c("#B84193", "#BF812D", "darkcyan"))
ggspider <- function(p_data, 
                          legend_title = "Group", 
                          background_color = "gray97", 
                          area_fill = TRUE, 
                          central_distance = 0.2, 
                          axis_name_offset = 0.2) {
  
  # Function to generate circular coordinates for each axis
  circle_coords <- function(r, n_axis = ncol(p_data) - 2){
    fi <- seq(0, 2*pi, length.out = n_axis + 1) + pi/2
    tibble(x = r * cos(fi), y = r * sin(fi), r)
  }
  
  # Background grid
  step_1 <- map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
    ggplot(aes(x, y)) +
    geom_polygon(data = circle_coords(1 + central_distance), alpha = 1, fill = background_color) +
    geom_path(aes(group = r), lty = 2, alpha = 0.1) +
    theme_void()

  polygon_sf <- map_df(1+central_distance, circle_coords) %>%
    summarise(geometry = st_sfc(st_polygon(list(
      rbind(cbind(x, y), c(x[1], y[1]))  # Ensure closure
    ))), .groups = "drop") %>%
    st_as_sf()
  
  # Calculate total area
  tot_polygon <- st_area(polygon_sf)
  
   # Define axis lines from center to 1
  axis_coords <- tibble(
    x = cos(seq(0, 2*pi, length.out = 2) + pi/2),
    y = sin(seq(0, 2*pi, length.out = 2) + pi/2),
    xend = x * (1 + central_distance),
    yend = y * (1 + central_distance)
  )
  
  text_data <- p_data %>%
    ungroup() %>%
    select(-group,-decade) %>%
    map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 0.25)) %>%
    mutate(r = seq(0, 1, 0.25)) %>%
    pivot_longer(-r, names_to = "parameter", values_to = "value")
  
  # Generate label positions along axes
  text_coords <- function(r, n_axis = ncol(p_data) - 2) {
    fi <- seq(0, 2 * pi, length.out = n_axis + 1)[- (n_axis + 1)] + pi / 2
    tibble(x = r * cos(fi), y = r * sin(fi), r = r - central_distance)
  }
  
  labels_data <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
    bind_cols(text_data %>% select(-r))
  
  label_var <- labels_data %>% 
    filter(r == 1) 
  
  label_scale <- labels_data %>%
    group_by(r) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  # rescaling data from 0 to 1
  rescaled_coords <- function(r, n_axis){
    fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
    tibble(r, fi) %>% mutate(x = r*cos(fi), y = r*sin(fi)) %>% select(-fi)
  }
  
  rescaled_data <- p_data %>% 
    ungroup() %>%
    mutate(across(-c(group,decade), ~ scales::rescale(.))) %>%
    mutate(copy = pull(., 3)) %>%
    pivot_longer(-c(group, decade), names_to = "parameter", values_to = "value") %>%
    group_by(group,decade) %>%
    mutate(coords = rescaled_coords(value + central_distance, ncol(p_data) - 2)) %>%
    unnest(cols = c(coords)) 
  
  # Create the plot
  plot <- step_1 + 
    geom_segment(data = axis_coords, aes(x = 0, y = 0, xend = xend, yend = yend), 
                 color = "black", size = 0.5) +
    geom_point(data = rescaled_data, aes(x, y, group = decade, col = decade), size = 3) +
    facet_wrap(~group) +
    geom_path(data = rescaled_data, aes(x, y, group = decade, col = decade), size = 1) +
    {if (area_fill) geom_polygon(data = rescaled_data, aes(x, y, group = decade, col = decade, fill = decade), 
                                 size = 1, alpha = 0.2, show.legend = FALSE)} +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    geom_text(data = label_scale, 
              aes(x, y, label = r), alpha = 0.65, size = 2.5, col = "black",
              nudge_x = 0.1) +
    geom_text(data = label_var, 
              aes(x, y, label = parameter), alpha = 0.65, size = 3, 
              nudge_x = label_var$x * 0.2, # Move further along x and y-axis
              nudge_y = label_var$y * 0.2) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  return(list(plot = plot, rescaled_data = rescaled_data, total_area = tot_polygon))
}


# for upper and lower safety margins (extreme values)
pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Mediterranean Sea", 
                              "Lusitanian", pol_ecoreg$PROVINCE)
pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Tropical Northwestern Atlantic", 
                              "Warm Temperate Northwest Atlantic", pol_ecoreg$PROVINCE)
pol_ecoreg$PROVINCE <- ifelse(pol_ecoreg$PROVINCE == "Warm Temperate Northeast Pacific", 
                              "Cold Temperate Northeast Pacific", pol_ecoreg$PROVINCE)

# plot spalding areas
pol_ecoreg <- st_as_sf(pol_ecoreg)
pol_ecoreg <- st_wrap_dateline(pol_ecoreg, options = c("WRAPDATELINE=YES","DATELINEOFFSET=180"), quiet = TRUE)

world_map <- map_data("world")

# world map as a background for fig 1.c
library(MetBrewer)
provmap <- pol_ecoreg %>%
  distinct(hex_geometry, .keep_all = TRUE)
ggplot(provmap) +
  geom_sf(aes(geometry = hex_geometry, fill = PROVINCE), alpha = .4, color = NA) +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "gray75") +
  coord_sf(xlim = c(-180, 60), ylim = c(20, 85), expand = FALSE) +
  theme_minimal() +
  scale_fill_manual(values = met.brewer("Manet", length(unique(pol_ecoreg$PROVINCE)))) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  labs(y="", x="")


# apply spider plot function
p_data1 <- pol_ecoreg %>%
  dplyr::group_by(decade, PROVINCE) %>%
  dplyr::summarise(across(where(is.numeric), function(x) mean(x, na.rm = TRUE)))
p_data2 <- p_data1[,c(1:6)]
p_data2 <- st_drop_geometry(p_data2)
names(p_data2)[2] <- "group"
p_data2 <- p_data2[,c(1,2,3,5,4,6)] # sort in the correct order for plot
colnames(p_data2) <- gsub("_mean", "", colnames(p_data2))
p_data2 <- p_data2 %>%
  mutate(across(where(is.numeric), ~ rescale(.)))

spider_safmar <- ggspider(p_data2)

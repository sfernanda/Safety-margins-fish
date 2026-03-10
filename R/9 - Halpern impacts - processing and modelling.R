##### IX - Processing Halpern impacts and modelling ##### 

library(tibble)
library(sf)
library(ggplot2)
library(patchwork)
library(dplyr)
library(exactextractr)
library(mgcv)
library(tidyr)
library(purrr)

##### Relevant info ####
# The impacts are divided in 4 categories: 
# •	Fishing: commercial demersal destructive, commercial demersal nondestructive high bycatch, commercial demersal nondestructive low bycatch, pelagic high bycatch, pelagic low bycatch, artisanal
# •	Ocean: shipping
# •	Land-based: nutrient pollution, organic chemical pollution, direct human, light
#
# Averages and slopes are extracted for models:
# Only fishing is extracted individually - average within each grid cell per year (then tot average across years and slope)
# All other groups are cumulative (sum) per yr - then avg across all years within each cell and slope of the trend 

##### 1. Load data: per groups and fish individually ####
load("Data/fishing_hex_df.rData")
load("Data/ocean_hex_df.rData")
load("Data/landbased_hex_df.rData")

load(file = "DataBases/fishing_halpern_means.rData") # individual fishing types (need to calculate slope)
load(file = "DataBases/all_imp_halpern_means.rData") # cumulative impacts without sst
load(file = "DataBases/all_imp_halpern_means_sst.rData") # cumulative impacts with sst


# reshape the categories df
reshape_hex_df <- function(df, prefix) {
  df %>%
    mutate(year_col = paste0(prefix, "_", year, "_impact")) %>%
    dplyr::select(cell, year_col, mean) %>%
    pivot_wider(names_from = year_col, values_from = mean)
}

fishing_wide     <- reshape_hex_df(fishing_hex_df, "fishing_tot")
ocean_wide       <- reshape_hex_df(ocean_hex_df, "ocean")
landbased_wide   <- reshape_hex_df(landbased_hex_df, "landbased")

# Function to extract slopes for fishing individual
function_slopes <- function(df) {
  # Get the years from the column names
  year_cols <- grep("_impact$", names(df), value = TRUE)
  years <- as.numeric(gsub(".*_(\\d{4})_impact$", "\\1", year_cols))
  
  # Select only the impact columns
  impacts <- df %>% dplyr::select(all_of(year_cols))
  impacts <- st_drop_geometry(impacts)
  
  # For each row (cell), fit lm(impact ~ year) and extract slope
  slopes <- apply(impacts, 1, function(x) {
    if (all(is.na(x))) return(NA)
    fit <- lm(x ~ years)
    coef(fit)[2]  # slope
  })
  
  df$slope_per_type <- slopes
  return(df)
}

fishing_wide <- function_slopes(fishing_wide)
names(fishing_wide)[14] <- "slope_fishing_tot"
fishing_wide <- fishing_wide %>%
  mutate(mean_fishing_tot = rowMeans(across(contains("_impact") & where(is.numeric)), na.rm = TRUE))

ocean_wide <- function_slopes(ocean_wide)
names(ocean_wide)[14] <- "slope_ocean_tot"
ocean_wide <- ocean_wide %>%
  mutate(mean_ocean_tot = rowMeans(across(contains("_impact") & where(is.numeric)), na.rm = TRUE))

landbased_wide <- function_slopes(landbased_wide)
names(landbased_wide)[14] <- "slope_landbased_tot"
landbased_wide <- landbased_wide %>%
  mutate(mean_landbased_tot = rowMeans(across(contains("_impact") & where(is.numeric)), na.rm = TRUE))

##### 2. Merge all impact data by 'cell' #### 
df_list <- list(
  fishing_wide %>% dplyr::select(cell, geometry, mean_fishing_tot, slope_fishing_tot),
  ocean_wide %>% dplyr::select(cell, geometry, mean_ocean_tot, slope_ocean_tot),
  landbased_wide %>% dplyr::select(cell, geometry, mean_landbased_tot, slope_landbased_tot))

halpern_impacts <- reduce(df_list, function(x, y) {
  inner_join(sf::st_drop_geometry(x),
             sf::st_drop_geometry(y),
             by = "cell")
}) # to avoid error due to df geometries

# add geometry back
halpern_impacts <- left_join(
  df_list[[1]] %>% dplyr::select(cell, geometry),
  halpern_impacts,
  by = "cell"
) %>% sf::st_as_sf()

##### 3. Open safety margins metrics after permutations to calculate cumulative safety ####
load(file = "Data/resultsbiomass_list_cwm_230cells.rData")

# Combine all results into one data frame
final_results <- bind_rows(results_list)

# Subset rows where 'cell_decade' contains '2003_2012'
final_results <- final_results %>%
  filter(grepl("2003_2012", cell_decade))

# mean of all metrics after randomization
final_results2 <- final_results %>%
  dplyr::group_by(cell_decade, hex_geometry) %>%
  dplyr::summarise(
    # preferred tolerances
    usm_temp_mean_pref = mean(usm_temp_pref, na.rm = TRUE),
    lsm_temp_mean_pref = mean(lsm_temp_pref, na.rm = TRUE),
    lsm_sal_mean_pref = mean(lsm_sal_pref, na.rm = TRUE),
    lsm_o2_mean_pref = mean(lsm_o2_pref, na.rm = TRUE),
    .groups = "drop")
str(final_results2)

final_results2 <- final_results2 %>%
  dplyr::mutate(cell = as.numeric(sub("_.*", "", cell_decade)))

# rescale safety margin results in a scale 0-1
final_results3 <- final_results2 %>%
  mutate(across(c(usm_temp_mean_pref, lsm_temp_mean_pref,
                  lsm_sal_mean_pref, lsm_o2_mean_pref), ~ scales::rescale(.)))
names(final_results3)[2] <- "geometry"

final_results3 <- final_results3 %>%
  group_by(geometry, cell) %>%
  dplyr::summarise(Safety = usm_temp_mean_pref + lsm_temp_mean_pref + lsm_sal_mean_pref + lsm_o2_mean_pref)

##### 4. Combine response var (Safety) with all Halpern impacts for modelling  ####
data_mod <- final_results3 %>%
  left_join(halpern_impacts, by = "cell") # note 2 types of geometry - 2 different projections
data_mod <- na.omit(data_mod)

# log transform explanatory variables
data_mod$mean_fishing_tot <- log(data_mod$mean_fishing_tot + 1)
data_mod$mean_ocean_tot <- log(data_mod$mean_ocean_tot + 1)
data_mod$mean_landbased_tot <- log(data_mod$mean_landbased_tot + 1)

# Function to fit GAM/LM with Safety ~ s(selected_variable, k=3)
fit_significant_lms_with_plot <- function(data, vars) {
  significant_models <- list()
  plot_list <- list()
  model_summaries <- list()
  residuals_checks <- list()
  
  for (var in vars) {
    # Build the formula
    form <- as.formula(paste("Safety ~", var))
    
    # Fit the linear model safely
    model <- tryCatch(lm(form, data = data), error = function(e) NULL)
    
    # Check significance of the slope
    if (!is.null(model)) {
      pval <- summary(model)$coefficients[2, "Pr(>|t|)"]
      
      if (!is.na(pval) && pval < 0.05) {
        significant_models[[var]] <- model
        
        # Summary of model
        model_summaries[[var]] <- capture.output(summary(model))
        
        # Residuals diagnostic
        residuals_data <- residuals(model)
        shapiro_test <- tryCatch(shapiro.test(residuals_data), error = function(e) NULL)
        residuals_checks[[var]] <- list(
          shapiro_p = if (!is.null(shapiro_test)) shapiro_test$p.value else NA,
          histogram = ggplot(data.frame(residuals = residuals_data), aes(x = residuals)) +
            geom_histogram(fill = "#69b3a2", color = "white") +
            theme_minimal() +
            labs(title = var)
        )
        
        # Clean data and plot
        clean_data <- data[!is.na(data[[var]]) & !is.na(data$Safety), ]
        
        p <- ggplot(clean_data, aes(x = .data[[var]], y = Safety)) +
          geom_point(alpha = 0.2) +
          geom_smooth(method = "lm", se = TRUE, color = "gray80") +
          theme_minimal() +
          labs(y = "Cumulative Safety Margin")
        
        plot_list[[var]] <- p
      }
    }
  }
  
  list(
    models = significant_models,
    plots = plot_list,
    summaries = model_summaries,
    residuals = residuals_checks
  )
}

variables <- names(data_mod)[c(5,7,9)] # select only fishing, landbased, and mean ocean (shipping)
signif_lms <- fit_significant_lms_with_plot(data_mod, variables)

# see which variables had significant models:
names(signif_lms$models)
ggpubr::ggarrange(signif_lms$plots[[1]],
                  signif_lms$plots[[2]],
                  signif_lms$plots[[3]])

hist_plots <- lapply(signif_lms$residuals, function(x) x$histogram)
ggpubr::ggarrange(plotlist = hist_plots, ncol = 2, nrow = ceiling(length(hist_plots)/2))

mod_pressures1 <- ggplot(data_mod, aes(x = mean_fishing_tot, y = Safety)) +
  geom_point(alpha = 0.5, col = "gray60", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "gray60", fill = "gray60") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80")
  ) +
  labs(
    x = "Log Fishing",
    y = "Cumulative Safety Margin")

mod_pressures2 <- ggplot(data_mod, aes(x = mean_ocean_tot, y = Safety)) +
  geom_point(alpha = 0.5, col = "gray60", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "gray60", fill = "gray60") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80")
  ) +
  labs(
    x = "Log Shipping",
    y = "Cumulative Safety Margin")

mod_pressures3 <- ggplot(data_mod, aes(x = mean_landbased_tot, y = Safety)) +
  geom_point(alpha = 0.5, col = "gray60", size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "gray60", fill = "gray60") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray80")
  ) +
  labs(
    x = "Log Land-based",
    y = "Cumulative Safety Margin")

ggpubr::ggarrange(mod_pressures1,mod_pressures2,mod_pressures3, ncol = 3)

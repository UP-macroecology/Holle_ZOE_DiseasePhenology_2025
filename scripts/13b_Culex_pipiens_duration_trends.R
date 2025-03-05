# ZOE project 
# Disease phenology analysis of Culex pipiens in Europe (primary transmitter of WNV)

# ---------------------------------------------------------------------- #
#       10b. Past and future duration trends of disease phenology        #
# ---------------------------------------------------------------------- #


# Load the needed packages
library(terra)
library(tidyverse)
library(ggplot2)
library(sf)


#-------------------------------------------------------------------------------

# 1. Past duration trends based on disease phenology ---------------------------

# Read in the monthly binary prediction data from 1970 to 2019
r_curr_preds_clim_landuse_ens_bin <- terra::rast("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif")

# Extract the years from the layer names
r_names_curr <- names(r_curr_preds_clim_landuse_ens_bin)
r_years_curr <- as.numeric(substr(r_names_curr, 4, 7))

# Calculate the number of months with a predicted presence per year (for each cell)
length_presence_year_curr <- terra::tapp(r_curr_preds_clim_landuse_ens_bin, index = r_years_curr, fun = sum, na.rm = TRUE)

# Rename yearly layers
names(length_presence_year_curr) <- paste0(unique(r_years_curr))

# Create a vector with the studied decades
r_names_curr <- names(length_presence_year_curr)
r_decades_curr <- floor(as.numeric(r_names_curr) / 10) * 10

# Compute the mean presence length per decade (per cell)
decadal_means_curr <- tapp(length_presence_year_curr, index = r_decades_curr, fun = mean, na.rm = TRUE)

# Rename decadal layers
names(decadal_means_curr) <- paste0(unique(r_decades_curr))

# Function to compute slope (trend) from 1970s to 2010s for each cell
calc_slope_curr <- function(x) {
  if (all(is.na(x))) return(NA)  # Handle NA pixels
  decade_values <- seq(1970, 2010, by = 10)  # Decade midpoint years
  lm_fit <- lm(x ~ decade_values)  # Fit linear model
  return(coef(lm_fit)[2])  # Extract slope coefficient
}

# Apply function to raster to calculate trend (slope)
slope_raster_curr <- app(decadal_means_curr, fun = calc_slope_curr)
names(slope_raster_curr) <- "Trend_Slope"

# Convert raster to a data frame for plotting
slope_df_curr <- as.data.frame(slope_raster_curr, xy = TRUE, na.rm = TRUE)
colnames(slope_df_curr) <- c("lon", "lat", "trend")

# Find out the highest and lowest slope
print(max(slope_df_curr$trend))
print(min(slope_df_curr$trend))  

# Plot
ggplot(slope_df_curr) +
  geom_tile(aes(x = lon, y = lat, fill = trend)) +
  scale_fill_gradient2(
    low = "midnightblue", mid = "grey93", high = "firebrick4",
    midpoint = 0, name = "Trend / Slope",
    limits = c(-0.35, 0.36) 
  ) +
  theme_minimal() +
  labs(title = "Culex pipiens - Trend in Presence Length (1970s-2010s)",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(face = "bold"))

ggsave("output_data/plots/duration_trends/C_pipiens_duration_trends_past.png", width = 8.5, height = 6)





#-------------------------------------------------------------------------------

# 2. Future duration trends based on disease phenology -------------------------

# Create a vector containing the three different environmental change scenarios
env_scenarios <- c("ssp126", "ssp370", "ssp585")



for (s in env_scenarios) { # Start the loop over the three environmental scenarios
  
  print(s)
  
  
  # Read in the monthly binary prediction data from 1970 to 2019
  r_fut_preds_clim_landuse_ens_bin <- terra::rast(paste0("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_",s,".tif"))
  
  # Extract the years from the layer names
  r_names_fut <- names(r_fut_preds_clim_landuse_ens_bin)
  r_years_fut <- as.numeric(substr(r_names_fut, 4, 7))
  
  # Calculate the number of months with a predicted presence per year (for each cell)
  length_presence_year_fut <- terra::tapp(r_fut_preds_clim_landuse_ens_bin, index = r_years_fut, fun = sum, na.rm = TRUE)
  
  # Rename yearly layers
  names(length_presence_year_fut) <- paste0(unique(r_years_fut))
  
  # Create a vector with the studied decades
  r_names_fut <- names(length_presence_year_fut)
  r_decades_fut <- floor(as.numeric(r_names_fut) / 10) * 10
  
  # Compute the mean presence length per decade (per cell)
  decadal_means_fut <- tapp(length_presence_year_fut, index = r_decades_fut, fun = mean, na.rm = TRUE)
  
  # Rename decadal layers
  names(decadal_means_fut) <- paste0(unique(r_decades_fut))
  
  decadal_means_curr_2010 <- subset(decadal_means_curr, "2010")
  decadal_means_fut_2030_2050 <- subset(decadal_means_fut, c("2030", "2040", "2050"))
  
  decadal_means_fut <- c(decadal_means_curr_2010, decadal_means_fut_2030_2050)
  
  
  # Function to compute slope (trend) from 1970s to 2010s for each cell
  calc_slope_fut <- function(x) {
    if (all(is.na(x))) return(NA)  # Handle NA pixels
    decade_values <- c(2010, 2030, 2040, 2050)  # Decade midpoint years
    lm_fit <- lm(x ~ decade_values)  # Fit linear model
    return(coef(lm_fit)[2])  # Extract slope coefficient
  }
  
  
  # Apply function to raster to calculate trend (slope)
  slope_raster_fut <- app(decadal_means_fut, fun = calc_slope_fut)
  names(slope_raster_fut) <- "Trend_Slope"
  
  # Convert raster to a data frame for plotting
  slope_df_fut <- as.data.frame(slope_raster_fut, xy = TRUE, na.rm = TRUE)
  colnames(slope_df_fut) <- c("lon", "lat", "trend")
  
  # Find out the highest and lowest slope
  print(max(slope_df_fut$trend))
  print(min(slope_df_fut$trend))
  
  # Plot
  ggplot(slope_df_fut) +
    geom_tile(aes(x = lon, y = lat, fill = trend)) +
    scale_fill_gradient2(
      low = "midnightblue", mid = "grey93", high = "firebrick4",
      midpoint = 0, name = "Trend / Slope",
      limits = c(-0.33, 0.36) 
    ) +
    theme_minimal() +
    labs(title = paste0("Culex pipiens - Trend in Presence Length (2010s-2050s,", s,")"),
         x = "Longitude", y = "Latitude") +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(paste0("output_data/plots/duration_trends/C_pipiens_duration_trends_future_",s,".png"), width = 8.5, height = 6)
  
  
} # Close the loop over the three environmental scenarios




#-------------------------------------------------------------------------------

# 3. Create a tile plot to visualize shifts in duration length over decades ----

# Read in the monthly binary prediction data from 1970 to 2019
r_curr_preds_clim_landuse_ens_bin <- terra::rast("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif")
r_fut_preds_clim_landuse_ens_bin_126 <- terra::rast(paste0("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp126.tif"))
r_fut_preds_clim_landuse_ens_bin_370 <- terra::rast(paste0("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp370.tif"))
r_fut_preds_clim_landuse_ens_bin_585 <- terra::rast(paste0("output_data/results/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp585.tif"))

# Extract layer names
rasters_names_curr <- names(r_curr_preds_clim_landuse_ens_bin)
rasters_names_fut <- names(r_fut_preds_clim_landuse_ens_bin_126)

# Extract years from layer names
rasters_years_curr <- as.numeric(sub(".*/", "", rasters_names_curr))
rasters_years_fut <- as.numeric(sub(".*/", "", rasters_names_fut))

# Extract layers by decade
rasters_1970s <- which(rasters_years_curr >= 1970 & rasters_years_curr < 1980)
rasters_1970s <- r_curr_preds_clim_landuse_ens_bin[[rasters_1970s]]
rasters_2010s <- which(rasters_years_curr >= 2010 & rasters_years_curr < 2020)
rasters_2010s <- r_curr_preds_clim_landuse_ens_bin[[rasters_2010s]]
rasters_2050s_126 <- which(rasters_years_fut >= 2050 & rasters_years_fut < 2060)
rasters_2050s_126 <- r_fut_preds_clim_landuse_ens_bin_126[[rasters_2050s_126]]
rasters_2050s_370 <- which(rasters_years_fut >= 2050 & rasters_years_fut < 2060)
rasters_2050s_370 <- r_fut_preds_clim_landuse_ens_bin_370[[rasters_2050s_370]]
rasters_2050s_585 <- which(rasters_years_fut >= 2050 & rasters_years_fut < 2060)
rasters_2050s_585 <- r_fut_preds_clim_landuse_ens_bin_585[[rasters_2050s_585]]

# Extract the months from layer names
months_1970s <- as.numeric(sub("^([0-9]{2})/([0-9]{4})$", "\\1", names(rasters_1970s)))
months_2010s <- as.numeric(sub("^([0-9]{2})/([0-9]{4})$", "\\1", names(rasters_2010s)))  
months_2050s_126 <- as.numeric(sub("^([0-9]{2})/([0-9]{4})$", "\\1", names(rasters_2050s_126)))  
months_2050s_370 <- as.numeric(sub("^([0-9]{2})/([0-9]{4})$", "\\1", names(rasters_2050s_370))) 
months_2050s_585 <- as.numeric(sub("^([0-9]{2})/([0-9]{4})$", "\\1", names(rasters_2050s_585))) 

# Calculate the mean number of cells with monthly presences per decade
mean_presence_by_month <- function(r_stack, months) {
  tapply(1:nlyr(r_stack), months, function(i) {
    # Compute the total presence per layer (global() returns a data frame)
    presence_counts <- sapply(i, function(layer) global(r_stack[[layer]], fun = "sum", na.rm = TRUE)[1, 1])
    mean(presence_counts, na.rm = TRUE)  # Take the mean of presence counts for the month
  })
}

mean_presences_1970s <- mean_presence_by_month(rasters_1970s, months_1970s)
mean_presences_2010s <- mean_presence_by_month(rasters_2010s, months_2010s)
mean_presences_2050s_126 <- mean_presence_by_month(rasters_2050s_126, months_2050s_126)
mean_presences_2050s_370 <- mean_presence_by_month(rasters_2050s_370, months_2050s_370)
mean_presences_2050s_585 <- mean_presence_by_month(rasters_2050s_585, months_2050s_585)

df_mean_presences_1970s <- data.frame(month = 1:12, mean_presence = mean_presences_1970s, decade = "1970s")
df_mean_presences_2010s <- data.frame(month = 1:12, mean_presence = mean_presences_2010s, decade = "2010s")
df_mean_presences_2050s_126 <- data.frame(month = 1:12, mean_presence = mean_presences_2050s_126, decade = "2050s; ssp126")
df_mean_presences_2050s_370 <- data.frame(month = 1:12, mean_presence = mean_presences_2050s_370, decade = "2050s; ssp370")
df_mean_presences_2050s_585 <- data.frame(month = 1:12, mean_presence = mean_presences_2050s_585, decade = "2050s; ssp585")

# Bind the results data frame of the different decades
df_mean_presences <- rbind(df_mean_presences_1970s, df_mean_presences_2010s, df_mean_presences_2050s_126,
                           df_mean_presences_2050s_370, df_mean_presences_2050s_585)

# Plot
ggplot(df_mean_presences, aes(x = month, y = factor(decade), fill = mean_presence)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma", name = "Mean count of\npresence cells") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Decadal trends in the activity season duration of Culex pipiens across Europe",
       x = "Month", y = "Decade") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("output_data/plots/duration_trends/C_pipiens_duration_trends_activityseason.png", width = 8.5, height = 3.5)

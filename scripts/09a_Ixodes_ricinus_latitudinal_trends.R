# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

# ---------------------------------------------------------------------- #
#     08a. Past and future latitudinal trends of disease phenology       #
# ---------------------------------------------------------------------- #


library(ggplot2)
library(terra)
library(tidyverse)
library(scico)


#-------------------------------------------------------------------------------

# 1. Past and future latitudinal trends of disease phenology -------------------
# for three different environmental change scenarios


# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("max", "mean")


for (o in operations) { # Loop over the max and mean functions
  
  print(o)
  
  # Read in the predictions data
  r_curr_preds_clim_landuse <- terra::rast("output_data/results/I_ricinus_preds_clim_landuse_1970_2019.tif") # under observed climate and land use change for the years 1970 to 2019
  r_curr_preds_clim <- terra::rast("output_data/results/I_ricinus_preds_clim_1970_2019.tif") # under observed climate change for the years 1970 to 2019 (mean land use values)
  
  r_fut_preds_clim_landuse <- terra::rast(paste0("output_data/results/I_ricinus_preds_clim_landuse_2030_2070_ssp370.tif")) # under scenario of future climate and land use change
  r_fut_preds_clim <- terra::rast(paste0("output_data/results/I_ricinus_preds_clim_2030_2070_ssp370.tif")) # under future climate change (mean land use values)
  
  # Generate past time information and assign dates as layer names
  dates_past <- seq(as.Date("1970-01-01"), as.Date("2019-12-01"), by = "month")
  names(r_curr_preds_clim_landuse) <- dates_past
  names(r_curr_preds_clim) <- dates_past
  
  # Generate future time information and assign dates as layer names
  dates_future <- seq(as.Date("2030-01-01"), as.Date("2070-12-01"), by = "month")
  names(r_fut_preds_clim_landuse) <- dates_future
  names(r_fut_preds_clim) <- dates_future
  
  # Convert raster stack to a data frame
  raster_df_curr_clim_landuse <- as.data.frame(r_curr_preds_clim_landuse, xy = TRUE)
  raster_df_curr_clim <- as.data.frame(r_curr_preds_clim, xy = TRUE)
  raster_df_fut_clim_landuse <- as.data.frame(r_fut_preds_clim_landuse, xy = TRUE)
  raster_df_fut_clim <- as.data.frame(r_fut_preds_clim, xy = TRUE)
  
  # Reshape data to a long format
  df_long_curr_clim_landuse <- raster_df_curr_clim_landuse %>%
    pivot_longer(
      cols = -c(x,y),
      names_to = "date",
      values_to = "probability"
    ) %>%
    mutate(year = as.numeric(substring(date, 1, 4)),
           month = substring(date, 6, 7),
           decade = paste0((year %/% 10) * 10, "s")) 
  
  df_long_curr_clim <- raster_df_curr_clim %>%
    pivot_longer(
      cols = -c(x,y),
      names_to = "date",
      values_to = "probability"
    ) %>%
    mutate(year = as.numeric(substring(date, 1, 4)),
           month = substring(date, 6, 7),
           decade = paste0((year %/% 10) * 10, "s")) 
  
  df_long_fut_clim_landuse <- raster_df_fut_clim_landuse %>%
    pivot_longer(
      cols = -c(x,y),
      names_to = "date",
      values_to = "probability"
    ) %>%
    mutate(year = as.numeric(substring(date, 1, 4)),
           month = substring(date, 6, 7),
           decade = paste0((year %/% 10) * 10, "s")) 
  
  df_long_fut_clim <- raster_df_fut_clim %>%
    pivot_longer(
      cols = -c(x,y),
      names_to = "date",
      values_to = "probability"
    ) %>%
    mutate(year = as.numeric(substring(date, 1, 4)),
           month = substring(date, 6, 7),
           decade = paste0((year %/% 10) * 10, "s")) 
  
  # Define 6 latitudinal bands within Europe
  lat_bands_curr_clim_landuse <- data.frame(lat_band = cut(df_long_curr_clim_landuse$y, breaks = seq(34, 72, length.out = 6), 
                                                           labels = paste0("Band ", 5:1), include.lowest = TRUE))
  
  lat_bands_curr_clim <- data.frame(lat_band = cut(df_long_curr_clim$y, breaks = seq(34, 72, length.out = 6), 
                                                   labels = paste0("Band ", 5:1), include.lowest = TRUE))
  
  lat_bands_fut_clim_landuse <- data.frame(lat_band = cut(df_long_fut_clim_landuse$y, breaks = seq(34, 72, length.out = 6), 
                                                          labels = paste0("Band ", 5:1), include.lowest = TRUE))
  
  lat_bands_fut_clim <- data.frame(lat_band = cut(df_long_fut_clim$y, breaks = seq(34, 72, length.out = 6), 
                                                  labels = paste0("Band ", 5:1), include.lowest = TRUE))
  
  # Add the latitudinal bands to the data frame
  df_long_curr_clim_landuse <- df_long_curr_clim_landuse %>%
    mutate(scenario = "Climate + Land use change",
           lat_band = lat_bands_curr_clim_landuse$lat_band)
  
  
  df_long_curr_clim <- df_long_curr_clim %>%
    mutate(scenario = "Climate change",
           lat_band = lat_bands_curr_clim$lat_band)
  
  df_long_fut_clim_landuse <- df_long_fut_clim_landuse %>%
    mutate(scenario = "Climate + Land use change",
           lat_band = lat_bands_fut_clim_landuse$lat)
  
  df_long_fut_clim <- df_long_fut_clim %>%
    mutate(scenario = "Climate change",
           lat_band = lat_bands_fut_clim$lat_band)
  
  # Bind the data frames 
  combined_df <- bind_rows(df_long_curr_clim_landuse, df_long_curr_clim,
                           df_long_fut_clim_landuse, df_long_fut_clim)
  
  
  # Summarise data by latitudinal band, month, and scenario
  df_summary_bands <- combined_df %>%
    group_by(lat_band, month, scenario, year) %>%
    summarise(occ_probability = match.fun(o)(probability, na.rm = TRUE), .groups = "drop")
  
  df_summary_bands <- df_summary_bands %>%
    mutate(
      decade = case_when(
        year >= 1970 & year < 1980 ~ "1970s",
        year >= 1980 & year < 1990 ~ "1980s",
        year >= 1990 & year < 2000 ~ "1990s",
        year >= 2000 & year < 2010 ~ "2000s",
        year >= 2010 & year < 2020 ~ "2010s",
        year >= 2030 & year < 2040 ~ "2030s",
        year >= 2040 & year < 2050 ~ "2040s",
        year >= 2050 & year < 2060 ~ "2050s",
        year >= 2060 & year < 2070 ~ "2060s",

      )
    )
  
  # keep summarising the data over the decade
  df_summary_bands <- df_summary_bands %>%
    group_by(lat_band, month, decade, scenario) %>%
    summarise(occ_probability = mean(occ_probability, na.rm = TRUE), .groups = "drop")
  
  # Convert month values to names
  df_summary_bands <- df_summary_bands %>%
    mutate(month = factor(month, levels = sprintf("%02d", 1:12), labels = month.abb))
  
  # Filter the data frame with summarised predictions from the 1970s, 2010s, and 2050s
  df_filtered <- df_summary_bands %>%
    filter(decade %in% c("1970s", "2010s", "2050s"))
  
  # Convert month values and latitudinal bands to factors for proper ordering in the plot
  df_filtered <- df_filtered %>%
    mutate(month = factor(month, levels = month.abb),
           lat_band = factor(lat_band, levels =  rev(paste0("Band ", 5:1))))
  
  # Create the plot
  ggplot(df_filtered, aes(x = month, y = occ_probability, color = lat_band, group = interaction(lat_band, scenario))) +
    geom_line(aes(linetype = scenario), linewidth = 1) +  
    facet_wrap(~decade, ncol = 1) +  
    scale_y_continuous(paste(o, "occurrence probability")) +  
    scale_x_discrete("month of a year") +  
    scale_color_viridis_d(
      option = "rocket",
      name = "latitudinal band"
    ) +
    theme_minimal() +  
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 12), 
      plot.title = element_text(face = "bold", size = 16)
    ) +
    labs(
      title = paste("Ixodes ricinus -", o, "occurrence probability\nby month, decade, and latitudinal band")
    )
  
  ggsave(paste0("output_data/plots/latitudinal_trends/I_ricinus_latitudinal_trends_",o,"occprob.png"), width = 9, height = 9)
                       
} # Close the loop over the mean and max functions





#-------------------------------------------------------------------------------

# 2. Plot the examined latitudinal bands over study plots ----------------------

# Read in the 50 km raster of Europe (used mask)
europe_mask <- terra::rast("input_data/spatial_data/europe_mask_50km.tif")

# Convert the raster to a data frame
europe_mask_df <- as.data.frame(europe_mask, xy = TRUE)

# Define breakpoints for continuous latitudinal bands
breakpoints <- seq(34, 72, length.out = 6)

# Create a data frame for latitudinal bands
lat_band_lines <- data.frame(
  lat_min = head(breakpoints, -1),  # All values except the last
  lat_max = tail(breakpoints, -1), # All values except the first
  lat_band = paste0("Band ", rev(seq_along(head(breakpoints, -1)))) # Reverse order for naming
)

# Convert 'lat_band' to a factor
lat_band_lines$lat_band <- factor(lat_band_lines$lat_band)

ggplot() +
  geom_tile(data = europe_mask_df, aes(x = x, y = y)) + 
  geom_rect(data = lat_band_lines, aes(xmin = -Inf, xmax = Inf, ymin = lat_min, ymax = lat_max, fill = lat_band), 
            alpha = 0.45) + 
  scale_fill_viridis_d(
    option = "rocket",
    name = "Latitudinal bands"
  ) +
  labs(title = "Latitudinal bands over Europe used for analysis") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  ) 

ggsave("output_data/plots/latitudinal_trends/latitudinal_bands_Europe.png", width = 7, height = 5)

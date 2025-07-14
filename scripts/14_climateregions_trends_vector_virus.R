# ZOE disease phenology analysis

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                   14. Trends of main climate regions                   #
# ---------------------------------------------------------------------- #

# Load needed packages
library(terra)
library(tidyverse)
library(ggplot2)
library(ggh4x)

# Load needed data
europe_mask_50km <- terra::rast("input_data/spatial_data/europe_mask_50km.tif")




#-------------------------------------------------------------------------------

# 1. Prepare raster of main climate regions in Europe --------------------------

# Source: 
# https://www.gloh2o.org/koppen/
# Beck, H.E., T.R. McVicar, N. Vergopolan, A. Berg, N.J. Lutsko, A. Dufour, Z. Zeng, 
# X. Jiang, A.I.J.M. van Dijk, D.G. Miralles High-resolution (1 km) Köppen-Geiger 
# maps for 1901–2099 based on constrained CMIP6 projectionsScientific Data 10, 724, 
# doi:10.1038/s41597-023–02549‑6 (2023)



# a) Historical climate regions (two time frames) ------------------------------

# Load the raster of main climate regions in Europe in a resolution of 50 km
# for two different time frames that cover our historical predictions
clim_world_1991_2020 <- terra::rast("input_data/spatial_data/climate_regions/koppen_geiger_0p5_1991_2020.tif")
clim_world_1961_1990 <- terra::rast("input_data/spatial_data/climate_regions/koppen_geiger_0p5_1961_1990.tif")

# Clip the map extent to Europe
clim_eur_1991_2020 <- terra::mask(crop(clim_world_1991_2020, europe_mask_50km), europe_mask_50km)
clim_eur_1961_1990 <- terra::mask(crop(clim_world_1961_1990, europe_mask_50km), europe_mask_50km)

# Check unique climate class values in Europe
unique_vals_1 <- unique(values(clim_eur_1991_2020))
unique_vals_2 <- unique(values(clim_eur_1961_1990))

# Define class labels based on the climate zone values
# 1991 - 2020
labels_1991_2020 <- data.frame(value = c(4, 6, 7, 8, 9, 14, 15, 16, 17, 18, 25, 26, 27, 29),
                     # main_climate = c("Dry", "Dry", "Dry", "Temperate", "Temperate",
                     #                  "Temperate", "Temperate", "Temperate", "Continental",
                     #                  "Continental", "Continental", "Continental", 
                     #                  "Continental", "Polar"),
                     # main_climate_sub = c("Dry; Hot desert", "Dry; Hot semi-arid", "Dry; Cold semi-arid",
                     #                      "Temperate; Hot-summer Mediterranean", "Temperate; Warm-summer Mediterranean",
                     #                      "Temperate; Humid subtropical", "Temperate; oceanic/subtropical highland",
                     #                      "Temperate; subpolar oceanic", "Continental; Hot-summer Mediterranean",
                     #                      "Continental; Warm-summer Mediterranean", "Continental; Hot-summer humid",
                     #                      "Continental: Warm-summer humid", "Continental; Subarctic", "Polar; Tundra"))
                     main_climate_sub_imp = c("Dry", "Dry", "Dry", "Temperate; Mediterranean", "Temperate; Mediterranean",
                                              "Temperate; Humid subtropical", "Temperate; Oceanic", "Temperate; Subpolar oceanic",
                                              "Continental; Mediterranean", "Continental; Mediterranean", "Continental; Humid",
                                              "Continental; Humid", "Continental; Subarctic", "Polar"))


# Assign labels 
levels(clim_eur_1991_2020) <- labels_1991_2020

# 1961 - 1990
labels_1961_1990 <- data.frame(value = c(5, 6, 7, 8, 9, 14, 15, 16, 17, 18, 25, 26, 27, 29),
                               # main_climate = c("Dry", "Dry", "Dry", "Temperate", "Temperate",
                               #                  "Temperate", "Temperate", "Temperate", "Continental",
                               #                  "Continental", "Continental", "Continental", 
                               #                  "Continental", "Polar"),
                               # main_climate_sub = c("Dry; Hot desert", "Dry; Hot semi-arid", "Dry; Cold semi-arid",
                               #                      "Temperate; Hot-summer Mediterranean", "Temperate; Warm-summer Mediterranean",
                               #                      "Temperate; Humid subtropical", "Temperate; oceanic/subtropical highland",
                               #                      "Temperate; subpolar oceanic", "Continental; Hot-summer Mediterranean",
                               #                      "Continental; Warm-summer Mediterranean", "Continental; Hot-summer humid",
                               #                      "Continental: Warm-summer humid", "Continental; Subarctic", "Polar; Tundra"))
                               main_climate_sub_imp = c("Dry", "Dry", "Dry", "Temperate; Mediterranean", "Temperate; Mediterranean",
                                                        "Temperate; Humid subtropical", "Temperate; Oceanic", "Temperate; Subpolar oceanic",
                                                        "Continental; Mediterranean", "Continental; Mediterranean", "Continental; Humid",
                                                        "Continental; Humid", "Continental; Subarctic", "Polar"))


# Assign labels 
levels(clim_eur_1961_1990) <- labels_1961_1990



# b) Future climate regions (one time frame, three ssp) ------------------------

# Load the raster of main climate regions in Europe in a resolution of 50 km
# for one time frame that covers our future predictions and for three different
# climate scenarios
clim_world_2041_2070_ssp126 <- terra::rast("input_data/spatial_data/climate_regions/koppen_geiger_0p5_2041_2070_ssp126.tif")
clim_world_2041_2070_ssp370 <- terra::rast("input_data/spatial_data/climate_regions/koppen_geiger_0p5_2041_2070_ssp370.tif")
clim_world_2041_2070_ssp585 <- terra::rast("input_data/spatial_data/climate_regions/koppen_geiger_0p5_2041_2070_ssp585.tif")

# Clip the map extent to Europe
clim_eur_2041_2070_ssp126 <- terra::mask(crop(clim_world_2041_2070_ssp126, europe_mask_50km), europe_mask_50km)
clim_eur_2041_2070_ssp370 <- terra::mask(crop(clim_world_2041_2070_ssp370, europe_mask_50km), europe_mask_50km)
clim_eur_2041_2070_ssp585 <- terra::mask(crop(clim_world_2041_2070_ssp585, europe_mask_50km), europe_mask_50km)

# Check unique climate class values in Europe
unique_vals_ssp126 <- unique(values(clim_eur_2041_2070_ssp126))
unique_vals_ssp370 <- unique(values(clim_eur_2041_2070_ssp370))
unique_vals_ssp585 <- unique(values(clim_eur_2041_2070_ssp585))

# Assign the same labels as historical data from 1991 to 2020 as the same
# climate zones are present
levels(clim_eur_2041_2070_ssp126) <- labels_1991_2020
levels(clim_eur_2041_2070_ssp370) <- labels_1991_2020
levels(clim_eur_2041_2070_ssp585) <- labels_1991_2020




#-------------------------------------------------------------------------------

# 2. Calculate past trends of main climatic regions ----------------------------
# for main vectors as well as viruses
# focusing on the decades 1970s and 2010s


# a) Load data -----------------------------------------------------------------

# Load the needed data - postprocessed monthly continuous predictions of the main 
# vectors and their associated viruses for the years 1970 to 2019 under factual 
# climate and land use change)
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_1970_2019.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")

# TBE
TBE_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_1970_2019.tif")

# WNV
WNV_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_1970_2019.tif")




# b) Calculate occurrence probabilities per climate region ---------------------

# Generate past time information and assign dates as layer names
dates_past <- seq(as.Date("1970-01-01"), as.Date("2019-12-01"), by = "month")
names(I_ricinus_occ_prob_clim_landuse_ens) <- dates_past
names(C_pipiens_occ_prob_clim_landuse_ens) <- dates_past
names(TBE_occ_prob_clim_ens) <- dates_past
names(WNV_occ_prob_clim_ens) <- dates_past

# Convert raster stack to a data frame
I_ricinus_occ_prob_clim_landuse_ens_df <- as.data.frame(I_ricinus_occ_prob_clim_landuse_ens, xy = TRUE)
C_pipiens_occ_prob_clim_landuse_ens_df <- as.data.frame(C_pipiens_occ_prob_clim_landuse_ens, xy = TRUE)
TBE_occ_prob_clim_ens_df <- as.data.frame(TBE_occ_prob_clim_ens, xy = TRUE)
WNV_occ_prob_clim_ens_df <- as.data.frame(WNV_occ_prob_clim_ens, xy = TRUE)

# Reshape data to a long format
I_ricinus_occ_prob_clim_landuse_ens_df_long <- I_ricinus_occ_prob_clim_landuse_ens_df %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

C_pipiens_occ_prob_clim_landuse_ens_df_long <- C_pipiens_occ_prob_clim_landuse_ens_df %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

TBE_occ_prob_clim_ens_df_long <- TBE_occ_prob_clim_ens_df %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

WNV_occ_prob_clim_ens_df_long <- WNV_occ_prob_clim_ens_df %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

# Filter the probability values by time to be able to assign climate classes
# according to the two time frames
# Ixodes ricinus
I_ricinus_1991_2019 <- I_ricinus_occ_prob_clim_landuse_ens_df_long %>%
  filter(year >= 1991 & year <= 2019)

I_ricinus_1970_1990 <- I_ricinus_occ_prob_clim_landuse_ens_df_long %>%
  filter(year >= 1970 & year <= 1990)

# Culex pipiens
C_pipiens_1991_2019 <- C_pipiens_occ_prob_clim_landuse_ens_df_long %>%
  filter(year >= 1991 & year <= 2019)

C_pipiens_1970_1990 <- C_pipiens_occ_prob_clim_landuse_ens_df_long %>%
  filter(year >= 1970 & year <= 1990)

# TBE
TBE_1991_2019 <- TBE_occ_prob_clim_ens_df_long %>%
  filter(year >= 1991 & year <= 2019)

TBE_1970_1990 <- TBE_occ_prob_clim_ens_df_long %>%
  filter(year >= 1970 & year <= 1990)

# WNV
WNV_1991_2019 <- WNV_occ_prob_clim_ens_df_long %>%
  filter(year >= 1991 & year <= 2019)

WNV_1970_1990 <- WNV_occ_prob_clim_ens_df_long %>%
  filter(year >= 1970 & year <= 1990)

# Convert the rasters containing the climate regions into a data frame
# with coordinates
df_clim_eur_1991_2020 <- as.data.frame(clim_eur_1991_2020, xy = TRUE, na.rm = TRUE)
colnames(df_clim_eur_1991_2020) <- c("x", "y", "clim_region")

df_clim_eur_1961_1990 <- as.data.frame(clim_eur_1961_1990, xy = TRUE, na.rm = TRUE)
colnames(df_clim_eur_1961_1990) <- c("x", "y", "clim_region")

# Join information on climate region with occurrence probabilities
# by x and y coordinate (for both time frames)
# Ixodes ricinus
I_ricinus_joined_1991_2019 <- I_ricinus_1991_2019 %>%
  inner_join(df_clim_eur_1991_2020, by = c("x", "y"))

I_ricinus_joined_1970_1990 <- I_ricinus_1970_1990 %>%
  inner_join(df_clim_eur_1961_1990, by = c("x", "y"))

I_ricinus_joined <- bind_rows(I_ricinus_joined_1991_2019, I_ricinus_joined_1970_1990)

# Culex pipiens
C_pipiens_joined_1991_2019 <- C_pipiens_1991_2019 %>%
  inner_join(df_clim_eur_1991_2020, by = c("x", "y"))

C_pipiens_joined_1970_1990 <- C_pipiens_1970_1990 %>%
  inner_join(df_clim_eur_1961_1990, by = c("x", "y"))

C_pipiens_joined <- bind_rows(C_pipiens_joined_1991_2019, C_pipiens_joined_1970_1990)

# TBE
TBE_joined_1991_2019 <- TBE_1991_2019 %>%
  inner_join(df_clim_eur_1991_2020, by = c("x", "y"))

TBE_joined_1970_1990 <- TBE_1970_1990 %>%
  inner_join(df_clim_eur_1961_1990, by = c("x", "y"))

TBE_joined <- bind_rows(TBE_joined_1991_2019, TBE_joined_1970_1990)

# WNV
WNV_joined_1991_2019 <- WNV_1991_2019 %>%
  inner_join(df_clim_eur_1991_2020, by = c("x", "y"))

WNV_joined_1970_1990 <- WNV_1970_1990 %>%
  inner_join(df_clim_eur_1961_1990, by = c("x", "y"))

WNV_joined <- bind_rows(WNV_joined_1991_2019, WNV_joined_1970_1990)



# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")


for (o in operations) { # Loop over the peak and mean functions
  
  print(o)
  
  # Calculate the mean habitat suitability per biogeographic region
  # for each month per year
  I_ricinus_summary <- I_ricinus_joined %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  C_pipiens_summary <- C_pipiens_joined %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  TBE_summary <- TBE_joined %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  WNV_summary <- WNV_joined %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  
  # Add the name of the vector species and virus to the data frame
  I_ricinus_summary <- I_ricinus_summary %>%
    mutate(species = "Ixodes ricinus",
           scenario = "Observed data")
  
  C_pipiens_summary <- C_pipiens_summary %>%
    mutate(species = "Culex pipiens",
           scenario = "Observed data")
  
  TBE_summary <- TBE_summary %>%
    mutate(species = "TBE",
           scenario = "Observed data")
  
  WNV_summary <- WNV_summary %>%
    mutate(species = "WNV",
           scenario = "Observed data")
  
  # Bind the data frames 
  combined_df <- bind_rows(I_ricinus_summary, C_pipiens_summary,
                           TBE_summary, WNV_summary)
  
  # Add the decade name as column
  combined_df <- combined_df %>%
    mutate(
      decade = case_when(
        year >= 1970 & year < 1980 ~ "1970s",
        year >= 1980 & year < 1990 ~ "1980s",
        year >= 1990 & year < 2000 ~ "1990s",
        year >= 2000 & year < 2010 ~ "2000s",
        year >= 2010 & year < 2020 ~ "2010s"))
  
  # Keep summarising the data over the decade, month, and vector/virus per biogeographic region
  aggregated_df_past <- combined_df %>%
    group_by(clim_region, decade, month, species, scenario) %>%
    summarise(occ_probability = if (o == "Peak") quantile(mean_probability, 0.95, na.rm = TRUE)
              else if (o == "Mean") mean(mean_probability, na.rm = TRUE),
              .groups = "drop")
  
  
  
  
# c) Prepare data frame for plotting -------------------------------------------
  
  # Convert month values to names
  aggregated_df_past  <- aggregated_df_past %>%
    mutate(month = factor(month, levels = sprintf("%02d", 1:12), labels = month.abb))
  
  # Filter the data frame with summarised predictions from the 1970s, and 2010s
  aggregated_df_past <- aggregated_df_past %>%
    filter(decade %in% c("1970s", "2010s"))
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_past <- aggregated_df_past %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Convert month values and latitudinal bands to factors for proper ordering in the plot
  aggregated_df_past <- aggregated_df_past %>%
    mutate(month = factor(month, levels = month.abb))
  
  
  # Save the resulting data frame
  save(aggregated_df_past, file = paste0("output_data/results/climateregions_trends/climateregions_trends_past_vector_virus_",o,".RData"))
  
} # Close the loop over the peak and mean functions



#-------------------------------------------------------------------------------

# 3. Calculate future trends of main climatic regions  -------------------------
# for main vectors as well as viruses
# focusing on the decade 2050s


# a) Load data -----------------------------------------------------------------

# Load needed data - postprocessed monthly predicted ensemble occurrence probabilities 
# of main vectors and the respective viruses, under climate and land use change for the
# future years 2030 to 2070 based and three studied environmental scenarios
# Ixodes ricinus; ssp126, ssp370, and ssp585
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp126.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp585.tif"))

# Culex pipiens; ssp126, ssp370, and ssp585
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp585.tif"))

# TBE; ssp126, ssp370, and ssp585
TBE_occ_prob_clim_ens_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp126.tif"))
TBE_occ_prob_clim_ens_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp370.tif"))
TBE_occ_prob_clim_ens_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp585.tif"))

# WNV; ssp126, ssp370, and ssp585
WNV_occ_prob_clim_ens_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp126.tif"))
WNV_occ_prob_clim_ens_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp370.tif"))
WNV_occ_prob_clim_ens_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp585.tif"))



# b) Calculate occurrence probabilities per climate region ---------------------

# Generate future time information and assign dates as layer names
dates_future <- seq(as.Date("2030-01-01"), as.Date("2070-12-01"), by = "month")
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126) <- dates_future
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370) <- dates_future
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585) <- dates_future
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126) <- dates_future
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370) <- dates_future
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585) <- dates_future
names(TBE_occ_prob_clim_ens_fut_ssp126) <- dates_future
names(TBE_occ_prob_clim_ens_fut_ssp370) <- dates_future
names(TBE_occ_prob_clim_ens_fut_ssp585) <- dates_future
names(WNV_occ_prob_clim_ens_fut_ssp126) <- dates_future
names(WNV_occ_prob_clim_ens_fut_ssp370) <- dates_future
names(WNV_occ_prob_clim_ens_fut_ssp585) <- dates_future


# Convert raster stack to a data frame
I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp126 <- as.data.frame(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, xy = TRUE)
I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp370 <- as.data.frame(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, xy = TRUE)
I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp585 <- as.data.frame(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, xy = TRUE)

C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp126 <- as.data.frame(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, xy = TRUE)
C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp370 <- as.data.frame(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, xy = TRUE)
C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp585 <- as.data.frame(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, xy = TRUE)

TBE_occ_prob_clim_ens_fut_df_ssp126 <- as.data.frame(TBE_occ_prob_clim_ens_fut_ssp126, xy = TRUE)
TBE_occ_prob_clim_ens_fut_df_ssp370 <- as.data.frame(TBE_occ_prob_clim_ens_fut_ssp370, xy = TRUE)
TBE_occ_prob_clim_ens_fut_df_ssp585 <- as.data.frame(TBE_occ_prob_clim_ens_fut_ssp585, xy = TRUE)

WNV_occ_prob_clim_ens_fut_df_ssp126 <- as.data.frame(WNV_occ_prob_clim_ens_fut_ssp126, xy = TRUE)
WNV_occ_prob_clim_ens_fut_df_ssp370 <- as.data.frame(WNV_occ_prob_clim_ens_fut_ssp370, xy = TRUE)
WNV_occ_prob_clim_ens_fut_df_ssp585 <- as.data.frame(WNV_occ_prob_clim_ens_fut_ssp585, xy = TRUE)


# Reshape data to a long format
I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp126 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp126 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp370 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp370 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp585 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_ssp585 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp126 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp126 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp370 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp370 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp585 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_ssp585 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

TBE_occ_prob_clim_ens_fut_df_long_ssp126 <- TBE_occ_prob_clim_ens_fut_df_ssp126 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

TBE_occ_prob_clim_ens_fut_df_long_ssp370 <- TBE_occ_prob_clim_ens_fut_df_ssp370 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

TBE_occ_prob_clim_ens_fut_df_long_ssp585 <- TBE_occ_prob_clim_ens_fut_df_ssp585 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

WNV_occ_prob_clim_ens_fut_df_long_ssp126 <- WNV_occ_prob_clim_ens_fut_df_ssp126 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

WNV_occ_prob_clim_ens_fut_df_long_ssp370 <- WNV_occ_prob_clim_ens_fut_df_ssp370 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 

WNV_occ_prob_clim_ens_fut_df_long_ssp585 <- WNV_occ_prob_clim_ens_fut_df_ssp585 %>%
  pivot_longer(
    cols = -c(x,y),
    names_to = "date",
    values_to = "probability"
  ) %>%
  mutate(year = as.numeric(substring(date, 1, 4)),
         month = substring(date, 6, 7),
         decade = paste0((year %/% 10) * 10, "s")) 


# Convert the rasters containing the climatic regions into a data frame
# with coordinates (based on each environmental ssp scenario)
df_clim_eur_2041_2070_ssp126 <- as.data.frame(clim_eur_2041_2070_ssp126, xy = TRUE, na.rm = TRUE)
colnames(df_clim_eur_2041_2070_ssp126) <- c("x", "y", "clim_region")

df_clim_eur_2041_2070_ssp370 <- as.data.frame(clim_eur_2041_2070_ssp370, xy = TRUE, na.rm = TRUE)
colnames(df_clim_eur_2041_2070_ssp370) <- c("x", "y", "clim_region")

df_clim_eur_2041_2070_ssp585 <- as.data.frame(clim_eur_2041_2070_ssp585, xy = TRUE, na.rm = TRUE)
colnames(df_clim_eur_2041_2070_ssp585) <- c("x", "y", "clim_region")

# Join information on biogeographic region with occurrence probabilities
# by x and y coordinate
I_ricinus_joined_ssp126 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp126 %>%
  inner_join(df_clim_eur_2041_2070_ssp126, by = c("x", "y"))

I_ricinus_joined_ssp370 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp370 %>%
  inner_join(df_clim_eur_2041_2070_ssp370, by = c("x", "y"))

I_ricinus_joined_ssp585 <- I_ricinus_occ_prob_clim_landuse_ens_fut_df_long_ssp585 %>%
  inner_join(df_clim_eur_2041_2070_ssp585, by = c("x", "y"))

C_pipiens_joined_ssp126 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp126 %>%
  inner_join(df_clim_eur_2041_2070_ssp126, by = c("x", "y"))

C_pipiens_joined_ssp370 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp370 %>%
  inner_join(df_clim_eur_2041_2070_ssp370, by = c("x", "y"))

C_pipiens_joined_ssp585 <- C_pipiens_occ_prob_clim_landuse_ens_fut_df_long_ssp585 %>%
  inner_join(df_clim_eur_2041_2070_ssp585, by = c("x", "y"))

TBE_joined_ssp126 <- TBE_occ_prob_clim_ens_fut_df_long_ssp126 %>%
  inner_join(df_clim_eur_2041_2070_ssp126, by = c("x", "y"))

TBE_joined_ssp370 <- TBE_occ_prob_clim_ens_fut_df_long_ssp370 %>%
  inner_join(df_clim_eur_2041_2070_ssp370, by = c("x", "y"))

TBE_joined_ssp585 <- TBE_occ_prob_clim_ens_fut_df_long_ssp585 %>%
  inner_join(df_clim_eur_2041_2070_ssp585, by = c("x", "y"))

WNV_joined_ssp126 <- WNV_occ_prob_clim_ens_fut_df_long_ssp126 %>%
  inner_join(df_clim_eur_2041_2070_ssp126, by = c("x", "y"))

WNV_joined_ssp370 <- WNV_occ_prob_clim_ens_fut_df_long_ssp370 %>%
  inner_join(df_clim_eur_2041_2070_ssp370, by = c("x", "y"))

WNV_joined_ssp585 <- WNV_occ_prob_clim_ens_fut_df_long_ssp585 %>%
  inner_join(df_clim_eur_2041_2070_ssp585, by = c("x", "y"))

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")


for (o in operations) { # Loop over the peak and mean functions
  
  print(o)
  
  # Calculate the mean habitat suitability per biogeographic region
  # for each month per year
  I_ricinus_summary_ssp126 <- I_ricinus_joined_ssp126 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  I_ricinus_summary_ssp370 <- I_ricinus_joined_ssp370 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  I_ricinus_summary_ssp585 <- I_ricinus_joined_ssp585 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  C_pipiens_summary_ssp126 <- C_pipiens_joined_ssp126 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  C_pipiens_summary_ssp370 <- C_pipiens_joined_ssp370 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  C_pipiens_summary_ssp585 <- C_pipiens_joined_ssp585 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  TBE_summary_ssp126 <- TBE_joined_ssp126 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  TBE_summary_ssp370 <- TBE_joined_ssp370 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  TBE_summary_ssp585 <- TBE_joined_ssp585 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  WNV_summary_ssp126 <- WNV_joined_ssp126 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  WNV_summary_ssp370 <- WNV_joined_ssp370 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  WNV_summary_ssp585 <- WNV_joined_ssp585 %>%
    group_by(clim_region, year, month) %>%
    summarise(mean_probability = mean(probability, na.rm = TRUE), .groups = "drop")
  
  # Add the name of the vector species and virus to the data frame
  I_ricinus_summary_ssp126 <- I_ricinus_summary_ssp126 %>%
    mutate(species = "Ixodes ricinus",
           scenario = "Projected future ssp126")
  
  I_ricinus_summary_ssp370 <- I_ricinus_summary_ssp370 %>%
    mutate(species = "Ixodes ricinus",
           scenario = "Projected future ssp370")
  
  I_ricinus_summary_ssp585 <- I_ricinus_summary_ssp585 %>%
    mutate(species = "Ixodes ricinus",
           scenario = "Projected future ssp585")
  
  C_pipiens_summary_ssp126 <- C_pipiens_summary_ssp126 %>%
    mutate(species = "Culex pipiens",
           scenario = "Projected future ssp126")
  
  C_pipiens_summary_ssp370 <- C_pipiens_summary_ssp370 %>%
    mutate(species = "Culex pipiens",
           scenario = "Projected future ssp370")
  
  C_pipiens_summary_ssp585 <- C_pipiens_summary_ssp585 %>%
    mutate(species = "Culex pipiens",
           scenario = "Projected future ssp585")
  
  TBE_summary_ssp126 <- TBE_summary_ssp126 %>%
    mutate(species = "TBE",
           scenario = "Projected future ssp126")
  
  TBE_summary_ssp370 <- TBE_summary_ssp370 %>%
    mutate(species = "TBE",
           scenario = "Projected future ssp370")
  
  TBE_summary_ssp585 <- TBE_summary_ssp585 %>%
    mutate(species = "TBE",
           scenario = "Projected future ssp585")
  
  WNV_summary_ssp126 <- WNV_summary_ssp126 %>%
    mutate(species = "WNV",
           scenario = "Projected future ssp126")
  
  WNV_summary_ssp370 <- WNV_summary_ssp370 %>%
    mutate(species = "WNV",
           scenario = "Projected future ssp370")
  
  WNV_summary_ssp585 <- WNV_summary_ssp585 %>%
    mutate(species = "WNV",
           scenario = "Projected future ssp585")
  
  # Bind the data frames 
  combined_df <- bind_rows(I_ricinus_summary_ssp126, C_pipiens_summary_ssp126, TBE_summary_ssp126, WNV_summary_ssp126,
                           I_ricinus_summary_ssp370, C_pipiens_summary_ssp370, TBE_summary_ssp370, WNV_summary_ssp370,
                           I_ricinus_summary_ssp585, C_pipiens_summary_ssp585, TBE_summary_ssp585, WNV_summary_ssp585)
  
  # Add the decade name as column
  combined_df <- combined_df %>%
    mutate(
      decade = case_when(
        year >= 2030 & year < 2040 ~ "2030s",
        year >= 2040 & year < 2050 ~ "2040s",
        year >= 2050 & year < 2060 ~ "2050s",
        year >= 2060 & year < 2070 ~ "2060s"))
  
  # Keep summarising the data over the decade, month, and vector/virus per biogeographic region
  aggregated_df_fut <- combined_df %>%
    group_by(clim_region, decade, month, species, scenario) %>%
    summarise(occ_probability = if (o == "Peak") quantile(mean_probability, 0.95, na.rm = TRUE)
              else if (o == "Mean") mean(mean_probability, na.rm = TRUE),
              .groups = "drop")
  
  
  
  
# c) Prepare data frame for plotting -------------------------------------------
  
  # Convert month values to names
  aggregated_df_fut  <- aggregated_df_fut %>%
    mutate(month = factor(month, levels = sprintf("%02d", 1:12), labels = month.abb))
  
  # Filter the data frame with summarised predictions from the 1970s, and 2010s
  aggregated_df_fut <- aggregated_df_fut %>%
    filter(decade %in% c("2050s"))
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_fut <- aggregated_df_fut %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Convert month values and latitudinal bands to factors for proper ordering in the plot
  aggregated_df_fut <- aggregated_df_fut %>%
    mutate(month = factor(month, levels = month.abb))
  
  
  
  # Save the resulting data frame
  save(aggregated_df_fut, file = paste0("output_data/results/climateregions_trends/climateregions_trends_fut_vector_virus_",o,".RData"))
  
} # Close the loop over the peak and mean functions







#-------------------------------------------------------------------------------

# 3. Visualise the main climate regions in Europe ------------------------------
# for three different time intervals

# a) Prepare data for plotting -------------------------------------------------

# Load needed data
# The masks containing the countries with a mandatory reporting system of 
# disease surveillance data to the ECDC; we load different masks depending on
# the virus as Austria did not report of NUTS3 level for TBE
eu_eea_mask_WNV <- terra::rast("input_data/spatial_data/eu_eea_mask_WNV.tif") 
eu_eea_mask_TBE <- terra::rast("input_data/spatial_data/eu_eea_mask_TBE.tif")

# Align the extent of Europe mask the extent of WNV/TBE mask 
# (containing the countries with mandatory reporting; it does not matter if we take
# TBE or WNV mask as both have the same extent)
# 1961 - 1990
clim_eur_1961_1990 <- crop(clim_eur_1961_1990, ext(eu_eea_mask_WNV))
clim_eur_1961_1990 <- mask(clim_eur_1961_1990, eu_eea_mask_WNV)

# 1991 - 2020
clim_eur_1991_2020 <- crop(clim_eur_1991_2020, ext(eu_eea_mask_WNV))
clim_eur_1991_2020 <- mask(clim_eur_1991_2020, eu_eea_mask_WNV)

# 2041 - 2070 (environmental scenario ssp126)
clim_eur_2041_2070_ssp126 <- crop(clim_eur_2041_2070_ssp126, ext(eu_eea_mask_WNV))
clim_eur_2041_2070_ssp126 <- mask(clim_eur_2041_2070_ssp126, eu_eea_mask_WNV)

# 2041 - 2070 (environmental scenario ssp370)
clim_eur_2041_2070_ssp370 <- crop(clim_eur_2041_2070_ssp370, ext(eu_eea_mask_WNV))
clim_eur_2041_2070_ssp370 <- mask(clim_eur_2041_2070_ssp370, eu_eea_mask_WNV)

# 2041 - 2070 (environmental scenario ssp585)
clim_eur_2041_2070_ssp585 <- crop(clim_eur_2041_2070_ssp585, ext(eu_eea_mask_WNV))
clim_eur_2041_2070_ssp585 <- mask(clim_eur_2041_2070_ssp585, eu_eea_mask_WNV)


# Create a data frame from rasters
# 1961 - 1990
df_clim_eur_1961_1990_masked <- as.data.frame(clim_eur_1961_1990, xy = TRUE)
colnames(df_clim_eur_1961_1990_masked) <- c("x", "y", "clim_region")

# 1991 - 2020
df_clim_eur_1991_2020_masked <- as.data.frame(clim_eur_1991_2020, xy = TRUE)
colnames(df_clim_eur_1991_2020_masked) <- c("x", "y", "clim_region")

# 2041 - 2070 (environmental scenario ssp126)
df_clim_eur_2041_2070_ssp126_masked <- as.data.frame(clim_eur_2041_2070_ssp126, xy = TRUE)
colnames(df_clim_eur_2041_2070_ssp126_masked) <- c("x", "y", "clim_region")

# 2041 - 2070 (environmental scenario ssp370)
df_clim_eur_2041_2070_ssp370_masked <- as.data.frame(clim_eur_2041_2070_ssp370, xy = TRUE)
colnames(df_clim_eur_2041_2070_ssp370_masked) <- c("x", "y", "clim_region")

# 2041 - 2070 (environmental scenario ssp585)
df_clim_eur_2041_2070_ssp585_masked <- as.data.frame(clim_eur_2041_2070_ssp585, xy = TRUE)
colnames(df_clim_eur_2041_2070_ssp585_masked) <- c("x", "y", "clim_region")

# Turn climatic region into factor
df_clim_eur_1961_1990_masked$clim_region <- factor(df_clim_eur_1961_1990_masked$clim_region)
df_clim_eur_1991_2020_masked$clim_region <- factor(df_clim_eur_1991_2020_masked$clim_region)
df_clim_eur_2041_2070_ssp126_masked$clim_region <- factor(df_clim_eur_2041_2070_ssp126_masked$clim_region)
df_clim_eur_2041_2070_ssp370_masked$clim_region <- factor(df_clim_eur_2041_2070_ssp370_masked$clim_region)
df_clim_eur_2041_2070_ssp585_masked$clim_region <- factor(df_clim_eur_2041_2070_ssp585_masked$clim_region)

# Add a column to the data frame indicating the time frame
df_clim_eur_1961_1990_masked$time_frame <- "1970s"
df_clim_eur_1991_2020_masked$time_frame <- "2010s"
df_clim_eur_2041_2070_ssp126_masked$time_frame <- "2050s; ssp126"
df_clim_eur_2041_2070_ssp370_masked$time_frame <- "2050s; ssp370"
df_clim_eur_2041_2070_ssp585_masked$time_frame <- "2050s; ssp585"

# Add a column indicating the purpose of the plot 
df_clim_eur_1961_1990_masked$purpose <- "Köppen-Geiger climates"
df_clim_eur_1991_2020_masked$purpose <- "Köppen-Geiger climates"
df_clim_eur_2041_2070_ssp126_masked$purpose <- "Köppen-Geiger climates"
df_clim_eur_2041_2070_ssp370_masked$purpose <- "Köppen-Geiger climates"
df_clim_eur_2041_2070_ssp585_masked$purpose <- "Köppen-Geiger climates"

# Combine all into one data frame
df_all_clim_ssp126 <- bind_rows(df_clim_eur_1961_1990_masked,
                                df_clim_eur_1991_2020_masked,
                                df_clim_eur_2041_2070_ssp126_masked)

df_all_clim_ssp370 <- bind_rows(df_clim_eur_1961_1990_masked,
                                df_clim_eur_1991_2020_masked,
                                df_clim_eur_2041_2070_ssp370_masked)

df_all_clim_ssp585 <- bind_rows(df_clim_eur_1961_1990_masked,
                                df_clim_eur_1991_2020_masked,
                                df_clim_eur_2041_2070_ssp585_masked)


# Convert Europe mask spatraster into a data frame
europe_mask_50km_df <- as.data.frame(europe_mask_50km, xy = TRUE)


# b) Check frequency of climate regions -----------------------------------------

# Extract the frequency of cell per climate region for each time frame
# to set cells to NA with very small climate cells (not considered in analysis)
# 1961 - 1990
freq(clim_eur_1961_1990)

# 1991 - 2020
freq(clim_eur_1991_2020)

# 2041 - 2070 (environmental scenario ssp126)
freq(clim_eur_2041_2070_ssp126)

# 2041 - 2070 (environmental scenario ssp370)
freq(clim_eur_2041_2070_ssp370)

# 2041 - 2070 (environmental scenario ssp585)
freq(clim_eur_2041_2070_ssp585)

# Remove rows/cells of climate regions belonging to regions with low amount of
# cells
df_all_clim_ssp126[df_all_clim_ssp126$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA
df_all_clim_ssp370[df_all_clim_ssp370$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA
df_all_clim_ssp585[df_all_clim_ssp585$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA

# Remove NA rows
df_all_clim_ssp126 <- na.omit(df_all_clim_ssp126)
df_all_clim_ssp370 <- na.omit(df_all_clim_ssp370)
df_all_clim_ssp585 <- na.omit(df_all_clim_ssp585)





# c) Visualise main climate regions in Europe  ---------------------------------
# for the two past time frames and the future time frame for three different
# ssp environmental scenarios


# Visualise maps with main climate regions (ssp126)
ggplot(df_all_clim_ssp126, aes(x = x, y = y, fill = clim_region)) +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster() +
  facet_grid2(rows = vars(purpose), cols = vars(time_frame),
    strip = strip_themed(
      background_x = element_rect(fill = "grey75", color = NA),
      background_y = element_rect(fill = "grey75", color = NA))) +
  scale_fill_manual(
    values = c(
      "Dry" = "burlywood2",
      "Temperate; Mediterranean" = "darkorchid4",
      "Temperate; Humid subtropical" = "mediumorchid",
      "Temperate; Oceanic" = "thistle3",
      "Continental; Humid" = "#4D84C4",
      "Continental; Subarctic" = "skyblue1",
      "Polar" = "grey55"
    ), name = "Climate classifications")  +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "grey75", color = NA)) +
  guides(
    fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE))

# Save the plot
ggsave("output_data/plots/climateregions_trends/mainclimate_regions_Europe_ssp126.png", width = 11.5, height = 5)


# Visualise maps with main climate regions (ssp370)
ggplot(df_all_clim_ssp370, aes(x = x, y = y, fill = clim_region)) +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster() +
  facet_grid2(rows = vars(purpose), cols = vars(time_frame),
              strip = strip_themed(
                background_x = element_rect(fill = "grey75", color = NA),
                background_y = element_rect(fill = "grey75", color = NA))) +
  scale_fill_manual(
    values = c(
      "Dry" = "burlywood2",
      "Temperate; Mediterranean" = "darkorchid4",
      "Temperate; Humid subtropical" = "mediumorchid",
      "Temperate; Oceanic" = "thistle3",
      "Continental; Humid" = "#4D84C4",
      "Continental; Subarctic" = "skyblue1",
      "Polar" = "grey55"
    ), name = "Climate classifications")  +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "grey75", color = NA)) +
  guides(
    fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE))

# Save the plot
ggsave("output_data/plots/climateregions_trends/mainclimate_regions_Europe_ssp370.png", width = 11.5, height = 5)


# Visualise maps with main climate regions (ssp585)
ggplot(df_all_clim_ssp585, aes(x = x, y = y, fill = clim_region)) +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster() +
  facet_grid2(rows = vars(purpose), cols = vars(time_frame),
              strip = strip_themed(
                background_x = element_rect(fill = "grey75", color = NA),
                background_y = element_rect(fill = "grey75", color = NA))) +
  scale_fill_manual(
    values = c(
      "Dry" = "burlywood2",
      "Temperate; Mediterranean" = "darkorchid4",
      "Temperate; Humid subtropical" = "mediumorchid",
      "Temperate; Oceanic" = "thistle3",
      "Continental; Humid" = "#4D84C4",
      "Continental; Subarctic" = "skyblue1",
      "Polar" = "grey55"
    ), name = "Climate classifications")  +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 12.5, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_rect(fill = "grey75", color = NA)) +
  guides(
    fill = guide_legend(title.position = "top", nrow = 2, byrow = TRUE))

# Save the plot
ggsave("output_data/plots/climateregions_trends/mainclimate_regions_Europe_ssp585.png", width = 11.5, height = 5)










#-------------------------------------------------------------------------------

# 4. Visualise trends for main climatic regions --------------------------------
# for main vectors as well as viruses


# a) Prepare data frames for plotting ------------------------------------------
  
# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")

for (o in operations) { # Loop over the peak and mean functions
  
  print(o)
  
  # Load the needed data of past and future decadal trends of main climate regions
  load(paste0("output_data/results/climateregions_trends/climateregions_trends_past_vector_virus_",o,".RData"))
  load(paste0("output_data/results/climateregions_trends/climateregions_trends_fut_vector_virus_",o,".RData"))
  
  # Subset data for single environmental scenarios
  aggregated_df_futssp126 <- aggregated_df_fut[aggregated_df_fut$scenario == "Projected future ssp126", ]
  aggregated_df_futssp370 <- aggregated_df_fut[aggregated_df_fut$scenario == "Projected future ssp370", ]
  aggregated_df_futssp585 <- aggregated_df_fut[aggregated_df_fut$scenario == "Projected future ssp585", ]
  
  
  # Bind the two data frames containing the information of past and future decadal
  # trends of main clmate regions for Ixodes ricinus, TBE, Culex pipiens, and WNV
  decadal_trends_past_futssp126 <- rbind(aggregated_df_past, aggregated_df_futssp126)
  decadal_trends_past_futssp370 <- rbind(aggregated_df_past, aggregated_df_futssp370)
  decadal_trends_past_futssp585 <- rbind(aggregated_df_past, aggregated_df_futssp585)
  
  # Make sure the months are correctly ordered from January to December
  decadal_trends_past_futssp126$month <- factor(decadal_trends_past_futssp126$month, 
                                                levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  decadal_trends_past_futssp370$month <- factor(decadal_trends_past_futssp370$month, 
                                                levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  decadal_trends_past_futssp585$month <- factor(decadal_trends_past_futssp585$month, 
                                                levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                           "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vectors and diseases appear in the correct order
  decadal_trends_past_futssp126$species <- factor(decadal_trends_past_futssp126$species, 
                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  decadal_trends_past_futssp370$species <- factor(decadal_trends_past_futssp370$species, 
                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  decadal_trends_past_futssp585$species <- factor(decadal_trends_past_futssp585$species, 
                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  # Make sure that the future decade of the 2050s includes the corresponding environmental scenario name
  decadal_trends_past_futssp126$decade <- as.character(decadal_trends_past_futssp126$decade)
  decadal_trends_past_futssp126$decade[decadal_trends_past_futssp126$decade == "2050s"] <- "2050s; ssp126"
  
  decadal_trends_past_futssp370$decade <- as.character(decadal_trends_past_futssp370$decade)
  decadal_trends_past_futssp370$decade[decadal_trends_past_futssp370$decade == "2050s"] <- "2050s; ssp370"
  
  decadal_trends_past_futssp585$decade <- as.character(decadal_trends_past_futssp585$decade)
  decadal_trends_past_futssp585$decade[decadal_trends_past_futssp585$decade == "2050s"] <- "2050s; ssp585"
  
  
  # Remove rows that belong to climate regions that have a very low number of cells
  # (these will not be considered in the analysis)
  decadal_trends_past_futssp126[decadal_trends_past_futssp126$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA
  decadal_trends_past_futssp370[decadal_trends_past_futssp370$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA
  decadal_trends_past_futssp585[decadal_trends_past_futssp585$clim_region %in% c("Temperate; Subpolar oceanic", "Continental; Mediterranean"), ] <- NA
  
  # Remove rows containing NA values
  decadal_trends_past_futssp126 <- na.omit(decadal_trends_past_futssp126)
  decadal_trends_past_futssp370 <- na.omit(decadal_trends_past_futssp370)
  decadal_trends_past_futssp585 <- na.omit(decadal_trends_past_futssp585)
  
  
# b) Visualise phenology for main climate regions (past + future ssp126) -------  
  
  # Create the plot for observed environmental data and environmental scenario ssp126
  ggplot(decadal_trends_past_futssp126, aes(x = month, y = occ_probability, color = clim_region, linetype = scenario, group = interaction(clim_region, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +  
    facet_grid2(rows = vars(species), cols = vars(decade), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3"),            
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Latitudinal band", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "Dry" = "burlywood2",
        "Temperate; Mediterranean" = "darkorchid4",
        "Temperate; Humid subtropical" = "mediumorchid",
        "Temperate; Oceanic" = "thistle3",
        "Continental; Humid" = "#4D84C4",
        "Continental; Subarctic" = "skyblue1",
        "Polar" = "grey55"
      ), name = "Climate classifications" 
    ) +
    scale_linetype_manual(
      values = c(
        "Observed data" = "solid",
        "Projected future ssp126" = "solid"
      )
    ) +
    theme_bw() +  
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 14),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12.5, face = "bold"),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA)) +
    guides(
      color = guide_legend(title.position = "top", nrow = 2, byrow = TRUE),  
      linetype = "none")
  
  ggsave(paste0("output_data/plots/climateregions_trends/climateregions_trends_vector_virus_ssp126_",o,".png"), width = 11.5, height = 11)
  
  
  
# c) Visualise phenology for main climate regions (past + future ssp370) -------
  
  # Create the plot for observed environmental data and environmental scenario ssp370
  ggplot(decadal_trends_past_futssp370, aes(x = month, y = occ_probability, color = clim_region, linetype = scenario, group = interaction(clim_region, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +  
    facet_grid2(rows = vars(species), cols = vars(decade), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3"),            
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Latitudinal band", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "Dry" = "burlywood2",
        "Temperate; Mediterranean" = "darkorchid4",
        "Temperate; Humid subtropical" = "mediumorchid",
        "Temperate; Oceanic" = "thistle3",
        "Continental; Humid" = "#4D84C4",
        "Continental; Subarctic" = "skyblue1",
        "Polar" = "grey55"
      ), name = "Climate classifications" 
      ) +
    scale_linetype_manual(
      values = c(
        "Observed data" = "solid",
        "Projected future ssp370" = "solid"
      )
    ) +
    theme_bw() +  
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 14),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12.5, face = "bold"),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA)) +
    guides(
      color = guide_legend(title.position = "top", nrow = 2, byrow = TRUE),  
      linetype = "none")
  
  ggsave(paste0("output_data/plots/climateregions_trends/climateregions_trends_vector_virus_ssp370_",o,".png"), width = 11.5, height = 11)
  
  
  
# d) Visualise phenology for main climate regions (past + future ssp585) -------
  
  # Create the plot for observed environmental data and environmental scenario ssp370
  ggplot(decadal_trends_past_futssp585, aes(x = month, y = occ_probability, color = clim_region, linetype = scenario, group = interaction(clim_region, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +  
    facet_grid2(rows = vars(species), cols = vars(decade), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3"),            
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Latitudinal band", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "Dry" = "burlywood2",
        "Temperate; Mediterranean" = "darkorchid4",
        "Temperate; Humid subtropical" = "mediumorchid",
        "Temperate; Oceanic" = "thistle3",
        "Continental; Humid" = "#4D84C4",
        "Continental; Subarctic" = "skyblue1",
        "Polar" = "grey55"
      ), name = "Climate classifications" 
    ) +
    scale_linetype_manual(
      values = c(
        "Observed data" = "solid",
        "Projected future ssp585" = "solid"
      )
    ) +
    theme_bw() +  
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 14),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12.5, face = "bold"),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 16, face = "bold"),
      strip.text = element_text(size = 16, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA)) +
    guides(
      color = guide_legend(title.position = "top", nrow = 2, byrow = TRUE),  
      linetype = "none")
  
  ggsave(paste0("output_data/plots/climateregions_trends/climateregions_trends_vector_virus_ssp585_",o,".png"), width = 11.5, height = 11)
  
  
} # Close the loop over the peak and mean functions







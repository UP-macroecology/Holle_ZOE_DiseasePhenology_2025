# ZOE disease phenology analysis

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                     ZOE extra script (case study regions)              #
# ---------------------------------------------------------------------- #

# Load needed packages
library(terra)
library(tidyverse)
library(ggplot2)
library(ggh4x)





#-------------------------------------------------------------------------------

# 4. Calculate past decadal trends  --------------------------------------------
# for main vectors as well as viruses for the case study regions
# Slovenia and Slovakia

# Read in the monthly predictions under factual climate and land use change
# from 1970 to 2019
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_1970_2019.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")

# TBE
TBE_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_1970_2019.tif")

# WNV
WNV_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_1970_2019.tif")



# Read in the masks containing the countries with a mandatory reporting system of 
# disease surveillance data to the ECDC; we load different masks depending on
# the virus as Austria did not report of NUTS3 level for TBE
eu_eea_mask_WNV <- terra::rast("input_data/spatial_data/eu_eea_mask_WNV.tif") 
eu_eea_mask_TBE <- terra::rast("input_data/spatial_data/eu_eea_mask_TBE.tif")

# Read in the shapefiles from the case study regions Slovenia and Slovakia
shp_slovenia <- terra::vect("input_data/spatial_data/si_1km.shp")
shp_slovakia <- terra::vect("input_data/spatial_data/sk_1km.shp")

# Reproject the shapefile to the used crs (WGS84)
shp_slovenia <- project(shp_slovenia, "EPSG:4326")
shp_slovakia <- project(shp_slovakia, "EPSG:4326")

# Convert these shapefiles into rasters with a 50km resolution
r_slovenia <- rasterize(shp_slovenia, eu_eea_mask_WNV)
r_slovakia <- rasterize(shp_slovakia, eu_eea_mask_WNV)

# Crop the study region rasters to match the extent of the prediction rasters
# based on one example raster
r_slovenia <- crop(r_slovenia, I_ricinus_occ_prob_clim_landuse_ens)
r_slovakia <- crop(r_slovakia, I_ricinus_occ_prob_clim_landuse_ens)


# Mask the values of the prediction rasters that lay outside of these countries
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_slovenia <- mask(I_ricinus_occ_prob_clim_landuse_ens, r_slovenia) # Case study region Slovenia
I_ricinus_occ_prob_clim_landuse_ens_slovakia <- mask(I_ricinus_occ_prob_clim_landuse_ens, r_slovakia) # Case study region Slovakia

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_slovenia <- mask(C_pipiens_occ_prob_clim_landuse_ens, r_slovenia)
C_pipiens_occ_prob_clim_landuse_ens_slovakia <- mask(C_pipiens_occ_prob_clim_landuse_ens, r_slovakia) 

# TBE
TBE_occ_prob_clim_ens_slovenia <- mask(TBE_occ_prob_clim_ens, r_slovenia) 
TBE_occ_prob_clim_ens_slovakia <- mask(TBE_occ_prob_clim_ens, r_slovakia)

# WNV
WNV_occ_prob_clim_ens_slovenia <- mask(WNV_occ_prob_clim_ens, r_slovenia) 
WNV_occ_prob_clim_ens_slovakia <- mask(WNV_occ_prob_clim_ens, r_slovakia)


# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")


for (o in operations) { # Loop over the 95th percentile and mean functions
  
  print(o)
  
  # Generate time information
  dates <- seq(as.Date("1970-01-01"), as.Date("2019-12-01"), by = "month")
  
  # Assign dates as layer names
  names(I_ricinus_occ_prob_clim_landuse_ens_slovenia) <- dates
  names(I_ricinus_occ_prob_clim_landuse_ens_slovakia) <- dates
  names(C_pipiens_occ_prob_clim_landuse_ens_slovenia) <- dates
  names(C_pipiens_occ_prob_clim_landuse_ens_slovakia) <- dates
  names(TBE_occ_prob_clim_ens_slovenia) <- dates
  names(TBE_occ_prob_clim_ens_slovakia) <- dates
  names(WNV_occ_prob_clim_ens_slovenia) <- dates
  names(WNV_occ_prob_clim_ens_slovakia) <- dates
  
  # Define function for 95th percentile and mean with removing NAs
  percentile_95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
  mean_na_rm <- function(x) mean(x, na.rm = TRUE)
  
  
  # Choose the appropriate function
  if (o == "Peak") { fun <- percentile_95 
  } else if (o == "Mean") { fun <- mean_na_rm
  }
  
  # Calculate mean or 95th percentile (peak) occurrence probability across the case study regions for each month
  # and create a data frame with the calculated information (under the factual scenario)
  I_ricinus_monthly_clim_landuse_slovenia <- global(I_ricinus_occ_prob_clim_landuse_ens_slovenia, fun = fun)
  I_ricinus_monthly_clim_landuse_slovenia <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_slovenia[,1], scenario = "Factual prediction")
  
  I_ricinus_monthly_clim_landuse_slovakia <- global(I_ricinus_occ_prob_clim_landuse_ens_slovakia, fun = fun)
  I_ricinus_monthly_clim_landuse_slovakia <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_slovakia[,1], scenario = "Factual prediction")
  
  C_pipiens_monthly_clim_landuse_slovenia <- global(C_pipiens_occ_prob_clim_landuse_ens_slovenia, fun = fun)
  C_pipiens_monthly_clim_landuse_slovenia <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_slovenia[,1], scenario = "Factual prediction")
  
  C_pipiens_monthly_clim_landuse_slovakia <- global(C_pipiens_occ_prob_clim_landuse_ens_slovakia, fun = fun)
  C_pipiens_monthly_clim_landuse_slovakia <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_slovakia[,1], scenario = "Factual prediction")
  
  TBE_monthly_clim_slovenia <- global(TBE_occ_prob_clim_ens_slovenia, fun = fun)
  TBE_monthly_clim_slovenia <- data.frame(date = dates, occurrence = TBE_monthly_clim_slovenia[,1], scenario = "Factual prediction")
  
  TBE_monthly_clim_slovakia <- global(TBE_occ_prob_clim_ens_slovakia, fun = fun)
  TBE_monthly_clim_slovakia <- data.frame(date = dates, occurrence = TBE_monthly_clim_slovakia[,1], scenario = "Factual prediction")
  
  WNV_monthly_clim_slovenia <- global(WNV_occ_prob_clim_ens_slovenia, fun = fun)
  WNV_monthly_clim_slovenia <- data.frame(date = dates, occurrence = WNV_monthly_clim_slovenia[,1], scenario = "Factual prediction")
  
  WNV_monthly_clim_slovakia <- global(WNV_occ_prob_clim_ens_slovakia, fun = fun)
  WNV_monthly_clim_slovakia <- data.frame(date = dates, occurrence = WNV_monthly_clim_slovakia[,1], scenario = "Factual prediction")
  
  
  # Add another column to the data frames indicating the species/pathogen
  I_ricinus_monthly_clim_landuse_slovenia$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_landuse_slovakia$species <- "Ixodes ricinus"
  C_pipiens_monthly_clim_landuse_slovenia$species <- "Culex pipiens"
  C_pipiens_monthly_clim_landuse_slovakia$species <- "Culex pipiens"
  TBE_monthly_clim_slovenia$species <- "TBE"
  TBE_monthly_clim_slovakia$species <- "TBE"
  WNV_monthly_clim_slovenia$species <- "WNV"
  WNV_monthly_clim_slovakia$species <- "WNV"
  
  
  # Combine the data frames for the two case study regions
  combined_df_slovenia <- bind_rows(I_ricinus_monthly_clim_landuse_slovenia,
                                    C_pipiens_monthly_clim_landuse_slovenia,
                                    TBE_monthly_clim_slovenia,
                                    WNV_monthly_clim_slovenia)
  
  combined_df_slovakia <- bind_rows(I_ricinus_monthly_clim_landuse_slovakia,
                                    C_pipiens_monthly_clim_landuse_slovakia,
                                    TBE_monthly_clim_slovakia,
                                    WNV_monthly_clim_slovakia)
  
  
  
  # Create columns for month (1:12) and decade (1970s, 1980s, etc.) to group data
  combined_df_slovenia <- combined_df_slovenia %>%
    mutate(
      month = month(date),
      decade = case_when(
        year(date) >= 1970 & year(date) < 1980 ~ "1970s",
        year(date) >= 1980 & year(date) < 1990 ~ "1980s",
        year(date) >= 1990 & year(date) < 2000 ~ "1990s",
        year(date) >= 2000 & year(date) < 2010 ~ "2000s",
        year(date) >= 2010 & year(date) < 2020 ~ "2010s"
      )
    )
  
  combined_df_slovakia <- combined_df_slovakia %>%
    mutate(
      month = month(date),
      decade = case_when(
        year(date) >= 1970 & year(date) < 1980 ~ "1970s",
        year(date) >= 1980 & year(date) < 1990 ~ "1980s",
        year(date) >= 1990 & year(date) < 2000 ~ "1990s",
        year(date) >= 2000 & year(date) < 2010 ~ "2000s",
        year(date) >= 2010 & year(date) < 2020 ~ "2010s"
      )
    )
  
  # Aggregate by decade and month to calculate the mean for each time window (decade), month,
  # scenario and species (vector, virus)
  aggregated_df_slovenia <- combined_df_slovenia %>%
    group_by(decade, month, scenario, species) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  aggregated_df_slovakia <- combined_df_slovakia %>%
    group_by(decade, month, scenario, species) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  # Change month numbers to month name abbreviations
  aggregated_df_past_slovenia <- aggregated_df_slovenia %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  aggregated_df_past_slovakia <- aggregated_df_slovakia %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_past_slovenia <- aggregated_df_past_slovenia %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  aggregated_df_past_slovakia <- aggregated_df_past_slovakia %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Make sure the months are correctly ordered from January to December
  aggregated_df_past_slovenia$month <- factor(aggregated_df_past_slovenia$month, 
                                              levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  aggregated_df_past_slovakia$month <- factor(aggregated_df_past_slovakia$month, 
                                              levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                         "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vector and diseases appear in the correct position
  aggregated_df_past_slovenia$species <- factor(aggregated_df_past_slovenia$species, 
                                                levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  aggregated_df_past_slovakia$species <- factor(aggregated_df_past_slovakia$species, 
                                                levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  
  # Save the resulting data frames
  save(aggregated_df_past_slovenia, file = paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,"_slovenia.RData"))
  save(aggregated_df_past_slovakia, file = paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,"_slovakia.RData"))
  
  
} # Close the loop over the 95th percentile and mean functions




#-------------------------------------------------------------------------------

# 5. Calculate future decadal trends  ------------------------------------------
# for main vectors as well as viruses for the case study regions
# Slovenia and Slovakia

# Read in the monthly predictions under climate and land use change
# from 2030 to 2070 for environmental scenario ssp370
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp370.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp370.tif")

# TBE
TBE_occ_prob_clim_ens_fut_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp370.tif")

# WNV
WNV_occ_prob_clim_ens_fut_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp370.tif")




# Mask the values of the prediction rasters that lay outside of these countries
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovenia <- mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, r_slovenia) # Case study region Slovenia
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovakia <- mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, r_slovakia) # Case study region Slovakia

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovenia <- mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, r_slovenia)
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovakia <- mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, r_slovakia)

# TBE
TBE_occ_prob_clim_ens_fut_ssp370_slovenia <- mask(TBE_occ_prob_clim_ens_fut_ssp370, r_slovenia)
TBE_occ_prob_clim_ens_fut_ssp370_slovakia <- mask(TBE_occ_prob_clim_ens_fut_ssp370, r_slovakia)

# WNV
WNV_occ_prob_clim_ens_fut_ssp370_slovenia <- mask(WNV_occ_prob_clim_ens_fut_ssp370, r_slovenia)
WNV_occ_prob_clim_ens_fut_ssp370_slovakia <- mask(WNV_occ_prob_clim_ens_fut_ssp370, r_slovakia)

# Generate time information
dates <- seq(as.Date("2030-01-01"), as.Date("2070-12-01"), by = "month")

# Assign dates as layer names
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovenia) <- dates
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovakia) <- dates
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovenia) <- dates
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovakia) <- dates
names(TBE_occ_prob_clim_ens_fut_ssp370_slovenia) <- dates
names(TBE_occ_prob_clim_ens_fut_ssp370_slovakia) <- dates
names(WNV_occ_prob_clim_ens_fut_ssp370_slovenia) <- dates
names(WNV_occ_prob_clim_ens_fut_ssp370_slovakia) <- dates



# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")

for (o in operations) { # Loop over the 95th percentile and mean functions
  
  print(o)
  
  # Define function for 95th percentile and mean with removing NAs
  percentile_95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
  mean_na_rm <- function(x) mean(x, na.rm = TRUE)
  
  # Choose the appropriate function
  if (o == "Peak") { fun <- percentile_95 
  } else if (o == "Mean") { fun <- mean_na_rm
  }
  
  I_ricinus_monthly_clim_landuse_ssp370_slovenia <- global(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovenia, fun = fun)
  I_ricinus_monthly_clim_landuse_ssp370_slovenia <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_ssp370_slovenia[,1], scenario = "Factual prediction")
  
  I_ricinus_monthly_clim_landuse_ssp370_slovakia <- global(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370_slovakia, fun = fun)
  I_ricinus_monthly_clim_landuse_ssp370_slovakia <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_ssp370_slovakia[,1], scenario = "Factual prediction")
  
  C_pipiens_monthly_clim_landuse_ssp370_slovenia <- global(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovenia, fun = fun)
  C_pipiens_monthly_clim_landuse_ssp370_slovenia <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_ssp370_slovenia[,1], scenario = "Factual prediction")
  
  C_pipiens_monthly_clim_landuse_ssp370_slovakia <- global(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370_slovakia, fun = fun)
  C_pipiens_monthly_clim_landuse_ssp370_slovakia <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_ssp370_slovakia[,1], scenario = "Factual prediction")
  
  TBE_monthly_clim_ssp370_slovenia <- global(TBE_occ_prob_clim_ens_fut_ssp370_slovenia, fun = fun)
  TBE_monthly_clim_ssp370_slovenia <- data.frame(date = dates, occurrence = TBE_monthly_clim_ssp370_slovenia[,1], scenario = "Factual prediction")
  
  TBE_monthly_clim_ssp370_slovakia <- global(TBE_occ_prob_clim_ens_fut_ssp370_slovakia, fun = fun)
  TBE_monthly_clim_ssp370_slovakia <- data.frame(date = dates, occurrence = TBE_monthly_clim_ssp370_slovakia[,1], scenario = "Factual prediction")
  
  WNV_monthly_clim_ssp370_slovenia <- global(WNV_occ_prob_clim_ens_fut_ssp370_slovenia, fun = fun)
  WNV_monthly_clim_ssp370_slovenia <- data.frame(date = dates, occurrence = WNV_monthly_clim_ssp370_slovenia[,1], scenario = "Factual prediction")
  
  WNV_monthly_clim_ssp370_slovakia <- global(WNV_occ_prob_clim_ens_fut_ssp370_slovakia, fun = fun)
  WNV_monthly_clim_ssp370_slovakia <- data.frame(date = dates, occurrence = WNV_monthly_clim_ssp370_slovakia[,1], scenario = "Factual prediction")
  
  # Add another column to the data frames indicating the species/pathogen
  # as well as a column indicating that we look at future disease phenologies
  I_ricinus_monthly_clim_landuse_ssp370_slovenia$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_landuse_ssp370_slovakia$species <- "Ixodes ricinus"
  C_pipiens_monthly_clim_landuse_ssp370_slovenia$species <- "Culex pipiens"
  C_pipiens_monthly_clim_landuse_ssp370_slovakia$species <- "Culex pipiens"
  TBE_monthly_clim_ssp370_slovenia$species <- "TBE"
  TBE_monthly_clim_ssp370_slovakia$species <- "TBE"
  WNV_monthly_clim_ssp370_slovenia$species <- "WNV"
  WNV_monthly_clim_ssp370_slovakia$species <- "WNV"
  
  # Combine the data frames
  combined_df_slovenia <- bind_rows(I_ricinus_monthly_clim_landuse_ssp370_slovenia,
                                    C_pipiens_monthly_clim_landuse_ssp370_slovenia,
                                    TBE_monthly_clim_ssp370_slovenia,
                                    WNV_monthly_clim_ssp370_slovenia)
  
  combined_df_slovakia <- bind_rows(I_ricinus_monthly_clim_landuse_ssp370_slovakia,
                                    C_pipiens_monthly_clim_landuse_ssp370_slovakia,
                                    TBE_monthly_clim_ssp370_slovakia,
                                    WNV_monthly_clim_ssp370_slovakia)
  
  
  # Create columns for month (1:12) and decade (2030s, 2040s, etc.) to group data
  combined_df_slovenia <- combined_df_slovenia %>%
    mutate(
      month = month(date),
      decade = case_when(
        year(date) >= 2030 & year(date) < 2040 ~ "2030s",
        year(date) >= 2040 & year(date) < 2050 ~ "2040s",
        year(date) >= 2050 & year(date) < 2060 ~ "2050s",
        year(date) >= 2060 & year(date) < 2070 ~ "2060s",
      )
    )
  
  combined_df_slovakia <- combined_df_slovakia %>%
    mutate(
      month = month(date),
      decade = case_when(
        year(date) >= 2030 & year(date) < 2040 ~ "2030s",
        year(date) >= 2040 & year(date) < 2050 ~ "2040s",
        year(date) >= 2050 & year(date) < 2060 ~ "2050s",
        year(date) >= 2060 & year(date) < 2070 ~ "2060s",
      )
    )
  
  # Remove rows that conatin NA values (exclude 2070s)
  combined_df_slovenia <- na.omit(combined_df_slovenia)
  combined_df_slovakia <- na.omit(combined_df_slovakia)
  
  # Aggregate by decade and month to calculate the mean for each time window (decade) and month
  aggregated_df_fut_slovenia <- combined_df_slovenia %>%
    group_by(decade, month, scenario, species) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  aggregated_df_fut_slovakia <- combined_df_slovakia %>%
    group_by(decade, month, scenario, species) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  # Change month numbers to month name abbreviations
  aggregated_df_fut_slovenia <- aggregated_df_fut_slovenia %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  aggregated_df_fut_slovakia <- aggregated_df_fut_slovakia %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_fut_slovenia <- aggregated_df_fut_slovenia %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  aggregated_df_fut_slovakia <- aggregated_df_fut_slovakia %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Make sure the months are correctly ordered from January to Decemver
  aggregated_df_fut_slovenia$month <- factor(aggregated_df_fut_slovenia$month, 
                                             levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  aggregated_df_fut_slovakia$month <- factor(aggregated_df_fut_slovakia$month, 
                                             levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vector and diseases appear in the correct position
  aggregated_df_fut_slovenia$species <- factor(aggregated_df_fut_slovenia$species, 
                                               levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  aggregated_df_fut_slovakia$species <- factor(aggregated_df_fut_slovakia$species, 
                                               levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  
  
  # Save the resulting data frames
  save(aggregated_df_fut_slovenia, file = paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,"_slovenia.RData"))
  save(aggregated_df_fut_slovakia, file = paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,"_slovakia.RData"))
  
  
} # Close the loop over the 95th percentile and mean functions





#-------------------------------------------------------------------------------

# 3. Visualise decadal trends  -------------------------------------------------
# for main vectors as well as viruses for the case study region Slovenia and
# Slovakia


# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")

for (o in operations) { # Loop over the peak and mean functions
  
  print(o)
  
  # Load the needed data of past and future decadal trends
  load(paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,"_slovenia.RData"))
  load(paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,"_slovakia.RData"))
  load(paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,"_slovenia.RData"))
  load(paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,"_slovakia.RData"))
  
  # Bind the two data frames containing the information of past and future decadal
  # trends of Ixodes ricinus, TBE, Culex pipiens, and WNV for each case study regions
  decadal_trends_past_fut_slovenia <- rbind(aggregated_df_past_slovenia, aggregated_df_fut_slovenia)
  decadal_trends_past_fut_slovakia <- rbind(aggregated_df_past_slovakia, aggregated_df_fut_slovakia)
  
  # Make sure the months are correctly ordered from January to December
  decadal_trends_past_fut_slovenia$month <- factor(decadal_trends_past_fut_slovenia$month, 
                                                   levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  decadal_trends_past_fut_slovakia$month <- factor(decadal_trends_past_fut_slovakia$month, 
                                                   levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vector and diseases appear in the correct position
  decadal_trends_past_fut_slovenia$species <- factor(decadal_trends_past_fut_slovenia$species, 
                                                     levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  decadal_trends_past_fut_slovakiaspecies <- factor(decadal_trends_past_fut_slovakia$species, 
                                                    levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  # Subset the data frames to only contain prediction results from the
  # decades of the 1970s, 2010s, 2050s
  decadal_trends_past_fut_slovenia <- decadal_trends_past_fut_slovenia[decadal_trends_past_fut_slovenia$decade %in% c("1970s", "2010s", "2050s"), ]
  decadal_trends_past_fut_slovakia <- decadal_trends_past_fut_slovakia[decadal_trends_past_fut_slovakia$decade %in% c("1970s", "2010s", "2050s"), ]
  
  # Add a column containing the region name
  decadal_trends_past_fut_slovenia$region <- "Slovenia"
  decadal_trends_past_fut_slovakia$region <- "Slovakia"
  
  
  
  # Visualize the data for future environmental scenario ssp370
  # for Ixodes ricinus and TBE in the case study region Slovenia
  ggplot(data = decadal_trends_past_fut_slovenia[decadal_trends_past_fut_slovenia$species %in% c("Ixodes ricinus", "TBE"), ], 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(region), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3")           
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        #"1990s" = "midnightblue",
        #"2000s" = "seagreen",
        "2010s" = "seagreen",
        #"2030s" = "#FB9B06",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 12),  
      axis.title = element_text(size = 12.5),  
      axis.text = element_text(size = 10),  
      legend.title = element_text(size = 10.5, face = "bold"),  
      legend.text = element_text(size = 10),  
      plot.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE),  
      linetype = "none" 
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
  
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Ixodes_TBE_ssp370_",o,"_slovenia.png"), width = 5.5, height = 5)
  
  
  # Visualize the data for future environmental scenario ssp370
  # for Ixodes ricinus and TBE in the case study region Slovakia
  ggplot(data = decadal_trends_past_fut_slovakia[decadal_trends_past_fut_slovakia$species %in% c("Ixodes ricinus", "TBE"), ], 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(region), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3")           
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        #"1990s" = "midnightblue",
        #"2000s" = "seagreen",
        "2010s" = "seagreen",
        #"2030s" = "#FB9B06",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 12),  
      axis.title = element_text(size = 12.5),  
      axis.text = element_text(size = 10),  
      legend.title = element_text(size = 10.5, face = "bold"),  
      legend.text = element_text(size = 10),  
      plot.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE),  
      linetype = "none" 
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
  
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Ixodes_TBE_ssp370_",o,"_slovakia.png"), width = 5.5, height = 5)
  
  
  # Visualize the data for future environmental scenario ssp370
  # for Culex pipiens and WNV in the case study region Slovenia
  ggplot(data = decadal_trends_past_fut_slovenia[decadal_trends_past_fut_slovenia$species %in% c("Culex pipiens", "WNV"), ], 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(region), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        #"1990s" = "midnightblue",
        #"2000s" = "seagreen",
        "2010s" = "seagreen",
        #"2030s" = "#FB9B06",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 12),  
      axis.title = element_text(size = 12.5),  
      axis.text = element_text(size = 10),  
      legend.title = element_text(size = 10.5, face = "bold"),  
      legend.text = element_text(size = 10),  
      plot.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE),  
      linetype = "none" 
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
  
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Culex_WNV_ssp370_",o,"_slovenia.png"), width = 5.5, height = 5)
  
  
  # Visualize the data for future environmental scenario ssp370
  # for Culex pipiens and WNV in the case study region Slovakia
  ggplot(data = decadal_trends_past_fut_slovakia[decadal_trends_past_fut_slovakia$species %in% c("Culex pipiens", "WNV"), ], 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(region), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        #"1990s" = "midnightblue",
        #"2000s" = "seagreen",
        "2010s" = "seagreen",
        #"2030s" = "#FB9B06",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 12),  
      axis.title = element_text(size = 12.5),  
      axis.text = element_text(size = 10),  
      legend.title = element_text(size = 10.5, face = "bold"),  
      legend.text = element_text(size = 10),  
      plot.title = element_text(size = 13, face = "bold"),
      strip.text = element_text(size = 13, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE),  
      linetype = "none" 
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
  
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Culex_WNV_ssp370_",o,"_slovakia.png"), width = 5.5, height = 5)
  
  
} # Close the loop over the peak and mean functions
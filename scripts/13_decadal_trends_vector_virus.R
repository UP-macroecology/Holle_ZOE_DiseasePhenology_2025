# ZOE disease phenology analysis

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                         13. Decadal trends                             #
# ---------------------------------------------------------------------- #

# Load needed packages
library(terra)
library(tidyverse)
library(ggplot2)
library(ggh4x)




#-------------------------------------------------------------------------------

# 1. Calculate past decadal trends  --------------------------------------------
# for main vectors as well as viruses


# a) Load data -----------------------------------------------------------------

# Load needed data - postprocessed monthly predicted ensemble occurrence probability
# of main vectors and the respective viruses, under climate and land use change, as 
# well as under the counterfactual scenarios for the years 1970 to 2019
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_1970_2019.tif")
I_ricinus_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_landuse_ens_1970_2019.tif")
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_nolanduse_ens_1970_2019.tif")
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_nolanduse_ens_1970_2019.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif")

# TBE
TBE_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_1970_2019.tif")
TBE_occ_prob_noclim_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_noclim_ens_1970_2019.tif")
TBE_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_nolanduse_ens_1970_2019.tif")
TBE_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_noclim_nolanduse_ens_1970_2019.tif")

# WNV
WNV_occ_prob_clim_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_1970_2019.tif")
WNV_occ_prob_noclim_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_noclim_ens_1970_2019.tif")
WNV_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_nolanduse_ens_1970_2019.tif")
WNV_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_noclim_nolanduse_ens_1970_2019.tif")




# b) Calculate decadal trends --------------------------------------------------

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")


for (o in operations) { # Loop over the 95th percentile and mean functions
  
  print(o)
  
  # Generate time information
  dates <- seq(as.Date("1970-01-01"), as.Date("2019-12-01"), by = "month")
  
  # Assign dates as layer names
  names(I_ricinus_occ_prob_clim_landuse_ens) <- dates
  names(I_ricinus_occ_prob_noclim_landuse_ens) <- dates
  names(I_ricinus_occ_prob_clim_nolanduse_ens) <- dates
  names(I_ricinus_occ_prob_noclim_nolanduse_ens) <- dates
  names(C_pipiens_occ_prob_clim_landuse_ens) <- dates
  names(C_pipiens_occ_prob_noclim_landuse_ens) <- dates
  names(C_pipiens_occ_prob_clim_nolanduse_ens) <- dates
  names(C_pipiens_occ_prob_noclim_nolanduse_ens) <- dates
  names(TBE_occ_prob_clim_ens) <- dates
  names(TBE_occ_prob_noclim_ens) <- dates
  names(TBE_occ_prob_clim_nolanduse_ens) <- dates
  names(TBE_occ_prob_noclim_nolanduse_ens) <- dates
  names(WNV_occ_prob_clim_ens) <- dates
  names(WNV_occ_prob_noclim_ens) <- dates
  names(WNV_occ_prob_clim_nolanduse_ens) <- dates
  names(WNV_occ_prob_noclim_nolanduse_ens) <- dates
  
  # Define function for 95th percentile and mean with removing NAs
  percentile_95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
  mean_na_rm <- function(x) mean(x, na.rm = TRUE)
  
  
  # Choose the appropriate function
  if (o == "Peak") { fun <- percentile_95 
  } else if (o == "Mean") { fun <- mean_na_rm
  }
  
  # Calculate mean or 95th percentile (peak) occurrence probability across the study area (Europe) for each month
  # and create a data frame with the calculated information (under the factual and counterfactual scenarios)
  print("Calculating - Ixodes ricinus; Factual prediction")
  I_ricinus_monthly_clim_landuse <- global(I_ricinus_occ_prob_clim_landuse_ens, fun = fun)
  I_ricinus_monthly_clim_landuse <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse[,1], scenario = "Factual prediction")
  
  print("Calculating - Ixodes ricinus; Counterfactual climate")
  I_ricinus_monthly_noclim_landuse <- global(I_ricinus_occ_prob_noclim_landuse_ens, fun = fun)
  I_ricinus_monthly_noclim_landuse <- data.frame(date = dates, occurrence = I_ricinus_monthly_noclim_landuse[,1], scenario = "Counterfactual climate")
  
  print("Calculating - Ixodes ricinus; Counterfactual land use")
  I_ricinus_monthly_clim_nolanduse <- global(I_ricinus_occ_prob_clim_nolanduse_ens, fun = fun)
  I_ricinus_monthly_clim_nolanduse <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_nolanduse[,1], scenario = "Counterfactual land use")

  print("Calculating - Ixodes ricinus; Counterfactual climate + land use")
  I_ricinus_monthly_noclim_nolanduse <- global(I_ricinus_occ_prob_noclim_nolanduse_ens, fun = fun)
  I_ricinus_monthly_noclim_nolanduse <- data.frame(date = dates, occurrence = I_ricinus_monthly_noclim_nolanduse[,1], scenario = "Counterfactual climate + land use")
  
  print("Calculating - Culex pipiens; Factual prediction")
  C_pipiens_monthly_clim_landuse <- global(C_pipiens_occ_prob_clim_landuse_ens, fun = fun)
  C_pipiens_monthly_clim_landuse <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse[,1], scenario = "Factual prediction")
  
  print("Calculating - Culex pipiens; Counterfactual climate")
  C_pipiens_monthly_noclim_landuse <- global(C_pipiens_occ_prob_noclim_landuse_ens, fun = fun)
  C_pipiens_monthly_noclim_landuse <- data.frame(date = dates, occurrence = C_pipiens_monthly_noclim_landuse[,1], scenario = "Counterfactual climate")
  
  print("Calculating - Culex pipiens; Counterfactual land use")
  C_pipiens_monthly_clim_nolanduse <- global(C_pipiens_occ_prob_clim_nolanduse_ens, fun = fun)
  C_pipiens_monthly_clim_nolanduse <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_nolanduse[,1], scenario = "Counterfactual land use")

  print("Calculating - Culex pipiens; Counterfactual climate + land use")
  C_pipiens_monthly_noclim_nolanduse <- global(C_pipiens_occ_prob_noclim_nolanduse_ens, fun = fun)
  C_pipiens_monthly_noclim_nolanduse <- data.frame(date = dates, occurrence = C_pipiens_monthly_noclim_nolanduse[,1], scenario = "Counterfactual climate + land use")
  
  print("Calculating - TBE; Factual prediction")
  TBE_monthly_clim <- global(TBE_occ_prob_clim_ens, fun = fun)
  TBE_monthly_clim <- data.frame(date = dates, occurrence = TBE_monthly_clim[,1], scenario = "Factual prediction")
  
  print("Calculating - TBE; Counterfactual climate")
  TBE_monthly_noclim <- global(TBE_occ_prob_noclim_ens, fun = fun)
  TBE_monthly_noclim <- data.frame(date = dates, occurrence = TBE_monthly_noclim[,1], scenario = "Counterfactual climate")
  
  print("Calculating - TBE; Counterfactual land use")
  TBE_monthly_clim_nolanduse <- global(TBE_occ_prob_clim_nolanduse_ens, fun = fun)
  TBE_monthly_clim_nolanduse <- data.frame(date = dates, occurrence = TBE_monthly_clim_nolanduse[,1], scenario = "Counterfactual land use")

  print("Calculating - TBE; Counterfactual climate + land use")
  TBE_monthly_noclim_nolanduse <- global(TBE_occ_prob_noclim_nolanduse_ens, fun = fun)
  TBE_monthly_noclim_nolanduse <- data.frame(date = dates, occurrence = TBE_monthly_noclim_nolanduse[,1], scenario = "Counterfactual climate + land use")
  
  print("Calculating - WNV; Factual prediction")
  WNV_monthly_clim <- global(WNV_occ_prob_clim_ens, fun = fun)
  WNV_monthly_clim <- data.frame(date = dates, occurrence = WNV_monthly_clim[,1], scenario = "Factual prediction")
  
  print("Calculating - WNV; Counterfactual climate")
  WNV_monthly_noclim <- global(WNV_occ_prob_noclim_ens, fun = fun)
  WNV_monthly_noclim <- data.frame(date = dates, occurrence = WNV_monthly_noclim[,1], scenario = "Counterfactual climate")
  
  print("Calculating - WNV; Counterfactual land use")
  WNV_monthly_clim_nolanduse <- global(WNV_occ_prob_clim_nolanduse_ens, fun = fun)
  WNV_monthly_clim_nolanduse <- data.frame(date = dates, occurrence = WNV_monthly_clim_nolanduse[,1], scenario = "Counterfactual land use")

  print("Calculating - WNV; Counterfactual climate + land use")
  WNV_monthly_noclim_nolanduse <- global(WNV_occ_prob_noclim_nolanduse_ens, fun = fun)
  WNV_monthly_noclim_nolanduse <- data.frame(date = dates, occurrence = WNV_monthly_noclim_nolanduse[,1], scenario = "Counterfactual climate + land use")
  
  # Add another column to the data frames indicating the species/pathogen
  I_ricinus_monthly_clim_landuse$species <- "Ixodes ricinus"
  I_ricinus_monthly_noclim_landuse$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_nolanduse$species <- "Ixodes ricinus"
  I_ricinus_monthly_noclim_nolanduse$species <- "Ixodes ricinus"
  C_pipiens_monthly_clim_landuse$species <- "Culex pipiens"
  C_pipiens_monthly_noclim_landuse$species <- "Culex pipiens"
  C_pipiens_monthly_clim_nolanduse$species <- "Culex pipiens"
  C_pipiens_monthly_noclim_nolanduse$species <- "Culex pipiens"
  TBE_monthly_clim$species <- "TBE"
  TBE_monthly_noclim$species <- "TBE"
  TBE_monthly_clim_nolanduse$species <- "TBE"
  TBE_monthly_noclim_nolanduse$species <- "TBE"
  WNV_monthly_clim$species <- "WNV"
  WNV_monthly_noclim$species <- "WNV"
  WNV_monthly_clim_nolanduse$species <- "WNV"
  WNV_monthly_noclim_nolanduse$species <- "WNV"
  
  # Combine the data frames
  combined_df <- bind_rows(I_ricinus_monthly_clim_landuse, I_ricinus_monthly_noclim_landuse,
                           I_ricinus_monthly_clim_nolanduse, I_ricinus_monthly_noclim_nolanduse,
                           C_pipiens_monthly_clim_landuse, C_pipiens_monthly_noclim_landuse,
                           C_pipiens_monthly_clim_nolanduse, C_pipiens_monthly_noclim_nolanduse,
                           TBE_monthly_clim, TBE_monthly_noclim,
                           TBE_monthly_clim_nolanduse, TBE_monthly_noclim_nolanduse,
                           WNV_monthly_clim, WNV_monthly_noclim,
                           WNV_monthly_clim_nolanduse, WNV_monthly_noclim_nolanduse)
  
  # Create columns for month (1:12) and decade (1970s, 1980s, etc.) to group data
  combined_df <- combined_df %>%
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
  aggregated_df <- combined_df %>%
    group_by(decade, month, scenario, species) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  
  
# c) Prepare data frame for plotting -------------------------------------------
  
  # Change month numbers to month name abbreviations
  aggregated_df_past <- aggregated_df %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  # Add a column that indicate that the predictions results are based on 
  # historical data
  aggregated_df_past$time <- "Historical phenology"
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_past <- aggregated_df_past %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Make sure the months are correctly ordered from January to December
  aggregated_df_past$month <- factor(aggregated_df_past$month, 
                                     levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vector and diseases appear in the correct position
  aggregated_df_past$species <- factor(aggregated_df_past$species, 
                                       levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  
  
  
# d) Save resulting data frame -------------------------------------------------
  
  # Save the data frame containing the monthly mean and peak occurrence probabilities
  # per decade
  save(aggregated_df_past, file = paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,".RData"))
  
  
} # Close the loop over the 95th percentile and mean functions




#-------------------------------------------------------------------------------

# 2. Calculate future decadal trends  ------------------------------------------
# for main vectors as well as viruses


# a) Load data -----------------------------------------------------------------

# Load needed data - postprocessed monthly predicted ensemble occurrence probability of main vectors and
# the respective viruses, under climate and land use change for the
# future years 2030 to 2070 based and three studied environmental scenarios (summarising 5 climate models)
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



# b) Calculate decadal trends --------------------------------------------------

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")

# Generate time information
dates <- seq(as.Date("2030-01-01"), as.Date("2070-12-01"), by = "month")

# Assign dates as layer names
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126) <- dates
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370) <- dates
names(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585) <- dates
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126) <- dates
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370) <- dates
names(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585) <- dates
names(TBE_occ_prob_clim_ens_fut_ssp126) <- dates
names(TBE_occ_prob_clim_ens_fut_ssp370) <- dates
names(TBE_occ_prob_clim_ens_fut_ssp585) <- dates
names(WNV_occ_prob_clim_ens_fut_ssp126) <- dates
names(WNV_occ_prob_clim_ens_fut_ssp370) <- dates
names(WNV_occ_prob_clim_ens_fut_ssp585) <- dates



for (o in operations) { # Loop over the 95th percentile and mean functions
  
  print(o)
  
  # Define function for 95th percentile and mean with removing NAs
  percentile_95 <- function(x) quantile(x, probs = 0.95, na.rm = TRUE)
  mean_na_rm <- function(x) mean(x, na.rm = TRUE)
  
  # Choose the appropriate function
  if (o == "Peak") { fun <- percentile_95 
  } else if (o == "Mean") { fun <- mean_na_rm
  }
  
  print("Calculating - Ixodes ricinus; Factual prediction ssp126")
  I_ricinus_monthly_clim_landuse_ssp126 <- global(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, fun = fun)
  I_ricinus_monthly_clim_landuse_ssp126 <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_ssp126[,1], scenario = "Factual prediction")
  
  print("Calculating - Ixodes ricinus; Factual prediction ssp370")
  I_ricinus_monthly_clim_landuse_ssp370 <- global(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, fun = fun)
  I_ricinus_monthly_clim_landuse_ssp370 <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_ssp370[,1], scenario = "Factual prediction")
  
  print("Calculating - Ixodes ricinus; Factual prediction ssp585")
  I_ricinus_monthly_clim_landuse_ssp585 <- global(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, fun = fun)
  I_ricinus_monthly_clim_landuse_ssp585 <- data.frame(date = dates, occurrence = I_ricinus_monthly_clim_landuse_ssp585[,1], scenario = "Factual prediction")
  
  print("Calculating - Culex pipiens; Factual prediction ssp126")
  C_pipiens_monthly_clim_landuse_ssp126 <- global(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, fun = fun)
  C_pipiens_monthly_clim_landuse_ssp126 <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_ssp126[,1], scenario = "Factual prediction")
  
  print("Calculating - Culex pipiens; Factual prediction ssp370")
  C_pipiens_monthly_clim_landuse_ssp370 <- global(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, fun = fun)
  C_pipiens_monthly_clim_landuse_ssp370 <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_ssp370[,1], scenario = "Factual prediction")
  
  print("Calculating - Culex pipiens; Factual prediction ssp585")
  C_pipiens_monthly_clim_landuse_ssp585 <- global(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, fun = fun)
  C_pipiens_monthly_clim_landuse_ssp585 <- data.frame(date = dates, occurrence = C_pipiens_monthly_clim_landuse_ssp585[,1], scenario = "Factual prediction")
  
  print("Calculating - TBE; Factual prediction ssp126")
  TBE_monthly_clim_ssp126 <- global(TBE_occ_prob_clim_ens_fut_ssp126, fun = fun)
  TBE_monthly_clim_ssp126 <- data.frame(date = dates, occurrence = TBE_monthly_clim_ssp126[,1], scenario = "Factual prediction")
  
  print("Calculating - TBE; Factual prediction ssp370")
  TBE_monthly_clim_ssp370 <- global(TBE_occ_prob_clim_ens_fut_ssp370, fun = fun)
  TBE_monthly_clim_ssp370 <- data.frame(date = dates, occurrence = TBE_monthly_clim_ssp370[,1], scenario = "Factual prediction")
  
  print("Calculating - TBE; Factual prediction ssp585")
  TBE_monthly_clim_ssp585 <- global(TBE_occ_prob_clim_ens_fut_ssp585, fun = fun)
  TBE_monthly_clim_ssp585 <- data.frame(date = dates, occurrence = TBE_monthly_clim_ssp585[,1], scenario = "Factual prediction")
  
  print("Calculating - WNV; Factual prediction ssp126")
  WNV_monthly_clim_ssp126 <- global(WNV_occ_prob_clim_ens_fut_ssp126, fun = fun)
  WNV_monthly_clim_ssp126 <- data.frame(date = dates, occurrence = WNV_monthly_clim_ssp126[,1], scenario = "Factual prediction")
  
  print("Calculating - WNV; Factual prediction ssp370")
  WNV_monthly_clim_ssp370 <- global(WNV_occ_prob_clim_ens_fut_ssp370, fun = fun)
  WNV_monthly_clim_ssp370 <- data.frame(date = dates, occurrence = WNV_monthly_clim_ssp370[,1], scenario = "Factual prediction")
  
  print("Calculating - WNV; Factual prediction ssp585")
  WNV_monthly_clim_ssp585 <- global(WNV_occ_prob_clim_ens_fut_ssp585, fun = fun)
  WNV_monthly_clim_ssp585 <- data.frame(date = dates, occurrence = WNV_monthly_clim_ssp585[,1], scenario = "Factual prediction")
  
  # Add another column to the data frames indicating the species/pathogen
  # as well as a column indicating that we look at future disease phenologies
  I_ricinus_monthly_clim_landuse_ssp126$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_landuse_ssp126$time <- "Future phenology; ssp126"
  I_ricinus_monthly_clim_landuse_ssp370$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_landuse_ssp370$time <- "Future phenology; ssp370"
  I_ricinus_monthly_clim_landuse_ssp585$species <- "Ixodes ricinus"
  I_ricinus_monthly_clim_landuse_ssp585$time <- "Future phenology; ssp585"
  C_pipiens_monthly_clim_landuse_ssp126$species <- "Culex pipiens"
  C_pipiens_monthly_clim_landuse_ssp126$time <- "Future phenology; ssp126"
  C_pipiens_monthly_clim_landuse_ssp370$species <- "Culex pipiens"
  C_pipiens_monthly_clim_landuse_ssp370$time <- "Future phenology; ssp370"
  C_pipiens_monthly_clim_landuse_ssp585$species <- "Culex pipiens"
  C_pipiens_monthly_clim_landuse_ssp585$time <- "Future phenology; ssp585"
  TBE_monthly_clim_ssp126$species <- "TBE"
  TBE_monthly_clim_ssp126$time <- "Future phenology; ssp126"
  TBE_monthly_clim_ssp370$species <- "TBE"
  TBE_monthly_clim_ssp370$time <- "Future phenology; ssp370"
  TBE_monthly_clim_ssp585$species <- "TBE"
  TBE_monthly_clim_ssp585$time <- "Future phenology; ssp585"
  WNV_monthly_clim_ssp126$species <- "WNV"
  WNV_monthly_clim_ssp126$time <- "Future phenology; ssp126"
  WNV_monthly_clim_ssp370$species <- "WNV"
  WNV_monthly_clim_ssp370$time <- "Future phenology; ssp370"
  WNV_monthly_clim_ssp585$species <- "WNV"
  WNV_monthly_clim_ssp585$time <- "Future phenology; ssp585"
  
  # Combine the data frames
  combined_df <- bind_rows(I_ricinus_monthly_clim_landuse_ssp126, I_ricinus_monthly_clim_landuse_ssp370, I_ricinus_monthly_clim_landuse_ssp585,
                           C_pipiens_monthly_clim_landuse_ssp126, C_pipiens_monthly_clim_landuse_ssp370, C_pipiens_monthly_clim_landuse_ssp585,
                           TBE_monthly_clim_ssp126, TBE_monthly_clim_ssp370, TBE_monthly_clim_ssp585, 
                           WNV_monthly_clim_ssp126, WNV_monthly_clim_ssp370, WNV_monthly_clim_ssp585)
  
  
  # Create columns for month (1:12) and decade (2030s, 2040s, etc.) to group data
  combined_df <- combined_df %>%
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
  combined_df <- na.omit(combined_df)
  
  
  
  # Aggregate by decade and month to calculate the mean for each time window (decade) and month
  aggregated_df_fut <- combined_df %>%
    group_by(decade, month, scenario, species, time) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  
# c) Prepare data frame for plotting -------------------------------------------
  
  # Change month numbers to month name abbreviations
  aggregated_df_fut <- aggregated_df_fut %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  # Filter for months May to December for WNV (months with actual infection occurrences)
  aggregated_df_fut <- aggregated_df_fut %>%
    mutate(month = if_else(species == "WNV" & !(month %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
                           NA_character_, month)) %>%
    filter(!is.na(month))
  
  # Make sure the months are correctly ordered from January to Decemver
  aggregated_df_fut$month <- factor(aggregated_df_fut$month, 
                                    levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                                               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  
  # Make sure the vector and diseases appear in the correct position
  aggregated_df_fut$species <- factor(aggregated_df_fut$species, 
                                      levels = c("Ixodes ricinus", "Culex pipiens", "TBE", "WNV"))
  
  
  
# d) Save resulting data frame -------------------------------------------------
  
  # Save the data frame containing the monthly mean and peak occurrence probabilities
  # per decade
  save(aggregated_df_fut, file = paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,".RData"))
  
  
} # Close the loop over the 95th percentile and mean functions





#-------------------------------------------------------------------------------

# 3. Visualise decadal trends  -------------------------------------------------
# for main vectors as well as viruses


# a) Prepare data frame for visualisation --------------------------------------

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("Peak", "Mean")

for (o in operations) { # Loop over the peak and mean functions
  
  print(o)
  
  # Load the needed data of past and future decadal trends
  load(paste0("output_data/results/decadal_trends/decadal_trends_past_vector_virus_",o,".RData"))
  load(paste0("output_data/results/decadal_trends/decadal_trends_fut_vector_virus_",o,".RData"))
  
  # Subset the summarised predictions of the decade of the 2010s
  aggregated_df_past_2010s <- aggregated_df_past[aggregated_df_past$decade == "2010s" &
                                                   aggregated_df_past$scenario == "Factual prediction", ]
  
  # Add these data to also be represented in the future phenology
  aggregated_df_past_2010sssp126 <- aggregated_df_past_2010s
  aggregated_df_past_2010sssp126$time <- "Future phenology; ssp126"
  aggregated_df_past_2010sssp370 <- aggregated_df_past_2010s
  aggregated_df_past_2010sssp370$time <- "Future phenology; ssp370"
  aggregated_df_past_2010sssp585 <- aggregated_df_past_2010s
  aggregated_df_past_2010sssp585$time <- "Future phenology; ssp585"
  
  # Retain future prediction results based on the three different future env.
  # scenarios
  aggregated_df_futssp126 <- aggregated_df_fut[aggregated_df_fut$time == "Future phenology; ssp126", ]
  aggregated_df_futssp370 <- aggregated_df_fut[aggregated_df_fut$time == "Future phenology; ssp370", ]
  aggregated_df_futssp585 <- aggregated_df_fut[aggregated_df_fut$time == "Future phenology; ssp585", ]
  
  # Add the summarised predictions of the decades of the 2010s to the
  # correpsonding data frames of the different env. scnearios
  aggregated_df_futssp126 <- rbind(aggregated_df_futssp126, aggregated_df_past_2010sssp126)
  aggregated_df_futssp370 <- rbind(aggregated_df_futssp370, aggregated_df_past_2010sssp370)
  aggregated_df_futssp585 <- rbind(aggregated_df_futssp585, aggregated_df_past_2010sssp585)
  
  
  # Bind the two data frames containing the information of past and future decadal
  # trends of Ixodes ricinus, TBE, Culex pipiens, and WNV
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
  
  # Make sure the vector and diseases appear in the correct position
  decadal_trends_past_futssp126$species <- factor(decadal_trends_past_futssp126$species, 
                                            levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  decadal_trends_past_futssp370$species <- factor(decadal_trends_past_futssp370$species, 
                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  decadal_trends_past_futssp585$species <- factor(decadal_trends_past_futssp585$species, 
                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
  
  # Make sure the examined time ranges appear in the correct order
  decadal_trends_past_futssp126$time <- factor(decadal_trends_past_futssp126$time, 
                                                  levels = c("Historical phenology", "Future phenology; ssp126"))
  
  decadal_trends_past_futssp370$time <- factor(decadal_trends_past_futssp370$time, 
                                               levels = c("Historical phenology", "Future phenology; ssp370"))
  
  decadal_trends_past_futssp585$time <- factor(decadal_trends_past_futssp585$time, 
                                               levels = c("Historical phenology", "Future phenology; ssp585"))
  
  # Make sure the used scenarios are named in the correct order
  decadal_trends_past_futssp126$scenario <- factor(decadal_trends_past_futssp126$scenario, 
                                               levels = c("Factual prediction", "Counterfactual climate", "Counterfactual land use", 
                                                          "Counterfactual climate + land use"))
  
  decadal_trends_past_futssp370$scenario <- factor(decadal_trends_past_futssp370$scenario, 
                                                   levels = c("Factual prediction", "Counterfactual climate", "Counterfactual land use", 
                                                              "Counterfactual climate + land use"))
  
  decadal_trends_past_futssp585$scenario <- factor(decadal_trends_past_futssp585$scenario, 
                                                   levels = c("Factual prediction", "Counterfactual climate", "Counterfactual land use", 
                                                              "Counterfactual climate + land use"))
  
  # Subset the data frames to only contain prediction results from the
  # decades of the 1970s, 1990s, 2010s, 2030s, and 2050s
  decadal_trends_past_futssp126 <- decadal_trends_past_futssp126[decadal_trends_past_futssp126$decade %in% c("1970s", "1990s", "2010s", "2030s", "2050s"), ]
  decadal_trends_past_futssp370 <- decadal_trends_past_futssp370[decadal_trends_past_futssp370$decade %in% c("1970s", "1990s", "2010s", "2030s", "2050s"), ]
  decadal_trends_past_futssp585 <- decadal_trends_past_futssp585[decadal_trends_past_futssp585$decade %in% c("1970s", "1990s", "2010s", "2030s", "2050s"), ]
  

  
    
# b) Visualisation of data for scenario ssp126 ---------------------------------
  
  print("Visualise ssp126")
  
  # Visualise the data for future environmental scenario ssp126
  # for Ixodes ricinus and TBE
  ggplot(data = decadal_trends_past_futssp126, 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
      strip = strip_themed(background_y = list(
          "Ixodes ricinus" = element_rect(fill = "steelblue3"),
          "TBE" = element_rect(fill = "steelblue3"),
          "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
          "WNV" = element_rect(fill = "lightsteelblue1")       
        ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        "1990s" = "royalblue4",
        #"2000s" = "seagreen",
        "2010s" = "lightblue3",
        "2030s" = "lightcoral",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid",
        "Counterfactual climate" = "dotdash",
        "Counterfactual land use" = "dashed",
        "Counterfactual climate + land use" = "dotted"
      )
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 14.5),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12.5, face = "bold"),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 15, face = "bold"),
      strip.text = element_text(size = 15, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(1.5, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
      linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
    )
  
  # Save the figure
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Ixodes_TBE_Culex_WNV_ssp126_",o,".png"), width = 9, height = 11)
  
  
  
  
  
# c) Visualisation of data for scenario ssp370 ---------------------------------
  
  print("Visualise ssp370")
  
  # Visualize the data for future environmental scenario ssp370
  # for Ixodes ricinus and TBE
  ggplot(data = decadal_trends_past_futssp370,
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3"),
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1")  
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        "1990s" = "royalblue4",
        #"2000s" = "seagreen",
        "2010s" = "lightblue3",
        "2030s" = "lightcoral",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid",
        "Counterfactual climate" = "dotdash",
        "Counterfactual land use" = "dashed",
        "Counterfactual climate + land use" = "dotted"
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
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(1.5, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
      linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
    ) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
  
  
# Save the figure
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Ixodes_TBE_Culex_WNV_ssp370_",o,".png"), width = 9, height = 11)
  
  
  
  
# d) Visualisation of data for scenario ssp585 ---------------------------------
  
  print("Visualise ssp585")
  
  # Visualize the data for future environmental scenario ssp585
  # for Ixodes ricinus and TBE
  ggplot(data = decadal_trends_past_futssp585, 
         aes(x = month, y = occurrence, color = decade, linetype = scenario, group = interaction(decade, scenario))) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
                strip = strip_themed(background_y = list(
                  "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                  "TBE" = element_rect(fill = "steelblue3"),
                  "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
                  "WNV" = element_rect(fill = "lightsteelblue1") 
                ))) +
    labs(x = "Month in a year", y = paste(o, "occurrence probability"), color = "Decade", linetype = "Prediction basis") +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        #"1980s" = "midnightblue",
        "1990s" = "royalblue4",
        #"2000s" = "seagreen",
        "2010s" = "lightblue3",
        "2030s" = "lightcoral",
        #"2040s" = "#CF4446",
        "2050s" = "#A52C60"
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Factual prediction" = "solid",
        "Counterfactual climate" = "dotdash",
        "Counterfactual land use" = "dashed",
        "Counterfactual climate + land use" = "dotted"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 14.5),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 12.5, face = "bold"),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 15, face = "bold"),
      strip.text = element_text(size = 15, face = "bold"), 
      strip.background = element_rect(fill = "grey75", color = NA),
      panel.grid.major = element_line(linewidth = 0.3, color = "gray90"),
      panel.grid.minor = element_blank(), 
      axis.ticks.length = unit(0.3, "cm"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width = unit(1.5, "cm")) +
    guides(
      color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
      linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
    )
  
  # Save the figure
  ggsave(paste0("output_data/plots/decadal_trends/decadal_trends_Ixodes_TBE_Culex_WNV_ssp585_",o,".png"), width = 9, height = 11)
  
  
  
} # Close the loop over the peak and mean functions





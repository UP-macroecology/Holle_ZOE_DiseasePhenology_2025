# ZOE disease phenology analysis


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                          15. Duration trends                           #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We calculate and visualise temporal trends in the duration of vector activity 
# and potential virus transmission periods throughout the year to examine how
# these patterns may be changing across time and space. Using the ensemble
# predictions based on the factual historical climate and land-use changes, 
# as well as on the three future socio-economic scenarios, we first determine 
# the number of months with predicted presence in each cell. Next, we compute the
# mean duration per decade for each cell. To assess temporal historical and 
# future trends per cell, we fit two separate linear models for each cell: one
# based on the decadal mean durations of the historical decades (1970s - 2010s), 
# and one based on those of the future decades (2010s - 2050s). The slopes of 
# these models serve as an indicator of change. 


# Load the needed packages
library(terra) # terra_1.7-55
library(tidyterra) # tidyterra_0.6.1
library(tidyverse) # tidyverse_2.0.0
library(ggplot2) # ggplot2_4.0.0
library(ggh4x) # ggh4x_0.3.1 
library(sf) # sf_1.0-16
library(ggnewscale) # ggnewscale_0.5.1
library(viridis) # viridis_0.6.4

# Load needed data
europe_mask <- terra::rast("input_data/spatial_data/europe_mask.tif")


#-------------------------------------------------------------------------------

# 1. Historical duration trends ------------------------------------------------


# a) Load data -----------------------------------------------------------------

# Load needed data
# Read in the post-processed monthly binary prediction data from 1970 to 2019
# under factual climate and land use
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_1970_2019.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif")

# TBE
TBE_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_1970_2019.tif")

# WNV
WNV_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_1970_2019.tif")



# b) Calculation of trend values -----------------------------------------------

# For the Spatraster containing the predictions of the WNV, remove all layers with
# predictions made for the months January to April, as these predictions are too
# unreliable given that there were no disease ocpastences reported
WNV_occ_prob_clim_landuse_ens_bin_names <- names(WNV_occ_prob_clim_landuse_ens_bin) # Extract layer names
WNV_occ_prob_clim_landuse_ens_bin_months <- substr(WNV_occ_prob_clim_landuse_ens_bin_names, 1, 2)
WNV_occ_prob_clim_landuse_ens_bin_keep <- which(!(WNV_occ_prob_clim_landuse_ens_bin_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_landuse_ens_bin <- WNV_occ_prob_clim_landuse_ens_bin[[WNV_occ_prob_clim_landuse_ens_bin_keep]]


# Extract the years from the layer names for the vectors and viruses
I_ricinus_occ_prob_clim_landuse_ens_bin_names <- names(I_ricinus_occ_prob_clim_landuse_ens_bin)
I_ricinus_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(I_ricinus_occ_prob_clim_landuse_ens_bin_names, 4, 7))

C_pipiens_occ_prob_clim_landuse_ens_bin_names <- names(C_pipiens_occ_prob_clim_landuse_ens_bin)
C_pipiens_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(C_pipiens_occ_prob_clim_landuse_ens_bin_names, 4, 7))

TBE_occ_prob_clim_landuse_ens_bin_names <- names(TBE_occ_prob_clim_landuse_ens_bin)
TBE_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(TBE_occ_prob_clim_landuse_ens_bin_names, 4, 7))

WNV_occ_prob_clim_landuse_ens_bin_names <- names(WNV_occ_prob_clim_landuse_ens_bin)
WNV_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(WNV_occ_prob_clim_landuse_ens_bin_names, 4, 7))


# Calculate the number of months with a predicted presence per year (for each cell)
I_ricinus_length_presence_year <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin, index = I_ricinus_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin, index = C_pipiens_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)
TBE_length_presence_year <- terra::tapp(TBE_occ_prob_clim_landuse_ens_bin, index = TBE_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)
WNV_length_presence_year <- terra::tapp(WNV_occ_prob_clim_landuse_ens_bin, index = WNV_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)

# Rename yearly layers
names(I_ricinus_length_presence_year) <- paste0(unique(I_ricinus_occ_prob_clim_landuse_ens_bin_years))
names(C_pipiens_length_presence_year) <- paste0(unique(C_pipiens_occ_prob_clim_landuse_ens_bin_years))
names(TBE_length_presence_year) <- paste0(unique(TBE_occ_prob_clim_landuse_ens_bin_years))
names(WNV_length_presence_year) <- paste0(unique(WNV_occ_prob_clim_landuse_ens_bin_years))

# Create a vector with the studied decades
I_ricinus_r_names_past <- names(I_ricinus_length_presence_year)
I_ricinus_r_decades_past <- floor(as.numeric(I_ricinus_r_names_past) / 10) * 10

C_pipiens_r_names_past <- names(C_pipiens_length_presence_year)
C_pipiens_r_decades_past <- floor(as.numeric(C_pipiens_r_names_past) / 10) * 10

TBE_r_names_past <- names(TBE_length_presence_year)
TBE_r_decades_past <- floor(as.numeric(TBE_r_names_past) / 10) * 10

WNV_r_names_past <- names(WNV_length_presence_year)
WNV_r_decades_past <- floor(as.numeric(WNV_r_names_past) / 10) * 10


# Compute the mean presence length per decade (per cell)
I_ricinus_decadal_means_past <- tapp(I_ricinus_length_presence_year, index = I_ricinus_r_decades_past, fun = mean, na.rm = TRUE)
C_pipiens_decadal_means_past <- tapp(C_pipiens_length_presence_year, index = C_pipiens_r_decades_past, fun = mean, na.rm = TRUE)
TBE_decadal_means_past <- tapp(TBE_length_presence_year, index = TBE_r_decades_past, fun = mean, na.rm = TRUE)
WNV_decadal_means_past <- tapp(WNV_length_presence_year, index = WNV_r_decades_past, fun = mean, na.rm = TRUE)

# Rename decadal layers
names(I_ricinus_decadal_means_past) <- paste0(unique(I_ricinus_r_decades_past))
names(C_pipiens_decadal_means_past) <- paste0(unique(C_pipiens_r_decades_past))
names(TBE_decadal_means_past) <- paste0(unique(TBE_r_decades_past))
names(WNV_decadal_means_past) <- paste0(unique(WNV_r_decades_past))

# Extract the rasters containing the mean presence length per cell in the 
# decade of the 2010s for the vectors and viruses - for later usage when
# calculating future duration trend
I_ricinus_decadal_means_past_2010 <- I_ricinus_decadal_means_past[["2010"]]
C_pipiens_decadal_means_past_2010 <- C_pipiens_decadal_means_past[["2010"]]
TBE_decadal_means_past_2010 <- TBE_decadal_means_past[["2010"]]
WNV_decadal_means_past_2010 <- WNV_decadal_means_past[["2010"]]

# Function to compute slope (trend) from 1970s to 2010s for each cell
calc_slope_past <- function(x) {
  if (all(is.na(x))) return(NA) 
  decade_values <- seq(1970, 2010, by = 10)
  lm_fit <- lm(x ~ decade_values)  # Fit linear model
  return(coef(lm_fit)[2])  # Extract slope coefficient
}

# Apply function to raster to calculate trend (slope)
I_ricinus_slope_raster_past <- app(I_ricinus_decadal_means_past , fun = calc_slope_past)
names(I_ricinus_slope_raster_past) <- "Trend_Slope"

C_pipiens_slope_raster_past <- app(C_pipiens_decadal_means_past , fun = calc_slope_past)
names(C_pipiens_slope_raster_past) <- "Trend_Slope"

TBE_slope_raster_past <- app(TBE_decadal_means_past , fun = calc_slope_past)
names(TBE_slope_raster_past) <- "Trend_Slope"

WNV_slope_raster_past <- app(WNV_decadal_means_past , fun = calc_slope_past)
names(WNV_slope_raster_past) <- "Trend_Slope"



# c) Prepare data frame for plotting -------------------------------------------

# Convert raster to a data frame for plotting and add a column indicating the species/pathogen name
I_ricinus_slope_df_past <- as.data.frame(I_ricinus_slope_raster_past, xy = TRUE, na.rm = TRUE)
colnames(I_ricinus_slope_df_past) <- c("lon", "lat", "trend")
I_ricinus_slope_df_past$species <- "Ixodes ricinus"
I_ricinus_slope_df_past$time <- "1970s - 2010s"

C_pipiens_slope_df_past <- as.data.frame(C_pipiens_slope_raster_past, xy = TRUE, na.rm = TRUE)
colnames(C_pipiens_slope_df_past) <- c("lon", "lat", "trend")
C_pipiens_slope_df_past$species <- "Culex pipiens"
C_pipiens_slope_df_past$time <- "1970s - 2010s"

TBE_slope_df_past <- as.data.frame(TBE_slope_raster_past, xy = TRUE, na.rm = TRUE)
colnames(TBE_slope_df_past) <- c("lon", "lat", "trend")
TBE_slope_df_past$species <- "TBE"
TBE_slope_df_past$time <- "1970s - 2010s"

WNV_slope_df_past <- as.data.frame(WNV_slope_raster_past, xy = TRUE, na.rm = TRUE)
colnames(WNV_slope_df_past) <- c("lon", "lat", "trend")
WNV_slope_df_past$species <- "WNV"
WNV_slope_df_past$time <- "1970s - 2010s"

# Bind the four different data frames
slope_df_past <- rbind(I_ricinus_slope_df_past, C_pipiens_slope_df_past,
                       TBE_slope_df_past, WNV_slope_df_past)




# d) Save resulting data frame --------------------------------------------------

save(slope_df_past, file = "output_data/results/duration_trends/duration_trends_past_vector_virus.RData")




#-------------------------------------------------------------------------------

# 2. Future duration trends ----------------------------------------------------


# a) Load data -----------------------------------------------------------------

# Load needed data - post-processed monthly predicted ensemble occurrence 
# probabilities of main vectors and the respective viruses, under climate and 
# land-use change scenarios for the future years 2020 to 2059 based on three 
# studied environmental scenarios
# Ixodes ricinus; ssp126, ssp370, and ssp585
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif"))
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif"))

# Culex pipiens; ssp126, ssp370, and ssp585
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif"))

# TBE; ssp126, ssp370, and ssp585
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif"))
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif"))
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif"))

# WNV; ssp126, ssp370, and ssp585
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif"))
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif"))
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif"))



# b) Calculation of trend values -----------------------------------------------

# For the Spatraster containing the predictions of the WNV, remove all layers with
# predictions made for the months January to April, as these predictions are too
# unreliable given that there were no disease occurrences reported
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_names <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126) # Extract layer names
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_months <- substr(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_names, 1, 2)
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_keep <- which(!(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126[[WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126_keep]]

WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_names <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370) # Extract layer names
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_months <- substr(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_names, 1, 2)
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_keep <- which(!(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370[[WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370_keep]]

WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_names <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585) # Extract layer names
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_months <- substr(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_names, 1, 2)
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_keep <- which(!(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126[[WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585_keep]]


# Extract the years from the layer names for the vectors and viruses
I_ricinus_r_names_fut_ssp126 <- names(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126)
I_ricinus_r_years_fut_ssp126 <- as.numeric(substr(I_ricinus_r_names_fut_ssp126, 4, 7))
I_ricinus_r_names_fut_ssp370 <- names(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370)
I_ricinus_r_years_fut_ssp370 <- as.numeric(substr(I_ricinus_r_names_fut_ssp370, 4, 7))
I_ricinus_r_names_fut_ssp585 <- names(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585)
I_ricinus_r_years_fut_ssp585 <- as.numeric(substr(I_ricinus_r_names_fut_ssp585, 4, 7))

C_pipiens_r_names_fut_ssp126 <- names(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126)
C_pipiens_r_years_fut_ssp126 <- as.numeric(substr(C_pipiens_r_names_fut_ssp126, 4, 7))
C_pipiens_r_names_fut_ssp370 <- names(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370)
C_pipiens_r_years_fut_ssp370 <- as.numeric(substr(C_pipiens_r_names_fut_ssp370, 4, 7))
C_pipiens_r_names_fut_ssp585 <- names(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585)
C_pipiens_r_years_fut_ssp585 <- as.numeric(substr(C_pipiens_r_names_fut_ssp585, 4, 7))

TBE_r_names_fut_ssp126 <- names(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126)
TBE_r_years_fut_ssp126 <- as.numeric(substr(TBE_r_names_fut_ssp126, 4, 7))
TBE_r_names_fut_ssp370 <- names(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370)
TBE_r_years_fut_ssp370 <- as.numeric(substr(TBE_r_names_fut_ssp370, 4, 7))
TBE_r_names_fut_ssp585 <- names(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585)
TBE_r_years_fut_ssp585 <- as.numeric(substr(TBE_r_names_fut_ssp585, 4, 7))

WNV_r_names_fut_ssp126 <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126)
WNV_r_years_fut_ssp126 <- as.numeric(substr(WNV_r_names_fut_ssp126, 4, 7))
WNV_r_names_fut_ssp370 <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370)
WNV_r_years_fut_ssp370 <- as.numeric(substr(WNV_r_names_fut_ssp370, 4, 7))
WNV_r_names_fut_ssp585 <- names(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585)
WNV_r_years_fut_ssp585 <- as.numeric(substr(WNV_r_names_fut_ssp585, 4, 7))

# Calculate the number of months with a predicted presence per year (for each cell)
I_ricinus_length_presence_year_fut_ssp126 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = I_ricinus_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
I_ricinus_length_presence_year_fut_ssp370 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = I_ricinus_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
I_ricinus_length_presence_year_fut_ssp585 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = I_ricinus_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

C_pipiens_length_presence_year_fut_ssp126 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = C_pipiens_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year_fut_ssp370 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = C_pipiens_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year_fut_ssp585 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = C_pipiens_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

TBE_length_presence_year_fut_ssp126 <- terra::tapp(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = TBE_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
TBE_length_presence_year_fut_ssp370 <- terra::tapp(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = TBE_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
TBE_length_presence_year_fut_ssp585 <- terra::tapp(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = TBE_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

WNV_length_presence_year_fut_ssp126 <- terra::tapp(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = WNV_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
WNV_length_presence_year_fut_ssp370 <- terra::tapp(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = WNV_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
WNV_length_presence_year_fut_ssp585 <- terra::tapp(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = WNV_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

# Rename yearly layers
names(I_ricinus_length_presence_year_fut_ssp126) <- paste0(unique(I_ricinus_r_years_fut_ssp126))
names(I_ricinus_length_presence_year_fut_ssp370) <- paste0(unique(I_ricinus_r_years_fut_ssp370))
names(I_ricinus_length_presence_year_fut_ssp585) <- paste0(unique(I_ricinus_r_years_fut_ssp585))

names(C_pipiens_length_presence_year_fut_ssp126) <- paste0(unique(C_pipiens_r_years_fut_ssp126))
names(C_pipiens_length_presence_year_fut_ssp370) <- paste0(unique(C_pipiens_r_years_fut_ssp370))
names(C_pipiens_length_presence_year_fut_ssp585) <- paste0(unique(C_pipiens_r_years_fut_ssp585))

names(TBE_length_presence_year_fut_ssp126) <- paste0(unique(TBE_r_years_fut_ssp126))
names(TBE_length_presence_year_fut_ssp370) <- paste0(unique(TBE_r_years_fut_ssp370))
names(TBE_length_presence_year_fut_ssp585) <- paste0(unique(TBE_r_years_fut_ssp585))

names(WNV_length_presence_year_fut_ssp126) <- paste0(unique(WNV_r_years_fut_ssp126))
names(WNV_length_presence_year_fut_ssp370) <- paste0(unique(WNV_r_years_fut_ssp370))
names(WNV_length_presence_year_fut_ssp585) <- paste0(unique(WNV_r_years_fut_ssp585))

# Create a vector with the studied decades
I_ricinus_r_names_fut_ssp126 <- names(I_ricinus_length_presence_year_fut_ssp126)
I_ricinus_r_decades_fut_ssp126 <- floor(as.numeric(I_ricinus_r_names_fut_ssp126) / 10) * 10
I_ricinus_r_names_fut_ssp370 <- names(I_ricinus_length_presence_year_fut_ssp370)
I_ricinus_r_decades_fut_ssp370 <- floor(as.numeric(I_ricinus_r_names_fut_ssp370) / 10) * 10
I_ricinus_r_names_fut_ssp585 <- names(I_ricinus_length_presence_year_fut_ssp585)
I_ricinus_r_decades_fut_ssp585 <- floor(as.numeric(I_ricinus_r_names_fut_ssp585) / 10) * 10

C_pipiens_r_names_fut_ssp126 <- names(C_pipiens_length_presence_year_fut_ssp126)
C_pipiens_r_decades_fut_ssp126 <- floor(as.numeric(C_pipiens_r_names_fut_ssp126) / 10) * 10
C_pipiens_r_names_fut_ssp370 <- names(C_pipiens_length_presence_year_fut_ssp370)
C_pipiens_r_decades_fut_ssp370 <- floor(as.numeric(C_pipiens_r_names_fut_ssp370) / 10) * 10
C_pipiens_r_names_fut_ssp585 <- names(C_pipiens_length_presence_year_fut_ssp585)
C_pipiens_r_decades_fut_ssp585 <- floor(as.numeric(C_pipiens_r_names_fut_ssp585) / 10) * 10

TBE_r_names_fut_ssp126 <- names(TBE_length_presence_year_fut_ssp126)
TBE_r_decades_fut_ssp126 <- floor(as.numeric(TBE_r_names_fut_ssp126) / 10) * 10
TBE_r_names_fut_ssp370 <- names(TBE_length_presence_year_fut_ssp370)
TBE_r_decades_fut_ssp370 <- floor(as.numeric(TBE_r_names_fut_ssp370) / 10) * 10
TBE_r_names_fut_ssp585 <- names(TBE_length_presence_year_fut_ssp585)
TBE_r_decades_fut_ssp585 <- floor(as.numeric(TBE_r_names_fut_ssp585) / 10) * 10

WNV_r_names_fut_ssp126 <- names(WNV_length_presence_year_fut_ssp126)
WNV_r_decades_fut_ssp126 <- floor(as.numeric(WNV_r_names_fut_ssp126) / 10) * 10
WNV_r_names_fut_ssp370 <- names(WNV_length_presence_year_fut_ssp370)
WNV_r_decades_fut_ssp370 <- floor(as.numeric(WNV_r_names_fut_ssp370) / 10) * 10
WNV_r_names_fut_ssp585 <- names(WNV_length_presence_year_fut_ssp585)
WNV_r_decades_fut_ssp585 <- floor(as.numeric(WNV_r_names_fut_ssp585) / 10) * 10



# Compute the mean presence length per decade (per cell)
I_ricinus_decadal_means_fut_ssp126 <- tapp(I_ricinus_length_presence_year_fut_ssp126, index = I_ricinus_r_decades_fut_ssp126, fun = mean, na.rm = TRUE)
I_ricinus_decadal_means_fut_ssp370 <- tapp(I_ricinus_length_presence_year_fut_ssp370, index = I_ricinus_r_decades_fut_ssp370, fun = mean, na.rm = TRUE)
I_ricinus_decadal_means_fut_ssp585 <- tapp(I_ricinus_length_presence_year_fut_ssp585, index = I_ricinus_r_decades_fut_ssp585, fun = mean, na.rm = TRUE)

C_pipiens_decadal_means_fut_ssp126 <- tapp(C_pipiens_length_presence_year_fut_ssp126, index = C_pipiens_r_decades_fut_ssp126, fun = mean, na.rm = TRUE)
C_pipiens_decadal_means_fut_ssp370 <- tapp(C_pipiens_length_presence_year_fut_ssp370, index = C_pipiens_r_decades_fut_ssp370, fun = mean, na.rm = TRUE)
C_pipiens_decadal_means_fut_ssp585 <- tapp(C_pipiens_length_presence_year_fut_ssp585, index = C_pipiens_r_decades_fut_ssp585, fun = mean, na.rm = TRUE)

TBE_decadal_means_fut_ssp126 <- tapp(TBE_length_presence_year_fut_ssp126, index = TBE_r_decades_fut_ssp126, fun = mean, na.rm = TRUE)
TBE_decadal_means_fut_ssp370 <- tapp(TBE_length_presence_year_fut_ssp370, index = TBE_r_decades_fut_ssp370, fun = mean, na.rm = TRUE)
TBE_decadal_means_fut_ssp585 <- tapp(TBE_length_presence_year_fut_ssp585, index = TBE_r_decades_fut_ssp585, fun = mean, na.rm = TRUE)

WNV_decadal_means_fut_ssp126 <- tapp(WNV_length_presence_year_fut_ssp126, index = WNV_r_decades_fut_ssp126, fun = mean, na.rm = TRUE)
WNV_decadal_means_fut_ssp370 <- tapp(WNV_length_presence_year_fut_ssp370, index = WNV_r_decades_fut_ssp370, fun = mean, na.rm = TRUE)
WNV_decadal_means_fut_ssp585 <- tapp(WNV_length_presence_year_fut_ssp585, index = WNV_r_decades_fut_ssp585, fun = mean, na.rm = TRUE)


# Rename decadal layers
names(I_ricinus_decadal_means_fut_ssp126) <- paste0(unique(I_ricinus_r_decades_fut_ssp126))
names(I_ricinus_decadal_means_fut_ssp370) <- paste0(unique(I_ricinus_r_decades_fut_ssp370))
names(I_ricinus_decadal_means_fut_ssp585) <- paste0(unique(I_ricinus_r_decades_fut_ssp585))

names(C_pipiens_decadal_means_fut_ssp126) <- paste0(unique(C_pipiens_r_decades_fut_ssp126))
names(C_pipiens_decadal_means_fut_ssp370) <- paste0(unique(C_pipiens_r_decades_fut_ssp370))
names(C_pipiens_decadal_means_fut_ssp585) <- paste0(unique(C_pipiens_r_decades_fut_ssp585))

names(TBE_decadal_means_fut_ssp126) <- paste0(unique(TBE_r_decades_fut_ssp126))
names(TBE_decadal_means_fut_ssp370) <- paste0(unique(TBE_r_decades_fut_ssp370))
names(TBE_decadal_means_fut_ssp585) <- paste0(unique(TBE_r_decades_fut_ssp585))

names(WNV_decadal_means_fut_ssp126) <- paste0(unique(WNV_r_decades_fut_ssp126))
names(WNV_decadal_means_fut_ssp370) <- paste0(unique(WNV_r_decades_fut_ssp370))
names(WNV_decadal_means_fut_ssp585) <- paste0(unique(WNV_r_decades_fut_ssp585))

# Add the mean presence length per cell of the decades of the 2010s as first layer
# to include these into future trend calculations
I_ricinus_decadal_means_fut_ssp126 <- c(I_ricinus_decadal_means_past_2010, I_ricinus_decadal_means_fut_ssp126)
I_ricinus_decadal_means_fut_ssp370 <- c(I_ricinus_decadal_means_past_2010, I_ricinus_decadal_means_fut_ssp370)
I_ricinus_decadal_means_fut_ssp585 <- c(I_ricinus_decadal_means_past_2010, I_ricinus_decadal_means_fut_ssp585)

C_pipiens_decadal_means_fut_ssp126 <- c(C_pipiens_decadal_means_past_2010, C_pipiens_decadal_means_fut_ssp126)
C_pipiens_decadal_means_fut_ssp370 <- c(C_pipiens_decadal_means_past_2010, C_pipiens_decadal_means_fut_ssp370)
C_pipiens_decadal_means_fut_ssp585 <- c(C_pipiens_decadal_means_past_2010, C_pipiens_decadal_means_fut_ssp585)

TBE_decadal_means_fut_ssp126 <- c(TBE_decadal_means_past_2010, TBE_decadal_means_fut_ssp126)
TBE_decadal_means_fut_ssp370 <- c(TBE_decadal_means_past_2010, TBE_decadal_means_fut_ssp370)
TBE_decadal_means_fut_ssp585 <- c(TBE_decadal_means_past_2010, TBE_decadal_means_fut_ssp585)

WNV_decadal_means_fut_ssp126 <- c(WNV_decadal_means_past_2010, WNV_decadal_means_fut_ssp126)
WNV_decadal_means_fut_ssp370 <- c(WNV_decadal_means_past_2010, WNV_decadal_means_fut_ssp370)
WNV_decadal_means_fut_ssp585 <- c(WNV_decadal_means_past_2010, WNV_decadal_means_fut_ssp585)


# Function to compute slope (trend) from 1910s to 2050s for each cell
calc_slope_fut <- function(x) {
  if (all(is.na(x))) return(NA)
  decade_values <- c(2010, 2020, 2030, 2040, 2050)
  lm_fit <- lm(x ~ decade_values)  # Fit linear model
  return(coef(lm_fit)[2])  # Extract slope coefficient
}

# Filter the Spatrasters to only contain decades from 2010s to 2050s
I_ricinus_decadal_means_fut_ssp126 <- subset(I_ricinus_decadal_means_fut_ssp126, names(I_ricinus_decadal_means_fut_ssp126) %in% c("2010", "2020", "2030", "2040", "2050"))
I_ricinus_decadal_means_fut_ssp370 <- subset(I_ricinus_decadal_means_fut_ssp370, names(I_ricinus_decadal_means_fut_ssp370) %in% c("2010", "2020", "2030", "2040", "2050"))
I_ricinus_decadal_means_fut_ssp585 <- subset(I_ricinus_decadal_means_fut_ssp585, names(I_ricinus_decadal_means_fut_ssp585) %in% c("2010", "2020", "2030", "2040", "2050"))

C_pipiens_decadal_means_fut_ssp126 <- subset(C_pipiens_decadal_means_fut_ssp126, names(C_pipiens_decadal_means_fut_ssp126) %in% c("2010", "2020", "2030", "2040", "2050"))
C_pipiens_decadal_means_fut_ssp370 <- subset(C_pipiens_decadal_means_fut_ssp370, names(C_pipiens_decadal_means_fut_ssp370) %in% c("2010", "2020", "2030", "2040", "2050"))
C_pipiens_decadal_means_fut_ssp585 <- subset(C_pipiens_decadal_means_fut_ssp585, names(C_pipiens_decadal_means_fut_ssp585) %in% c("2010", "2020", "2030", "2040", "2050"))

TBE_decadal_means_fut_ssp126 <- subset(TBE_decadal_means_fut_ssp126, names(TBE_decadal_means_fut_ssp126) %in% c("2010", "2020", "2030", "2040", "2050"))
TBE_decadal_means_fut_ssp370 <- subset(TBE_decadal_means_fut_ssp370, names(TBE_decadal_means_fut_ssp370) %in% c("2010", "2020", "2030", "2040", "2050"))
TBE_decadal_means_fut_ssp585 <- subset(TBE_decadal_means_fut_ssp585, names(TBE_decadal_means_fut_ssp585) %in% c("2010", "2020", "2030", "2040", "2050"))

WNV_decadal_means_fut_ssp126 <- subset(WNV_decadal_means_fut_ssp126, names(WNV_decadal_means_fut_ssp126) %in% c("2010", "2020", "2030", "2040", "2050"))
WNV_decadal_means_fut_ssp370 <- subset(WNV_decadal_means_fut_ssp370, names(WNV_decadal_means_fut_ssp370) %in% c("2010", "2020", "2030", "2040", "2050"))
WNV_decadal_means_fut_ssp585 <- subset(WNV_decadal_means_fut_ssp585, names(WNV_decadal_means_fut_ssp585) %in% c("2010", "2020", "2030", "2040", "2050"))

# Apply function to raster to calculate trend (slope)
I_ricinus_slope_raster_fut_ssp126 <- app(I_ricinus_decadal_means_fut_ssp126, fun = calc_slope_fut)
names(I_ricinus_slope_raster_fut_ssp126) <- "Trend_Slope"

I_ricinus_slope_raster_fut_ssp370 <- app(I_ricinus_decadal_means_fut_ssp370, fun = calc_slope_fut)
names(I_ricinus_slope_raster_fut_ssp370) <- "Trend_Slope"

I_ricinus_slope_raster_fut_ssp585 <- app(I_ricinus_decadal_means_fut_ssp585, fun = calc_slope_fut)
names(I_ricinus_slope_raster_fut_ssp585) <- "Trend_Slope"

C_pipiens_slope_raster_fut_ssp126 <- app(C_pipiens_decadal_means_fut_ssp126, fun = calc_slope_fut)
names(C_pipiens_slope_raster_fut_ssp126) <- "Trend_Slope"

C_pipiens_slope_raster_fut_ssp370 <- app(C_pipiens_decadal_means_fut_ssp370, fun = calc_slope_fut)
names(C_pipiens_slope_raster_fut_ssp370) <- "Trend_Slope"

C_pipiens_slope_raster_fut_ssp585 <- app(C_pipiens_decadal_means_fut_ssp585, fun = calc_slope_fut)
names(C_pipiens_slope_raster_fut_ssp585) <- "Trend_Slope"

TBE_slope_raster_fut_ssp126 <- app(TBE_decadal_means_fut_ssp126, fun = calc_slope_fut)
names(TBE_slope_raster_fut_ssp126) <- "Trend_Slope"

TBE_slope_raster_fut_ssp370 <- app(TBE_decadal_means_fut_ssp370, fun = calc_slope_fut)
names(TBE_slope_raster_fut_ssp370) <- "Trend_Slope"

TBE_slope_raster_fut_ssp585 <- app(TBE_decadal_means_fut_ssp585, fun = calc_slope_fut)
names(TBE_slope_raster_fut_ssp585) <- "Trend_Slope"

WNV_slope_raster_fut_ssp126 <- app(WNV_decadal_means_fut_ssp126, fun = calc_slope_fut)
names(WNV_slope_raster_fut_ssp126) <- "Trend_Slope"

WNV_slope_raster_fut_ssp370 <- app(WNV_decadal_means_fut_ssp370, fun = calc_slope_fut)
names(WNV_slope_raster_fut_ssp370) <- "Trend_Slope"

WNV_slope_raster_fut_ssp585 <- app(WNV_decadal_means_fut_ssp585, fun = calc_slope_fut)
names(WNV_slope_raster_fut_ssp585) <- "Trend_Slope"



# c) Prepare data frame for plotting -------------------------------------------

# Convert raster to a data frame for plotting
I_ricinus_slope_df_fut_ssp126 <- as.data.frame(I_ricinus_slope_raster_fut_ssp126, xy = TRUE, na.rm = TRUE)
colnames(I_ricinus_slope_df_fut_ssp126) <- c("lon", "lat", "trend")
I_ricinus_slope_df_fut_ssp126$species <- "Ixodes ricinus"
I_ricinus_slope_df_fut_ssp126$time <- "2010s - 2050s; ssp126"

I_ricinus_slope_df_fut_ssp370 <- as.data.frame(I_ricinus_slope_raster_fut_ssp370, xy = TRUE, na.rm = TRUE)
colnames(I_ricinus_slope_df_fut_ssp370) <- c("lon", "lat", "trend")
I_ricinus_slope_df_fut_ssp370$species <- "Ixodes ricinus"
I_ricinus_slope_df_fut_ssp370$time <- "2010s - 2050s; ssp370"

I_ricinus_slope_df_fut_ssp585 <- as.data.frame(I_ricinus_slope_raster_fut_ssp585, xy = TRUE, na.rm = TRUE)
colnames(I_ricinus_slope_df_fut_ssp585) <- c("lon", "lat", "trend")
I_ricinus_slope_df_fut_ssp585$species <- "Ixodes ricinus"
I_ricinus_slope_df_fut_ssp585$time <- "2010s - 2050s; ssp585"

C_pipiens_slope_df_fut_ssp126 <- as.data.frame(C_pipiens_slope_raster_fut_ssp126, xy = TRUE, na.rm = TRUE)
colnames(C_pipiens_slope_df_fut_ssp126) <- c("lon", "lat", "trend")
C_pipiens_slope_df_fut_ssp126$species <- "Culex pipiens"
C_pipiens_slope_df_fut_ssp126$time <- "2010s - 2050s; ssp126"

C_pipiens_slope_df_fut_ssp370 <- as.data.frame(C_pipiens_slope_raster_fut_ssp370, xy = TRUE, na.rm = TRUE)
colnames(C_pipiens_slope_df_fut_ssp370) <- c("lon", "lat", "trend")
C_pipiens_slope_df_fut_ssp370$species <- "Culex pipiens"
C_pipiens_slope_df_fut_ssp370$time <- "2010s - 2050s; ssp370"

C_pipiens_slope_df_fut_ssp585 <- as.data.frame(C_pipiens_slope_raster_fut_ssp585, xy = TRUE, na.rm = TRUE)
colnames(C_pipiens_slope_df_fut_ssp585) <- c("lon", "lat", "trend")
C_pipiens_slope_df_fut_ssp585$species <- "Culex pipiens"
C_pipiens_slope_df_fut_ssp585$time <- "2010s - 2050s; ssp585"

TBE_slope_df_fut_ssp126 <- as.data.frame(TBE_slope_raster_fut_ssp126, xy = TRUE, na.rm = TRUE)
colnames(TBE_slope_df_fut_ssp126) <- c("lon", "lat", "trend")
TBE_slope_df_fut_ssp126$species <- "TBE"
TBE_slope_df_fut_ssp126$time <- "2010s - 2050s; ssp126"

TBE_slope_df_fut_ssp370 <- as.data.frame(TBE_slope_raster_fut_ssp370, xy = TRUE, na.rm = TRUE)
colnames(TBE_slope_df_fut_ssp370) <- c("lon", "lat", "trend")
TBE_slope_df_fut_ssp370$species <- "TBE"
TBE_slope_df_fut_ssp370$time <- "2010s - 2050s; ssp370"

TBE_slope_df_fut_ssp585 <- as.data.frame(TBE_slope_raster_fut_ssp585, xy = TRUE, na.rm = TRUE)
colnames(TBE_slope_df_fut_ssp585) <- c("lon", "lat", "trend")
TBE_slope_df_fut_ssp585$species <- "TBE"
TBE_slope_df_fut_ssp585$time <- "2010s - 2050s; ssp585"

WNV_slope_df_fut_ssp126 <- as.data.frame(WNV_slope_raster_fut_ssp126, xy = TRUE, na.rm = TRUE)
colnames(WNV_slope_df_fut_ssp126) <- c("lon", "lat", "trend")
WNV_slope_df_fut_ssp126$species <- "WNV"
WNV_slope_df_fut_ssp126$time <- "2010s - 2050s; ssp126"

WNV_slope_df_fut_ssp370 <- as.data.frame(WNV_slope_raster_fut_ssp370, xy = TRUE, na.rm = TRUE)
colnames(WNV_slope_df_fut_ssp370) <- c("lon", "lat", "trend")
WNV_slope_df_fut_ssp370$species <- "WNV"
WNV_slope_df_fut_ssp370$time <- "2010s - 2050s; ssp370"

WNV_slope_df_fut_ssp585 <- as.data.frame(WNV_slope_raster_fut_ssp585, xy = TRUE, na.rm = TRUE)
colnames(WNV_slope_df_fut_ssp585) <- c("lon", "lat", "trend")
WNV_slope_df_fut_ssp585$species <- "WNV"
WNV_slope_df_fut_ssp585$time <- "2010s - 2050s; ssp585"

# Bind the different data frames for each ssp scenario we look at
slope_df_futssp126 <- rbind(I_ricinus_slope_df_fut_ssp126, C_pipiens_slope_df_fut_ssp126,
                            TBE_slope_df_fut_ssp126, WNV_slope_df_fut_ssp126)

slope_df_futssp370 <- rbind(I_ricinus_slope_df_fut_ssp370, C_pipiens_slope_df_fut_ssp370,
                            TBE_slope_df_fut_ssp370, WNV_slope_df_fut_ssp370)

slope_df_futssp585 <- rbind(I_ricinus_slope_df_fut_ssp585, C_pipiens_slope_df_fut_ssp585,
                            TBE_slope_df_fut_ssp585, WNV_slope_df_fut_ssp585)



# d) Save resulting data frames ------------------------------------------------

save(slope_df_futssp126, file = "output_data/results/duration_trends/duration_trends_fut_vector_virus_ssp126.RData")
save(slope_df_futssp370, file = "output_data/results/duration_trends/duration_trends_fut_vector_virus_ssp370.RData")
save(slope_df_futssp585, file = "output_data/results/duration_trends/duration_trends_fut_vector_virus_ssp585.RData")






#-------------------------------------------------------------------------------

# 3. Visualise duration trends  ------------------------------------------------
# for main vectors as well as viruses
# showing duration of 2010s and future trends (2010s - 2050s)


# a) Prepare resulting past and future trend data frame ------------------------

# Convert raster with mean duration trends in the 2010s into a data frame
# and add important columns
I_ricinus_decadal_means_past_2010_df <- as.data.frame(I_ricinus_decadal_means_past_2010, xy = TRUE)
colnames(I_ricinus_decadal_means_past_2010_df) <- c("lon", "lat", "trend")
I_ricinus_decadal_means_past_2010_df$species <- "Ixodes ricinus"
I_ricinus_decadal_means_past_2010_df$time <- "2010s"
I_ricinus_decadal_means_past_2010_df$legend <- "map1"

C_pipiens_decadal_means_past_2010_df <- as.data.frame(C_pipiens_decadal_means_past_2010, xy = TRUE)
colnames(C_pipiens_decadal_means_past_2010_df) <- c("lon", "lat", "trend")
C_pipiens_decadal_means_past_2010_df$species <- "Culex pipiens"
C_pipiens_decadal_means_past_2010_df$time <- "2010s"
C_pipiens_decadal_means_past_2010_df$legend <- "map1"

TBE_decadal_means_past_2010_df <- as.data.frame(TBE_decadal_means_past_2010, xy = TRUE)
colnames(TBE_decadal_means_past_2010_df) <- c("lon", "lat", "trend")
TBE_decadal_means_past_2010_df$species <- "TBE"
TBE_decadal_means_past_2010_df$time <- "2010s"
TBE_decadal_means_past_2010_df$legend <- "map1"

WNV_decadal_means_past_2010_df <- as.data.frame(WNV_decadal_means_past_2010, xy = TRUE)
colnames(WNV_decadal_means_past_2010_df) <- c("lon", "lat", "trend")
WNV_decadal_means_past_2010_df$species <- "WNV"
WNV_decadal_means_past_2010_df$time <- "2010s"
WNV_decadal_means_past_2010_df$legend <- "map1"

# Bind the four data frames of the different vectors and viruses
decadal_means_2010s <- rbind(I_ricinus_decadal_means_past_2010_df,
                             C_pipiens_decadal_means_past_2010_df,
                             TBE_decadal_means_past_2010_df,
                             WNV_decadal_means_past_2010_df)

# Add a legend column for the data frames containing future predictions
slope_df_futssp126$legend <- "map2"
slope_df_futssp370$legend <- "map2"
slope_df_futssp585$legend <- "map2"

# Bind the data frames
duration_trends_past_futssp126 <- rbind(decadal_means_2010s, slope_df_futssp126)
duration_trends_past_futssp370 <- rbind(decadal_means_2010s, slope_df_futssp370)
duration_trends_past_futssp585 <- rbind(decadal_means_2010s, slope_df_futssp585)


# Make sure the vector and viruses appear in the correct order
duration_trends_past_futssp126$species <- factor(duration_trends_past_futssp126$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))

duration_trends_past_futssp370$species <- factor(duration_trends_past_futssp370$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))

duration_trends_past_futssp585$species <- factor(duration_trends_past_futssp585$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))


# Convert Europe mask spatraster into a data frame
europe_mask_df <- as.data.frame(europe_mask, xy = TRUE)





# b) Visualisation of duration trends for ssp126 -------------------------------

# Calculate the range of the future trend values to symmetrise the make the
# colour legend symmetric around 0
range_values <- range(duration_trends_past_futssp126$trend[duration_trends_past_futssp126$legend == "map2"], na.rm = TRUE)
max_abs <- max(abs(range_values))
  
  # Plot the duration trend per cell based on the factual historical data and the 
  # environmental scenario ssp126
  ggplot() +
  geom_raster(data = europe_mask_df, aes(x = x, y = y), fill = "grey30") +
  geom_raster(data = duration_trends_past_futssp126 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp126 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        limits = c(-max_abs, max_abs), breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(species ~ time, scales = "free_y",
              strip = strip_themed(background_y = elem_list_rect(fill = c("steelblue3", "steelblue3", "lightsteelblue1", "lightsteelblue1")),
                                   text_y = elem_list_text(face = c("bold.italic", NA, "bold.italic", NA)))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
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
    strip.background = element_rect(fill = "grey75", color = NA)
  )

# Save the figure
ggsave("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp126.png", width = 7.5, height = 11)




# c) Visualisation of duration trends for ssp370 -------------------------------

# Calculate the range of the future trend values to symmetrise the make the
# colour legend symmetric around 0
range_values <- range(duration_trends_past_futssp370$trend[duration_trends_past_futssp370$legend == "map2"], na.rm = TRUE)
max_abs <- max(abs(range_values))

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp370
ggplot() +
  geom_raster(data = europe_mask_df, aes(x = x, y = y), fill = "gray30") +
  geom_raster(data = duration_trends_past_futssp370 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp370 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        limits = c(-max_abs, max_abs), breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(species ~ time, scales = "free_y",
              strip = strip_themed(background_y = elem_list_rect(fill = c("steelblue3", "steelblue3", "lightsteelblue1", "lightsteelblue1")),
                                   text_y = elem_list_text(face = c("bold.italic", NA, "bold.italic", NA)))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
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
    strip.background = element_rect(fill = "grey75", color = NA)
  )

# Save the figure
ggsave("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp370.png", width = 7.5, height = 11)




# d) Visualisation of duration trends for ssp585 -------------------------------

# Calculate the range of the future trend values to symmetrise the make the
# colour legend symmetric around 0
range_values <- range(duration_trends_past_futssp585$trend[duration_trends_past_futssp585$legend == "map2"], na.rm = TRUE)
max_abs <- max(abs(range_values))

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp585
ggplot() +
  geom_raster(data = europe_mask_df, aes(x = x, y = y), fill = "gray30") +
  geom_raster(data = duration_trends_past_futssp585 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp585 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        limits = c(-max_abs, max_abs), breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(species ~ time, scales = "free_y",
              strip = strip_themed(background_y = elem_list_rect(fill = c("steelblue3", "steelblue3", "lightsteelblue1", "lightsteelblue1")),
                                   text_y = elem_list_text(face = c("bold.italic", NA, "bold.italic", NA)))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
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
    strip.background = element_rect(fill = "grey75", color = NA)
  )

# Save the figure
ggsave("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp585.png", width = 7.5, height = 11)




#-------------------------------------------------------------------------------
# 
# 4. Descriptive statistics  ---------------------------------------------------


# a) Extract maximum values of duration trend ----------------------------------

# Ixodes ricinus
slope_df_futssp370_I_ricinus <- slope_df_futssp370[slope_df_futssp370$species == "Ixodes ricinus", ]

max_trend_I_ricinus <- max(slope_df_futssp370_I_ricinus$trend, na.rm = TRUE)

max_trend_I_ricinus_rn <- which(slope_df_futssp370$species == "Ixodes ricinus" & slope_df_futssp370$trend == max_trend_I_ricinus)
slope_df_futssp370[max_trend_I_ricinus_rn, ]

# TBE
slope_df_futssp370_TBE <- slope_df_futssp370[slope_df_futssp370$species == "TBE", ]

max_trend_TBE <- max(slope_df_futssp370_TBE$trend, na.rm = TRUE)

max_trend_TBE_rn <- which(slope_df_futssp370$species == "TBE" & slope_df_futssp370$trend == max_trend_TBE)
slope_df_futssp370[max_trend_TBE_rn, ]

# Culex pipiens
slope_df_futssp370_C_pipiens <- slope_df_futssp370[slope_df_futssp370$species == "Culex pipiens", ]

max_trend_C_pipiens <- max(slope_df_futssp370_C_pipiens$trend, na.rm = TRUE)

max_trend_C_pipiens_rn <- which(slope_df_futssp370$species == "Culex pipiens" & slope_df_futssp370$trend == max_trend_C_pipiens)
slope_df_futssp370[max_trend_C_pipiens_rn, ]

# WNV
slope_df_futssp370_WNV <- slope_df_futssp370[slope_df_futssp370$species == "WNV", ]

max_trend_WNV <- max(slope_df_futssp370_WNV$trend, na.rm = TRUE)

max_trend_WNV_rn <- which(slope_df_futssp370$species == "WNV" & slope_df_futssp370$trend == max_trend_WNV)
slope_df_futssp370[max_trend_WNV_rn, ]




# b) Extract percentage of study area with positive duration trend -------------

# SSP370
# Ixodes ricinus 
total_rows_I_ricinus <- sum(slope_df_futssp370$species == "Ixodes ricinus")

positive_trend_I_ricinus <- sum(slope_df_futssp370$species == "Ixodes ricinus" & slope_df_futssp370$trend > 0)

percentage_positive_I_ricinus <- (positive_trend_I_ricinus / total_rows_I_ricinus) * 100

mean(slope_df_futssp370$species == "Ixodes ricinus" & slope_df_futssp370$trend > 0)

# TBE
total_rows_TBE <- sum(slope_df_futssp370$species == "TBE")

positive_trend_TBE <- sum(slope_df_futssp370$species == "TBE" & slope_df_futssp370$trend > 0)

percentage_positive_TBE <- (positive_trend_TBE / total_rows_TBE) * 100

mean(slope_df_futssp370$species == "TBE" & slope_df_futssp370$trend > 0)

# Culex pipiens
total_rows_C_pipiens <- sum(slope_df_futssp370$species == "Culex pipiens")

positive_trend_C_pipiens <- sum(slope_df_futssp370$species == "Culex pipiens" & slope_df_futssp370$trend > 0)

percentage_positive_C_pipiens <- (positive_trend_C_pipiens / total_rows_C_pipiens) * 100

mean(slope_df_futssp370$species == "Culex pipiens" & slope_df_futssp370$trend > 0)

# WNV
total_rows_WNV <- sum(slope_df_futssp370$species == "WNV")

positive_trend_WNV <- sum(slope_df_futssp370$species == "WNV" & slope_df_futssp370$trend > 0)

percentage_positive_WNV <- (positive_trend_WNV / total_rows_WNV) * 100

mean(slope_df_futssp370$species == "WNV" & slope_df_futssp370$trend > 0)

# SSP126
# Ixodes ricinus
total_rows_I_ricinus <- sum(slope_df_futssp126$species == "Ixodes ricinus")

positive_trend_I_ricinus <- sum(slope_df_futssp126$species == "Ixodes ricinus" & slope_df_futssp126$trend > 0)

percentage_positive_I_ricinus <- (positive_trend_I_ricinus / total_rows_I_ricinus) * 100


# TBE
total_rows_TBE <- sum(slope_df_futssp126$species == "TBE")

positive_trend_TBE <- sum(slope_df_futssp126$species == "TBE" & slope_df_futssp126$trend > 0)

percentage_positive_TBE <- (positive_trend_TBE / total_rows_TBE) * 100


# Culex pipiens
total_rows_C_pipiens <- sum(slope_df_futssp126$species == "Culex pipiens")

positive_trend_C_pipiens <- sum(slope_df_futssp126$species == "Culex pipiens" & slope_df_futssp126$trend > 0)

percentage_positive_C_pipiens <- (positive_trend_C_pipiens / total_rows_C_pipiens) * 100


# WNV
total_rows_WNV <- sum(slope_df_futssp126$species == "WNV")

positive_trend_WNV <- sum(slope_df_futssp126$species == "WNV" & slope_df_futssp126$trend > 0)

percentage_positive_WNV <- (positive_trend_WNV / total_rows_WNV) * 100


# SSP585
# Ixodes ricinus
total_rows_I_ricinus <- sum(slope_df_futssp585$species == "Ixodes ricinus")

positive_trend_I_ricinus <- sum(slope_df_futssp585$species == "Ixodes ricinus" & slope_df_futssp585$trend > 0)

percentage_positive_I_ricinus <- (positive_trend_I_ricinus / total_rows_I_ricinus) * 100


# TBE
total_rows_TBE <- sum(slope_df_futssp585$species == "TBE")

positive_trend_TBE <- sum(slope_df_futssp585$species == "TBE" & slope_df_futssp585$trend > 0)

percentage_positive_TBE <- (positive_trend_TBE / total_rows_TBE) * 100


# Culex pipiens
total_rows_C_pipiens <- sum(slope_df_futssp585$species == "Culex pipiens")

positive_trend_C_pipiens <- sum(slope_df_futssp585$species == "Culex pipiens" & slope_df_futssp585$trend > 0)

percentage_positive_C_pipiens <- (positive_trend_C_pipiens / total_rows_C_pipiens) * 100


# WNV
total_rows_WNV <- sum(slope_df_futssp585$species == "WNV")

positive_trend_WNV <- sum(slope_df_futssp585$species == "WNV" & slope_df_futssp585$trend > 0)

percentage_positive_WNV <- (positive_trend_WNV / total_rows_WNV) * 100






# Ixodes ricinus
decadal_means_2010s_I_ricinus <- decadal_means_2010s[decadal_means_2010s$species == "Ixodes ricinus", ]

mean(decadal_means_2010s_I_ricinus$trend)

# TBE
decadal_means_2010s_TBE <- decadal_means_2010s[decadal_means_2010s$species == "TBE", ]

mean(decadal_means_2010s_TBE$trend)

# Culex pipiens
decadal_means_2010s_C_pipiens <- decadal_means_2010s[decadal_means_2010s$species == "Culex pipiens", ]

total_rows_C_pipiens <- sum(decadal_means_2010s$species == "Culex pipiens")

predicted_absence_C_pipiens <- sum(decadal_means_2010s_C_pipiens$trend == 0)

percentage_absence_C_pipiens <- (predicted_absence_C_pipiens / total_rows_C_pipiens) * 100

# WNV
decadal_means_2010s_WNV <- decadal_means_2010s[decadal_means_2010s$species == "WNV", ]

total_rows_WNV <- sum(decadal_means_2010s$species == "WNV")

predicted_absence_WNV <- sum(decadal_means_2010s_WNV$trend == 0)

percentage_absence_WNV <- (predicted_absence_WNV / total_rows_WNV) * 100





#-------------------------------------------------------------------------------
# 
# for poster 
# c) Visualisation of duration trends for ssp370 -------------------------------

duration_trends_past_futssp370_1 <- duration_trends_past_futssp370[duration_trends_past_futssp370$species %in% c("Ixodes ricinus", "TBE"), ]

range_values <- range(duration_trends_past_futssp370_1$trend[duration_trends_past_futssp370_1$legend == "map2"], na.rm = TRUE)
max_abs <- max(abs(range_values))

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp370
ggplot() +
  geom_raster(data = europe_mask_df, aes(x = x, y = y), fill = "gray30") +
  geom_raster(data = duration_trends_past_futssp370_1 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp370_1 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        limits = c(-max_abs, max_abs), breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(species ~ time, scales = "free_y",
              strip = strip_themed(background_y = elem_list_rect(fill = c("steelblue3", "steelblue3")),
                                   text_y = elem_list_text(face = c("bold.italic", NA)))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
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
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

# Save the figure
ggsave(paste0("output_data/plots/duration_trends_2010s_Ixodes_TBE_ssp370_poster.png"), width = 7.5, height = 7, dpi = 400, bg = "transparent", device = grDevices::png,     # use base R's PNG function
       type = "cairo")



duration_trends_past_futssp370_2 <- duration_trends_past_futssp370[duration_trends_past_futssp370$species %in% c("Culex pipiens", "WNV"), ]

range_values <- range(duration_trends_past_futssp370_2$trend[duration_trends_past_futssp370_2$legend == "map2"], na.rm = TRUE)
max_abs <- max(abs(range_values))

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp370
ggplot() +
  geom_raster(data = europe_mask_df, aes(x = x, y = y), fill = "gray30") +
  geom_raster(data = duration_trends_past_futssp370_2 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp370_2 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        limits = c(-max_abs, max_abs), breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(species ~ time, scales = "free_y",
              strip = strip_themed(background_y = elem_list_rect(fill = c("lightsteelblue1", "lightsteelblue1")),
                                   text_y = elem_list_text(face = c("bold.italic", NA)))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
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
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

# Save the figure
ggsave(paste0("output_data/plots/duration_trends_2010s_Culex_WNV_ssp370_poster.png"), width = 7.5, height = 7, dpi = 400, bg = "transparent", device = grDevices::png,     # use base R's PNG function
       type = "cairo")

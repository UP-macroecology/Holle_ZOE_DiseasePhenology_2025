# ZOE disease phenology analysis


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#         15. Past and future duration trends of disease phenology       #
# ---------------------------------------------------------------------- #

# Load the needed packages
library(terra)
library(tidyterra)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(sf)
library(ggnewscale)
library(viridis)

# Load needed data
europe_mask_50km <- terra::rast("input_data/spatial_data/europe_mask_50km.tif")


#-------------------------------------------------------------------------------

# 1. Past duration trends based on disease phenology ---------------------------

# a) Load data -----------------------------------------------------------------

# Load needed data
# Read in the postprocessed monthly binary prediction data from 1970 to 2019
# under factual climate and land use
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_1970_2019.tif")

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif")

# TBE
TBE_occ_prob_clim_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_1970_2019.tif")

# WNV
WNV_occ_prob_clim_ens_bin <- terra::rast("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_1970_2019.tif")



# b) Calculation of trend values -----------------------------------------------

# For the Spatraster containing the predictions of the WNV, remove all layers with
# predictions made for the months January to April, as these predictions are too
# unreliable given that there were no disease occurrences reported
WNV_occ_prob_clim_ens_bin_names <- names(WNV_occ_prob_clim_ens_bin) # Extract layer names
WNV_occ_prob_clim_ens_bin_months <- substr(WNV_occ_prob_clim_ens_bin_names, 1, 2)
WNV_occ_prob_clim_ens_bin_keep <- which(!(WNV_occ_prob_clim_ens_bin_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_ens_bin <- WNV_occ_prob_clim_ens_bin[[WNV_occ_prob_clim_ens_bin_keep]]


# Extract the years from the layer names
I_ricinus_occ_prob_clim_landuse_ens_bin_names <- names(I_ricinus_occ_prob_clim_landuse_ens_bin)
I_ricinus_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(I_ricinus_occ_prob_clim_landuse_ens_bin_names, 4, 7))

C_pipiens_occ_prob_clim_landuse_ens_bin_names <- names(C_pipiens_occ_prob_clim_landuse_ens_bin)
C_pipiens_occ_prob_clim_landuse_ens_bin_years <- as.numeric(substr(C_pipiens_occ_prob_clim_landuse_ens_bin_names, 4, 7))

TBE_occ_prob_clim_ens_bin_names <- names(TBE_occ_prob_clim_ens_bin)
TBE_occ_prob_clim_ens_bin_years <- as.numeric(substr(TBE_occ_prob_clim_ens_bin_names, 4, 7))

WNV_occ_prob_clim_ens_bin_names <- names(WNV_occ_prob_clim_ens_bin)
WNV_occ_prob_clim_ens_bin_years <- as.numeric(substr(WNV_occ_prob_clim_ens_bin_names, 4, 7))



# Calculate the number of months with a predicted presence per year (for each cell)
I_ricinus_length_presence_year <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin, index = I_ricinus_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin, index = C_pipiens_occ_prob_clim_landuse_ens_bin_years, fun = sum, na.rm = TRUE)
TBE_length_presence_year <- terra::tapp(TBE_occ_prob_clim_ens_bin, index = TBE_occ_prob_clim_ens_bin_years, fun = sum, na.rm = TRUE)
WNV_length_presence_year <- terra::tapp(WNV_occ_prob_clim_ens_bin, index = WNV_occ_prob_clim_ens_bin_years, fun = sum, na.rm = TRUE)

# Rename yearly layers
names(I_ricinus_length_presence_year) <- paste0(unique(I_ricinus_occ_prob_clim_landuse_ens_bin_years))
names(C_pipiens_length_presence_year) <- paste0(unique(C_pipiens_occ_prob_clim_landuse_ens_bin_years))
names(TBE_length_presence_year) <- paste0(unique(TBE_occ_prob_clim_ens_bin_years))
names(WNV_length_presence_year) <- paste0(unique(WNV_occ_prob_clim_ens_bin_years))

# Create a vector with the studied decades
I_ricinus_r_names_curr <- names(I_ricinus_length_presence_year)
I_ricinus_r_decades_curr <- floor(as.numeric(I_ricinus_r_names_curr) / 10) * 10

C_pipiens_r_names_curr <- names(C_pipiens_length_presence_year)
C_pipiens_r_decades_curr <- floor(as.numeric(C_pipiens_r_names_curr) / 10) * 10

TBE_r_names_curr <- names(TBE_length_presence_year)
TBE_r_decades_curr <- floor(as.numeric(TBE_r_names_curr) / 10) * 10

WNV_r_names_curr <- names(WNV_length_presence_year)
WNV_r_decades_curr <- floor(as.numeric(WNV_r_names_curr) / 10) * 10


# Compute the mean presence length per decade (per cell)
I_ricinus_decadal_means_curr <- tapp(I_ricinus_length_presence_year, index = I_ricinus_r_decades_curr, fun = mean, na.rm = TRUE)
C_pipiens_decadal_means_curr <- tapp(C_pipiens_length_presence_year, index = C_pipiens_r_decades_curr, fun = mean, na.rm = TRUE)
TBE_decadal_means_curr <- tapp(TBE_length_presence_year, index = TBE_r_decades_curr, fun = mean, na.rm = TRUE)
WNV_decadal_means_curr <- tapp(WNV_length_presence_year, index = WNV_r_decades_curr, fun = mean, na.rm = TRUE)

# Rename decadal layers
names(I_ricinus_decadal_means_curr) <- paste0(unique(I_ricinus_r_decades_curr))
names(C_pipiens_decadal_means_curr) <- paste0(unique(C_pipiens_r_decades_curr))
names(TBE_decadal_means_curr) <- paste0(unique(TBE_r_decades_curr))
names(WNV_decadal_means_curr) <- paste0(unique(WNV_r_decades_curr))

# Extract the rasters containing the mean presence length per cell in the 
# decade of the 2010s for the species and diseases for later usage when
# calculating future trend
I_ricinus_decadal_means_curr_2010 <- I_ricinus_decadal_means_curr[["2010"]]
C_pipiens_decadal_means_curr_2010 <- C_pipiens_decadal_means_curr[["2010"]]
TBE_decadal_means_curr_2010 <- TBE_decadal_means_curr[["2010"]]
WNV_decadal_means_curr_2010 <- WNV_decadal_means_curr[["2010"]]

# Function to compute slope (trend) from 1970s to 2010s for each cell
calc_slope_curr <- function(x) {
  if (all(is.na(x))) return(NA) 
  decade_values <- seq(1970, 2010, by = 10)  # Decade midpoint years
  lm_fit <- lm(x ~ decade_values)  # Fit linear model
  return(coef(lm_fit)[2])  # Extract slope coefficient
}

# Apply function to raster to calculate trend (slope)
I_ricinus_slope_raster_curr <- app(I_ricinus_decadal_means_curr , fun = calc_slope_curr)
names(I_ricinus_slope_raster_curr) <- "Trend_Slope"

C_pipiens_slope_raster_curr <- app(C_pipiens_decadal_means_curr , fun = calc_slope_curr)
names(C_pipiens_slope_raster_curr) <- "Trend_Slope"

TBE_slope_raster_curr <- app(TBE_decadal_means_curr , fun = calc_slope_curr)
names(TBE_slope_raster_curr) <- "Trend_Slope"

WNV_slope_raster_curr <- app(WNV_decadal_means_curr , fun = calc_slope_curr)
names(WNV_slope_raster_curr) <- "Trend_Slope"



# c) Prepare data frame for plotting -------------------------------------------

# Convert raster to a data frame for plotting and add a column indicating the species/pathogen name
I_ricinus_slope_df_curr <- as.data.frame(I_ricinus_slope_raster_curr, xy = TRUE, na.rm = TRUE)
colnames(I_ricinus_slope_df_curr) <- c("lon", "lat", "trend")
I_ricinus_slope_df_curr$species <- "Ixodes ricinus"
I_ricinus_slope_df_curr$time <- "1970s - 2010s"

C_pipiens_slope_df_curr <- as.data.frame(C_pipiens_slope_raster_curr, xy = TRUE, na.rm = TRUE)
colnames(C_pipiens_slope_df_curr) <- c("lon", "lat", "trend")
C_pipiens_slope_df_curr$species <- "Culex pipiens"
C_pipiens_slope_df_curr$time <- "1970s - 2010s"

TBE_slope_df_curr <- as.data.frame(TBE_slope_raster_curr, xy = TRUE, na.rm = TRUE)
colnames(TBE_slope_df_curr) <- c("lon", "lat", "trend")
TBE_slope_df_curr$species <- "TBE"
TBE_slope_df_curr$time <- "1970s - 2010s"

WNV_slope_df_curr <- as.data.frame(WNV_slope_raster_curr, xy = TRUE, na.rm = TRUE)
colnames(WNV_slope_df_curr) <- c("lon", "lat", "trend")
WNV_slope_df_curr$species <- "WNV"
WNV_slope_df_curr$time <- "1970s - 2010s"

# Bind the four different data frames
slope_df_past <- rbind(I_ricinus_slope_df_curr, C_pipiens_slope_df_curr,
                       TBE_slope_df_curr, WNV_slope_df_curr)




# d) Save resulting data frame --------------------------------------------------

save(slope_df_past, file = "output_data/results/duration_trends/duration_trends_past_vector_virus.RData")




#-------------------------------------------------------------------------------

# 2. Future duration trends based on disease phenology -------------------------


# a) Load data -----------------------------------------------------------------

# Load needed data - postprocessed monthly predicted ensemble occurrence probability of main vectors and
# the respective viruses, under climate and land use change for the
# future years 2030 to 2070 based and three studied environmental scenarios
# Ixodes ricinus; ssp126, ssp370, and ssp585
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp126.tif"))
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp585.tif"))

# Culex pipiens; ssp126, ssp370, and ssp585
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp585.tif"))

# TBE; ssp126, ssp370, and ssp585
TBE_occ_prob_clim_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp126.tif"))
TBE_occ_prob_clim_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp370.tif"))
TBE_occ_prob_clim_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp585.tif"))

# WNV; ssp126, ssp370, and ssp585
WNV_occ_prob_clim_ens_bin_fut_ssp126 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp126.tif"))
WNV_occ_prob_clim_ens_bin_fut_ssp370 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp370.tif"))
WNV_occ_prob_clim_ens_bin_fut_ssp585 <- terra::rast(paste0("output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp585.tif"))



# b) Calculation of trend values -----------------------------------------------

# For the Spatraster containing the predictions of the WNV, remove all layers with
# predictions made for the months January to May, as these predictions are too
# unreliable given that there were no disease occurrences reported
WNV_occ_prob_clim_ens_bin_fut_ssp126_names <- names(WNV_occ_prob_clim_ens_bin_fut_ssp126) # Extract layer names
WNV_occ_prob_clim_ens_bin_fut_ssp126_months <- substr(WNV_occ_prob_clim_ens_bin_fut_ssp126_names, 1, 2)
WNV_occ_prob_clim_ens_bin_fut_ssp126_keep <- which(!(WNV_occ_prob_clim_ens_bin_fut_ssp126_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_ens_bin_fut_ssp126[[WNV_occ_prob_clim_ens_bin_fut_ssp126_keep]]

WNV_occ_prob_clim_ens_bin_fut_ssp370_names <- names(WNV_occ_prob_clim_ens_bin_fut_ssp370) # Extract layer names
WNV_occ_prob_clim_ens_bin_fut_ssp370_months <- substr(WNV_occ_prob_clim_ens_bin_fut_ssp370_names, 1, 2)
WNV_occ_prob_clim_ens_bin_fut_ssp370_keep <- which(!(WNV_occ_prob_clim_ens_bin_fut_ssp370_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_ens_bin_fut_ssp370[[WNV_occ_prob_clim_ens_bin_fut_ssp370_keep]]

WNV_occ_prob_clim_ens_bin_fut_ssp585_names <- names(WNV_occ_prob_clim_ens_bin_fut_ssp585) # Extract layer names
WNV_occ_prob_clim_ens_bin_fut_ssp585_months <- substr(WNV_occ_prob_clim_ens_bin_fut_ssp585_names, 1, 2)
WNV_occ_prob_clim_ens_bin_fut_ssp585_keep <- which(!(WNV_occ_prob_clim_ens_bin_fut_ssp585_months %in% c("01", "02", "03", "04")))
WNV_occ_prob_clim_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_ens_bin_fut_ssp126[[WNV_occ_prob_clim_ens_bin_fut_ssp585_keep]]


# Extract the years from the layer names
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

TBE_r_names_fut_ssp126 <- names(TBE_occ_prob_clim_ens_bin_fut_ssp126)
TBE_r_years_fut_ssp126 <- as.numeric(substr(TBE_r_names_fut_ssp126, 4, 7))
TBE_r_names_fut_ssp370 <- names(TBE_occ_prob_clim_ens_bin_fut_ssp370)
TBE_r_years_fut_ssp370 <- as.numeric(substr(TBE_r_names_fut_ssp370, 4, 7))
TBE_r_names_fut_ssp585 <- names(TBE_occ_prob_clim_ens_bin_fut_ssp585)
TBE_r_years_fut_ssp585 <- as.numeric(substr(TBE_r_names_fut_ssp585, 4, 7))

WNV_r_names_fut_ssp126 <- names(WNV_occ_prob_clim_ens_bin_fut_ssp126)
WNV_r_years_fut_ssp126 <- as.numeric(substr(WNV_r_names_fut_ssp126, 4, 7))
WNV_r_names_fut_ssp370 <- names(WNV_occ_prob_clim_ens_bin_fut_ssp370)
WNV_r_years_fut_ssp370 <- as.numeric(substr(WNV_r_names_fut_ssp370, 4, 7))
WNV_r_names_fut_ssp585 <- names(WNV_occ_prob_clim_ens_bin_fut_ssp585)
WNV_r_years_fut_ssp585 <- as.numeric(substr(WNV_r_names_fut_ssp585, 4, 7))

# Calculate the number of months with a predicted presence per year (for each cell)
I_ricinus_length_presence_year_fut_ssp126 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = I_ricinus_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
I_ricinus_length_presence_year_fut_ssp370 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = I_ricinus_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
I_ricinus_length_presence_year_fut_ssp585 <- terra::tapp(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = I_ricinus_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

C_pipiens_length_presence_year_fut_ssp126 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126, index = C_pipiens_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year_fut_ssp370 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370, index = C_pipiens_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
C_pipiens_length_presence_year_fut_ssp585 <- terra::tapp(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585, index = C_pipiens_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

TBE_length_presence_year_fut_ssp126 <- terra::tapp(TBE_occ_prob_clim_ens_bin_fut_ssp126, index = TBE_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
TBE_length_presence_year_fut_ssp370 <- terra::tapp(TBE_occ_prob_clim_ens_bin_fut_ssp370, index = TBE_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
TBE_length_presence_year_fut_ssp585 <- terra::tapp(TBE_occ_prob_clim_ens_bin_fut_ssp585, index = TBE_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

WNV_length_presence_year_fut_ssp126 <- terra::tapp(WNV_occ_prob_clim_ens_bin_fut_ssp126, index = WNV_r_years_fut_ssp126, fun = sum, na.rm = TRUE)
WNV_length_presence_year_fut_ssp370 <- terra::tapp(WNV_occ_prob_clim_ens_bin_fut_ssp370, index = WNV_r_years_fut_ssp370, fun = sum, na.rm = TRUE)
WNV_length_presence_year_fut_ssp585 <- terra::tapp(WNV_occ_prob_clim_ens_bin_fut_ssp585, index = WNV_r_years_fut_ssp585, fun = sum, na.rm = TRUE)

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
I_ricinus_decadal_means_fut_ssp126 <- c(I_ricinus_decadal_means_curr_2010, I_ricinus_decadal_means_fut_ssp126)
I_ricinus_decadal_means_fut_ssp370 <- c(I_ricinus_decadal_means_curr_2010, I_ricinus_decadal_means_fut_ssp370)
I_ricinus_decadal_means_fut_ssp585 <- c(I_ricinus_decadal_means_curr_2010, I_ricinus_decadal_means_fut_ssp585)

C_pipiens_decadal_means_fut_ssp126 <- c(C_pipiens_decadal_means_curr_2010, C_pipiens_decadal_means_fut_ssp126)
C_pipiens_decadal_means_fut_ssp370 <- c(C_pipiens_decadal_means_curr_2010, C_pipiens_decadal_means_fut_ssp370)
C_pipiens_decadal_means_fut_ssp585 <- c(C_pipiens_decadal_means_curr_2010, C_pipiens_decadal_means_fut_ssp585)

TBE_decadal_means_fut_ssp126 <- c(TBE_decadal_means_curr_2010, TBE_decadal_means_fut_ssp126)
TBE_decadal_means_fut_ssp370 <- c(TBE_decadal_means_curr_2010, TBE_decadal_means_fut_ssp370)
TBE_decadal_means_fut_ssp585 <- c(TBE_decadal_means_curr_2010, TBE_decadal_means_fut_ssp585)

WNV_decadal_means_fut_ssp126 <- c(WNV_decadal_means_curr_2010, WNV_decadal_means_fut_ssp126)
WNV_decadal_means_fut_ssp370 <- c(WNV_decadal_means_curr_2010, WNV_decadal_means_fut_ssp370)
WNV_decadal_means_fut_ssp585 <- c(WNV_decadal_means_curr_2010, WNV_decadal_means_fut_ssp585)


# Function to compute slope (trend) from 1910s to 2050s for each cell
calc_slope_fut <- function(x) {
  if (all(is.na(x))) return(NA)  # Handle NA pixels
  decade_values <- c(2010, 2030, 2040, 2050)  # Decade midpoint years
  lm_fit <- lm(x ~ decade_values)  # Fit linear model
  return(coef(lm_fit)[2])  # Extract slope coefficient
}

# Filter the Spatrasters to only contain decades from 2030 to 2060
I_ricinus_decadal_means_fut_ssp126 <- subset(I_ricinus_decadal_means_fut_ssp126, names(I_ricinus_decadal_means_fut_ssp126) %in% c("2010", "2030", "2040", "2050"))
I_ricinus_decadal_means_fut_ssp370 <- subset(I_ricinus_decadal_means_fut_ssp370, names(I_ricinus_decadal_means_fut_ssp370) %in% c("2010", "2030", "2040", "2050"))
I_ricinus_decadal_means_fut_ssp585 <- subset(I_ricinus_decadal_means_fut_ssp585, names(I_ricinus_decadal_means_fut_ssp585) %in% c("2010", "2030", "2040", "2050"))

C_pipiens_decadal_means_fut_ssp126 <- subset(C_pipiens_decadal_means_fut_ssp126, names(C_pipiens_decadal_means_fut_ssp126) %in% c("2010", "2030", "2040", "2050"))
C_pipiens_decadal_means_fut_ssp370 <- subset(C_pipiens_decadal_means_fut_ssp370, names(C_pipiens_decadal_means_fut_ssp370) %in% c("2010", "2030", "2040", "2050"))
C_pipiens_decadal_means_fut_ssp585 <- subset(C_pipiens_decadal_means_fut_ssp585, names(C_pipiens_decadal_means_fut_ssp585) %in% c("2010", "2030", "2040", "2050"))

TBE_decadal_means_fut_ssp126 <- subset(TBE_decadal_means_fut_ssp126, names(TBE_decadal_means_fut_ssp126) %in% c("2010", "2030", "2040", "2050"))
TBE_decadal_means_fut_ssp370 <- subset(TBE_decadal_means_fut_ssp370, names(TBE_decadal_means_fut_ssp370) %in% c("2010", "2030", "2040", "2050"))
TBE_decadal_means_fut_ssp585 <- subset(TBE_decadal_means_fut_ssp585, names(TBE_decadal_means_fut_ssp585) %in% c("2010", "2030", "2040", "2050"))

WNV_decadal_means_fut_ssp126 <- subset(WNV_decadal_means_fut_ssp126, names(WNV_decadal_means_fut_ssp126) %in% c("2010", "2030", "2040", "2050"))
WNV_decadal_means_fut_ssp370 <- subset(WNV_decadal_means_fut_ssp370, names(WNV_decadal_means_fut_ssp370) %in% c("2010", "2030", "2040", "2050"))
WNV_decadal_means_fut_ssp585 <- subset(WNV_decadal_means_fut_ssp585, names(WNV_decadal_means_fut_ssp585) %in% c("2010", "2030", "2040", "2050"))

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
# showing duration of 2010s and future trends


# a) Prepare resulting past and future trend data frame ------------------------

# Convert raster with mean duration trends in the 2010s into a data frame
# and add important columns
I_ricinus_decadal_means_curr_2010_df <- as.data.frame(I_ricinus_decadal_means_curr_2010, xy = TRUE)
colnames(I_ricinus_decadal_means_curr_2010_df) <- c("lon", "lat", "trend")
I_ricinus_decadal_means_curr_2010_df$species <- "Ixodes ricinus"
I_ricinus_decadal_means_curr_2010_df$time <- "2010s"
I_ricinus_decadal_means_curr_2010_df$legend <- "map1"


C_pipiens_decadal_means_curr_2010_df <- as.data.frame(C_pipiens_decadal_means_curr_2010, xy = TRUE)
colnames(C_pipiens_decadal_means_curr_2010_df) <- c("lon", "lat", "trend")
C_pipiens_decadal_means_curr_2010_df$species <- "Culex pipiens"
C_pipiens_decadal_means_curr_2010_df$time <- "2010s"
C_pipiens_decadal_means_curr_2010_df$legend <- "map1"

TBE_decadal_means_curr_2010_df <- as.data.frame(TBE_decadal_means_curr_2010, xy = TRUE)
colnames(TBE_decadal_means_curr_2010_df) <- c("lon", "lat", "trend")
TBE_decadal_means_curr_2010_df$species <- "TBE"
TBE_decadal_means_curr_2010_df$time <- "2010s"
TBE_decadal_means_curr_2010_df$legend <- "map1"

WNV_decadal_means_curr_2010_df <- as.data.frame(WNV_decadal_means_curr_2010, xy = TRUE)
colnames(WNV_decadal_means_curr_2010_df) <- c("lon", "lat", "trend")
WNV_decadal_means_curr_2010_df$species <- "WNV"
WNV_decadal_means_curr_2010_df$time <- "2010s"
WNV_decadal_means_curr_2010_df$legend <- "map1"

# Bind the four data frames of the different species/pathogens
decadal_means_2010s <- rbind(I_ricinus_decadal_means_curr_2010_df,
                             C_pipiens_decadal_means_curr_2010_df,
                             TBE_decadal_means_curr_2010_df,
                             WNV_decadal_means_curr_2010_df)

# Add a legend column for the data frames containing future predictions
slope_df_futssp126$legend <- "map2"
slope_df_futssp370$legend <- "map2"
slope_df_futssp585$legend <- "map2"

# Bind the data frames
duration_trends_past_futssp126 <- rbind(decadal_means_2010s, slope_df_futssp126)
duration_trends_past_futssp370 <- rbind(decadal_means_2010s, slope_df_futssp370)
duration_trends_past_futssp585 <- rbind(decadal_means_2010s, slope_df_futssp585)


# Make sure the vector and diseases appear in the correct position
duration_trends_past_futssp126$species <- factor(duration_trends_past_futssp126$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))

duration_trends_past_futssp370$species <- factor(duration_trends_past_futssp370$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))

duration_trends_past_futssp585$species <- factor(duration_trends_past_futssp585$species, 
                                                 levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))


# Convert Europe mask spatraster into a data frame
europe_mask_50km_df <- as.data.frame(europe_mask_50km, xy = TRUE)





# b) Visualisation of duration trends for ssp126 -------------------------------
  
  # Plot the duration trend per cell based on the factual historical data and the 
  # environmental scenario ssp126
  ggplot() +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster(data = duration_trends_past_futssp126 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp126 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
              strip = strip_themed(background_y = list(
                "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                "TBE" = element_rect(fill = "steelblue3"),
                "Culex pipiens" = element_rect(fill = "lightsteelblue1"),
                "WNV" = element_rect(fill = "lightsteelblue1")
              ))) +
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
ggsave(paste0("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp126.png"), width = 7.5, height = 11)




# c) Visualisation of duration trends for ssp370 -------------------------------

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp370
ggplot() +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster(data = duration_trends_past_futssp370 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp370 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
              strip = strip_themed(background_y = list(
                "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                "TBE" = element_rect(fill = "steelblue3"),
                "Culex pipiens" = element_rect(fill = "lightsteelblue1"),
                "WNV" = element_rect(fill = "lightsteelblue1")
              ))) +
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
ggsave(paste0("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp370.png"), width = 7.5, height = 11)




# d) Visualisation of duration trends for ssp585 -------------------------------

# Plot the duration trend per cell based on the factual historical data and the 
# environmental scenario ssp585
ggplot() +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster(data = duration_trends_past_futssp585 %>% filter(legend == "map1"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_viridis_c(option = "F", direction = -1, name = "Duration (months)",
                       guide = guide_colorbar(order = 1)) +
  ggnewscale::new_scale_fill() +
  geom_raster(data = duration_trends_past_futssp585 %>% filter(legend == "map2"),
              aes(x = lon, y = lat, fill = trend)) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        breaks = c(-0.2, 0, 0.2), guide = guide_colorbar(order = 2)) +
  facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
              strip = strip_themed(background_y = list(
                "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                "TBE" = element_rect(fill = "steelblue3"),
                "Culex pipiens" = element_rect(fill = "lightsteelblue1"),
                "WNV" = element_rect(fill = "lightsteelblue1")
              ))) +
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
ggsave(paste0("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_Culex_WNV_ssp585.png"), width = 7.5, height = 11)







# #-------------------------------------------------------------------------------
# 
# # 3. Visualise duration trends  ------------------------------------------------
# # for main vectors as well as viruses
# 
# duration_trends_past_futssp126 <- rbind(slope_df_past, slope_df_futssp126)
# duration_trends_past_futssp370 <- rbind(slope_df_past, slope_df_futssp370)
# duration_trends_past_futssp585 <- rbind(slope_df_past, slope_df_futssp585)
# 
# # Make sure the vector and diseases appear in the correct position
# duration_trends_past_futssp126$species <- factor(duration_trends_past_futssp126$species, 
#                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
# 
# duration_trends_past_futssp370$species <- factor(duration_trends_past_futssp370$species, 
#                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
# 
# duration_trends_past_futssp585$species <- factor(duration_trends_past_futssp585$species, 
#                                                  levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))
# 
# # Convert Europe mask spatraster into a data frame
# europe_mask_50km_df <- as.data.frame(europe_mask_50km, xy = TRUE)
# 
# # Plot the duration trend per cell for the observed historical data and the 
# # environmental scenario ssp126
# ggplot(duration_trends_past_futssp126) +
#   geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
#   geom_raster(data = duration_trends_past_futssp126, aes(x = lon, y = lat, fill = trend)) +
#   scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)") +
#   facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
#               strip = strip_themed(background_y = list(
#                 "Ixodes ricinus" = element_rect(fill = "steelblue3"),
#                 "TBE" = element_rect(fill = "steelblue3"),            
#                 "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
#                 "WNV" = element_rect(fill = "lightsteelblue1")  
#               ))) +
#   theme_bw() +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     text = element_text(size = 14),  
#     axis.title = element_text(size = 14),  
#     axis.text = element_text(size = 12),  
#     legend.title = element_text(size = 12.5, face = "bold"),  
#     legend.text = element_text(size = 12),  
#     plot.title = element_text(size = 16, face = "bold"),
#     strip.text = element_text(size = 16, face = "bold"), 
#     strip.background = element_rect(fill = "grey75", color = NA)) +
#   guides(
#     color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
#     linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
#   )
# 
# ggsave(paste0("output_data/plots/duration_trends/duration_trends_vector_virus_ssp126.png"), width = 8.5, height = 12)
# 
# 
# # Plot the duration trend per cell for the observed historical data and the 
# # environmental scenario ssp1370
# ggplot(duration_trends_past_futssp370) +
#   geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
#   geom_raster(data = duration_trends_past_futssp370, aes(x = lon, y = lat, fill = trend)) +
#   scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)") +
#   facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
#               strip = strip_themed(background_y = list(
#                 "Ixodes ricinus" = element_rect(fill = "steelblue3"),
#                 "TBE" = element_rect(fill = "steelblue3"),            
#                 "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
#                 "WNV" = element_rect(fill = "lightsteelblue1")  
#               ))) +
#   theme_bw() +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     text = element_text(size = 14),  
#     axis.title = element_text(size = 14),  
#     axis.text = element_text(size = 12),  
#     legend.title = element_text(size = 12.5, face = "bold"),  
#     legend.text = element_text(size = 12),  
#     plot.title = element_text(size = 16, face = "bold"),
#     strip.text = element_text(size = 16, face = "bold"), 
#     strip.background = element_rect(fill = "grey75", color = NA)) +
#   guides(
#     color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
#     linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
#   )
# 
# ggsave(paste0("output_data/plots/duration_trends/duration_trends_vector_virus_ssp370.png"), width = 8.5, height = 12)
# 
# 
# 
# # Plot the duration trend per cell for the observed historical data and the 
# # environmental scenario ssp585
# ggplot(duration_trends_past_futssp585) +
#   geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
#   geom_raster(data = duration_trends_past_futssp585, aes(x = lon, y = lat, fill = trend)) +
#   scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)") +
#   facet_grid2(rows = vars(species), cols = vars(time), scales = "free_y",
#               strip = strip_themed(background_y = list(
#                 "Ixodes ricinus" = element_rect(fill = "steelblue3"),
#                 "TBE" = element_rect(fill = "steelblue3"),            
#                 "Culex pipiens" = element_rect(fill = "lightsteelblue1"),   
#                 "WNV" = element_rect(fill = "lightsteelblue1")  
#               ))) +
#   theme_bw() +
#   labs(x = "Longitude", y = "Latitude") +
#   theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     text = element_text(size = 14),  
#     axis.title = element_text(size = 14),  
#     axis.text = element_text(size = 12),  
#     legend.title = element_text(size = 12.5, face = "bold"),  
#     legend.text = element_text(size = 12),  
#     plot.title = element_text(size = 16, face = "bold"),
#     strip.text = element_text(size = 16, face = "bold"), 
#     strip.background = element_rect(fill = "grey75", color = NA)) +
#   guides(
#     color = guide_legend(title.position = "top", nrow = 3, byrow = TRUE),  
#     linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE) 
#   )
# 
# ggsave(paste0("output_data/plots/duration_trends/duration_trends_vector_virus_ssp585.png"), width = 8.5, height = 12)

# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                   13. Decadal trends in distribution                   #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We visualise the distribution trends of the vectors and viruses for one target
# month - May for Ixodes ricinus and TBE, and July for Culex pipiens and WNV - 
# across three different target decades: 1970s, 2010, 2050s. For the future 
# target decade (2050s), we look at predictions the were derived from the three 
# different socio-economic scenarios. Only grid cells with vector or virus 
# suitability values indicating at least one predicted presence within the 
# selected month of a given decade are displayed. 


# Load needed packages:
library(terra) # terra_1.7-55
library(ggplot2) # ggplot2_4.0.0
library(tidyverse) # tidyverse_2.0.0
library(ggnewscale) # ggnewscale_0.5.1
library(ggh4x) # ggh4x_0.3.1


# Load needed data
load("output_data/validation/I_ricinus_validation.RData") # Load validation results to obtain binary predictions
I_ricinus_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"] # Extract ensemble threshold

load("output_data/validation/C_pipiens_validation.RData")
C_pipiens_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]

load("output_data/validation/TBE_validation.RData") 
TBE_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]

load("output_data/validation/WNV_validation.RData")
WNV_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]




#-------------------------------------------------------------------------------

# 1. Historical distribution trends --------------------------------------------

# a) Load data -----------------------------------------------------------------


# Read in historical ensemble prediction outputs for the vector and viruses
# based on factual climate and land-use change
I_ricinus_r_past_preds_clim_landuse <- terra::rast("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_1970_2019.tif")
C_pipiens_r_past_preds_clim_landuse <- terra::rast("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")
TBE_r_past_preds_clim_landuse <- terra::rast("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_1970_2019.tif")
WNV_r_past_preds_clim_landuse <- terra::rast("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_1970_2019.tif")




# b) Prepare historical prediction rasters for plotting ------------------------
# 1970s and 2010s as target decades

# Write a function that extracts the layers for the years 1970 to 1979 (1970s)
get_decade_layers_1970s <- function(r){
  yrs <- as.numeric(sub(".*/", "", names(r))) 
  sel <- yrs >= 1970 & yrs <= 1979
  r[[sel]]
}

# Write a function that extracts the layers for the years 2010 to 2019 (2010s)
get_decade_layers_2010s <- function(r){
  yrs <- as.numeric(sub(".*/", "", names(r)))
  sel <- yrs >= 2010 & yrs <= 2019
  r[[sel]]
}

# Extract the layers for the two different decades
# Ixodes ricinus
I_ricinus_r_past_preds_clim_landuse_1970s <- get_decade_layers_1970s(I_ricinus_r_past_preds_clim_landuse)
I_ricinus_r_past_preds_clim_landuse_2010s <- get_decade_layers_2010s(I_ricinus_r_past_preds_clim_landuse)

# Culex pipiens
C_pipiens_r_past_preds_clim_landuse_1970s <- get_decade_layers_1970s(C_pipiens_r_past_preds_clim_landuse)
C_pipiens_r_past_preds_clim_landuse_2010s <- get_decade_layers_2010s(C_pipiens_r_past_preds_clim_landuse)

# TBE
TBE_r_past_preds_clim_landuse_1970s <- get_decade_layers_1970s(TBE_r_past_preds_clim_landuse)
TBE_r_past_preds_clim_landuse_2010s <- get_decade_layers_2010s(TBE_r_past_preds_clim_landuse)

# WNV
WNV_r_past_preds_clim_landuse_1970s <- get_decade_layers_1970s(WNV_r_past_preds_clim_landuse)
WNV_r_past_preds_clim_landuse_2010s <- get_decade_layers_2010s(WNV_r_past_preds_clim_landuse)




# Write a function that extracts the decade layers for only one specific month
keep_month <- function(r, month_code){
  layer_names <- names(r)
  months <- substr(layer_names, 1, 2)
  r[[months == month_code]]
}

# Apply the month filtering (May for Ixodes ricinus and TBE; July for Culex
# pipiens and WNV)
I_ricinus_r_past_preds_clim_landuse_1970s <- keep_month(I_ricinus_r_past_preds_clim_landuse_1970s, "05") # May
I_ricinus_r_past_preds_clim_landuse_2010s <- keep_month(I_ricinus_r_past_preds_clim_landuse_2010s, "05") # May
TBE_r_past_preds_clim_landuse_1970s <- keep_month(TBE_r_past_preds_clim_landuse_1970s, "05") # May
TBE_r_past_preds_clim_landuse_2010s <- keep_month(TBE_r_past_preds_clim_landuse_2010s, "05") # May 

C_pipiens_r_past_preds_clim_landuse_1970s <- keep_month(C_pipiens_r_past_preds_clim_landuse_1970s, "07") # July
C_pipiens_r_past_preds_clim_landuse_2010s <- keep_month(C_pipiens_r_past_preds_clim_landuse_2010s, "07") # July
WNV_r_past_preds_clim_landuse_1970s <- keep_month(WNV_r_past_preds_clim_landuse_1970s, "07") # July
WNV_r_past_preds_clim_landuse_2010s <- keep_month(WNV_r_past_preds_clim_landuse_2010s, "07") # July


# Write a function that applies the threshold for binarising predictions, then 
# averages the predicted values for each cell over the decade, this is only done
# for cells that show at least one presence for the selected month within the
# decade ("ever present" cells)
threshold_filter_mean <- function(r, threshold){
  r_above_thr <- r >= threshold # 0/1 presence for selected month
  ever_present <- max(r_above_thr, na.rm = TRUE) == 1 # Cells with at least one presence
  r_mean <- mean(r, na.rm = TRUE) # Mean prediction across month
  r_mean[!ever_present] <- NA # Mask cells never above threshold (NA)
  r_mean
}

# Apply function for the vector and viruses for the two different decades
I_ricinus_mean_thr_1970s <- threshold_filter_mean(I_ricinus_r_past_preds_clim_landuse_1970s, I_ricinus_thresh)
I_ricinus_mean_thr_2010s <- threshold_filter_mean(I_ricinus_r_past_preds_clim_landuse_2010s, I_ricinus_thresh)

C_pipiens_mean_thr_1970s <- threshold_filter_mean(C_pipiens_r_past_preds_clim_landuse_1970s, C_pipiens_thresh)
C_pipiens_mean_thr_2010s <- threshold_filter_mean(C_pipiens_r_past_preds_clim_landuse_2010s, C_pipiens_thresh)

TBE_mean_thr_1970s <- threshold_filter_mean(TBE_r_past_preds_clim_landuse_1970s, TBE_thresh)
TBE_mean_thr_2010s <- threshold_filter_mean(TBE_r_past_preds_clim_landuse_2010s, TBE_thresh)

WNV_mean_thr_1970s <- threshold_filter_mean(WNV_r_past_preds_clim_landuse_1970s, WNV_thresh)
WNV_mean_thr_2010s <- threshold_filter_mean(WNV_r_past_preds_clim_landuse_2010s, WNV_thresh)



# Convert rasters to data frames for plotting
I_ricinus_df_1970s <- as.data.frame(I_ricinus_mean_thr_1970s, xy = TRUE)
I_ricinus_df_2010s <- as.data.frame(I_ricinus_mean_thr_2010s, xy = TRUE)

C_pipiens_df_1970s <- as.data.frame(C_pipiens_mean_thr_1970s, xy = TRUE)
C_pipiens_df_2010s <- as.data.frame(C_pipiens_mean_thr_2010s, xy = TRUE)

TBE_df_1970s <- as.data.frame(TBE_mean_thr_1970s, xy = TRUE)
TBE_df_2010s <- as.data.frame(TBE_mean_thr_2010s, xy = TRUE)

WNV_df_1970s  <- as.data.frame(WNV_mean_thr_1970s,  xy = TRUE)
WNV_df_2010s  <- as.data.frame(WNV_mean_thr_2010s,  xy = TRUE)

# Add a column indicating the studied decade, the studied vector/virus, the 
# vector/virus combination (as we plot them in one figure)
I_ricinus_df_1970s$time <- "1970s"
I_ricinus_df_1970s$species <- "Ixodes ricinus"
I_ricinus_df_1970s$combination <- "Ixodes ricinus & TBE"

I_ricinus_df_2010s$time <- "2010s"
I_ricinus_df_2010s$species <- "Ixodes ricinus"
I_ricinus_df_2010s$combination <- "Ixodes ricinus & TBE"

C_pipiens_df_1970s$time <- "1970s"
C_pipiens_df_1970s$species <- "Culex pipiens"
C_pipiens_df_1970s$combination <- "Culex pipiens & WNV"

C_pipiens_df_2010s$time <- "2010s"
C_pipiens_df_2010s$species <- "Culex pipiens"
C_pipiens_df_2010s$combination <- "Culex pipiens & WNV"

TBE_df_1970s$time <- "1970s"
TBE_df_1970s$species <- "TBE"
TBE_df_1970s$combination <- "Ixodes ricinus & TBE"

TBE_df_2010s$time <- "2010s"
TBE_df_2010s$species <- "TBE"
TBE_df_2010s$combination <- "Ixodes ricinus & TBE"

WNV_df_1970s$time <- "1970s"
WNV_df_1970s$species <- "WNV"
WNV_df_1970s$combination <- "Culex pipiens & WNV"

WNV_df_2010s$time <- "2010s"
WNV_df_2010s$species <- "WNV"
WNV_df_2010s$combination <- "Culex pipiens & WNV"


# Combine the prepared data frames into one
aggregated_df_past <- bind_rows(I_ricinus_df_1970s, I_ricinus_df_2010s,
                                C_pipiens_df_1970s, C_pipiens_df_2010s,
                                TBE_df_1970s, TBE_df_2010s,
                                WNV_df_1970s, WNV_df_2010s)

# Save the data frame
save(aggregated_df_past, file = "output_data/results/distribution_trends/distribution_trends_past_vector_virus.RData")




#-------------------------------------------------------------------------------

# 2. Future distribution trends ------------------------------------------------

# a) Load data -----------------------------------------------------------------

# Read in future ensemble prediction outputs for the vector and viruses
# based on the three different socio-economic scenarios
# Ixodes ricinus
I_ricinus_r_fut_preds_clim_landuse_ssp126 <- terra::rast("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp126.tif")
I_ricinus_r_fut_preds_clim_landuse_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp370.tif")
I_ricinus_r_fut_preds_clim_landuse_ssp585 <- terra::rast("output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp585.tif")

# Culex pipiens
C_pipiens_r_fut_preds_clim_landuse_ssp126 <- terra::rast("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp126.tif")
C_pipiens_r_fut_preds_clim_landuse_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp370.tif")
C_pipiens_r_fut_preds_clim_landuse_ssp585 <- terra::rast("output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp585.tif")

# TBE
TBE_r_fut_preds_clim_landuse_ssp126 <- terra::rast("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp126.tif")
TBE_r_fut_preds_clim_landuse_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp370.tif")
TBE_r_fut_preds_clim_landuse_ssp585 <- terra::rast("output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp585.tif")

# WNV
WNV_r_fut_preds_clim_landuse_ssp126 <- terra::rast("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp126.tif")
WNV_r_fut_preds_clim_landuse_ssp370 <- terra::rast("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp370.tif")
WNV_r_fut_preds_clim_landuse_ssp585 <- terra::rast("output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp585.tif")




# b) Prepare future prediction rasters for plotting ----------------------------
# 2050s as target decade

# Write a function that extracts the layers for the years 2050 to 2059 (2050s)
get_decade_layers_2050s <- function(r){
  yrs <- as.numeric(sub(".*/", "", names(r))) 
  sel <- yrs >= 2050 & yrs <= 2059
  r[[sel]]
}

# Extract the layers for the future target decade
# Ixodes ricinus
I_ricinus_r_fut_preds_clim_landuse_ssp126_2050s <- get_decade_layers_2050s(I_ricinus_r_fut_preds_clim_landuse_ssp126)
I_ricinus_r_fut_preds_clim_landuse_ssp370_2050s <- get_decade_layers_2050s(I_ricinus_r_fut_preds_clim_landuse_ssp370)
I_ricinus_r_fut_preds_clim_landuse_ssp585_2050s <- get_decade_layers_2050s(I_ricinus_r_fut_preds_clim_landuse_ssp585)

# Culex pipiens
C_pipiens_r_fut_preds_clim_landuse_ssp126_2050s <- get_decade_layers_2050s(C_pipiens_r_fut_preds_clim_landuse_ssp126)
C_pipiens_r_fut_preds_clim_landuse_ssp370_2050s <- get_decade_layers_2050s(C_pipiens_r_fut_preds_clim_landuse_ssp370)
C_pipiens_r_fut_preds_clim_landuse_ssp585_2050s <- get_decade_layers_2050s(C_pipiens_r_fut_preds_clim_landuse_ssp585)

# TBE
TBE_r_fut_preds_clim_landuse_ssp126_2050s <- get_decade_layers_2050s(TBE_r_fut_preds_clim_landuse_ssp126)
TBE_r_fut_preds_clim_landuse_ssp370_2050s <- get_decade_layers_2050s(TBE_r_fut_preds_clim_landuse_ssp370)
TBE_r_fut_preds_clim_landuse_ssp585_2050s <- get_decade_layers_2050s(TBE_r_fut_preds_clim_landuse_ssp585)

# WNV
WNV_r_fut_preds_clim_landuse_ssp126_2050s <- get_decade_layers_2050s(WNV_r_fut_preds_clim_landuse_ssp126)
WNV_r_fut_preds_clim_landuse_ssp370_2050s <- get_decade_layers_2050s(WNV_r_fut_preds_clim_landuse_ssp370)
WNV_r_fut_preds_clim_landuse_ssp585_2050s <- get_decade_layers_2050s(WNV_r_fut_preds_clim_landuse_ssp585)


# Write a function that extracts the decade layers for only one specific month
keep_month <- function(r, month_code){
  layer_names <- names(r)
  months <- substr(layer_names, 1, 2)
  r[[months == month_code]]
}

# Apply the month filtering (May for Ixodes ricinus and TBE; July for Culex
# pipiens and WNV)
I_ricinus_r_fut_preds_clim_landuse_ssp126_2050s <- keep_month(I_ricinus_r_fut_preds_clim_landuse_ssp126_2050s, "05") # May
I_ricinus_r_fut_preds_clim_landuse_ssp370_2050s <- keep_month(I_ricinus_r_fut_preds_clim_landuse_ssp370_2050s, "05")
I_ricinus_r_fut_preds_clim_landuse_ssp585_2050s <- keep_month(I_ricinus_r_fut_preds_clim_landuse_ssp585_2050s, "05")

TBE_r_fut_preds_clim_landuse_ssp126_2050s <- keep_month(TBE_r_fut_preds_clim_landuse_ssp126_2050s, "05") # May
TBE_r_fut_preds_clim_landuse_ssp370_2050s <- keep_month(TBE_r_fut_preds_clim_landuse_ssp370_2050s, "05")
TBE_r_fut_preds_clim_landuse_ssp585_2050s <- keep_month(TBE_r_fut_preds_clim_landuse_ssp585_2050s, "05")

C_pipiens_r_fut_preds_clim_landuse_ssp126_2050s <- keep_month(C_pipiens_r_fut_preds_clim_landuse_ssp126_2050s, "07") # July
C_pipiens_r_fut_preds_clim_landuse_ssp370_2050s <- keep_month(C_pipiens_r_fut_preds_clim_landuse_ssp370_2050s, "07")
C_pipiens_r_fut_preds_clim_landuse_ssp585_2050s <- keep_month(C_pipiens_r_fut_preds_clim_landuse_ssp585_2050s, "07")

WNV_r_fut_preds_clim_landuse_ssp126_2050s <- keep_month(WNV_r_fut_preds_clim_landuse_ssp126_2050s, "07") # July
WNV_r_fut_preds_clim_landuse_ssp370_2050s <- keep_month(WNV_r_fut_preds_clim_landuse_ssp370_2050s, "07")
WNV_r_fut_preds_clim_landuse_ssp585_2050s <- keep_month(WNV_r_fut_preds_clim_landuse_ssp585_2050s, "07")



# Write a function that applies the threshold for binarising predictions, then 
# averages the predicted values for each cell over the decade, this is only done
# for cells that show at least one presence for the respective month within the
# decade ("ever present" cells)
threshold_filter_mean <- function(r, threshold){
  r_above_thr <- r >= threshold # 0/1 presence for selected month
  ever_present <- max(r_above_thr, na.rm = TRUE) == 1 # Cells with at least one presence
  r_mean <- mean(r, na.rm = TRUE) # Mean prediction across month
  r_mean[!ever_present] <- NA # Mask cells never above threshold
  r_mean
}



# Apply function for the vector and viruses for the three different socio-economic
# scenarios
I_ricinus_mean_thr_2050s_ssp126 <- threshold_filter_mean(I_ricinus_r_fut_preds_clim_landuse_ssp126_2050s, I_ricinus_thresh)
I_ricinus_mean_thr_2050s_ssp370 <- threshold_filter_mean(I_ricinus_r_fut_preds_clim_landuse_ssp370_2050s, I_ricinus_thresh)
I_ricinus_mean_thr_2050s_ssp585 <- threshold_filter_mean(I_ricinus_r_fut_preds_clim_landuse_ssp585_2050s, I_ricinus_thresh)

C_pipiens_mean_thr_2050s_ssp126 <- threshold_filter_mean(C_pipiens_r_fut_preds_clim_landuse_ssp126_2050s, C_pipiens_thresh)
C_pipiens_mean_thr_2050s_ssp370 <- threshold_filter_mean(C_pipiens_r_fut_preds_clim_landuse_ssp370_2050s, C_pipiens_thresh)
C_pipiens_mean_thr_2050s_ssp585 <- threshold_filter_mean(C_pipiens_r_fut_preds_clim_landuse_ssp585_2050s, C_pipiens_thresh)

TBE_mean_thr_2050s_ssp126 <- threshold_filter_mean(TBE_r_fut_preds_clim_landuse_ssp126_2050s, TBE_thresh)
TBE_mean_thr_2050s_ssp370 <- threshold_filter_mean(TBE_r_fut_preds_clim_landuse_ssp370_2050s, TBE_thresh)
TBE_mean_thr_2050s_ssp585 <- threshold_filter_mean(TBE_r_fut_preds_clim_landuse_ssp585_2050s, TBE_thresh)

WNV_mean_thr_2050s_ssp126 <- threshold_filter_mean(WNV_r_fut_preds_clim_landuse_ssp126_2050s, WNV_thresh)
WNV_mean_thr_2050s_ssp370 <- threshold_filter_mean(WNV_r_fut_preds_clim_landuse_ssp370_2050s, WNV_thresh)
WNV_mean_thr_2050s_ssp585 <- threshold_filter_mean(WNV_r_fut_preds_clim_landuse_ssp585_2050s, WNV_thresh)


# Convert rasters to data frames for plotting
I_ricinus_df_2050s_ssp126 <- as.data.frame(I_ricinus_mean_thr_2050s_ssp126, xy = TRUE)
I_ricinus_df_2050s_ssp370 <- as.data.frame(I_ricinus_mean_thr_2050s_ssp370, xy = TRUE)
I_ricinus_df_2050s_ssp585 <- as.data.frame(I_ricinus_mean_thr_2050s_ssp585, xy = TRUE)

C_pipiens_df_2050s_ssp126 <- as.data.frame(C_pipiens_mean_thr_2050s_ssp126, xy = TRUE)
C_pipiens_df_2050s_ssp370 <- as.data.frame(C_pipiens_mean_thr_2050s_ssp370, xy = TRUE)
C_pipiens_df_2050s_ssp585 <- as.data.frame(C_pipiens_mean_thr_2050s_ssp585, xy = TRUE)
                                           
TBE_df_2050s_ssp126 <- as.data.frame(TBE_mean_thr_2050s_ssp126, xy = TRUE)                                         
TBE_df_2050s_ssp370 <- as.data.frame(TBE_mean_thr_2050s_ssp370, xy = TRUE)                                           
TBE_df_2050s_ssp585 <- as.data.frame(TBE_mean_thr_2050s_ssp585, xy = TRUE)  

WNV_df_2050s_ssp126 <- as.data.frame(WNV_mean_thr_2050s_ssp126, xy = TRUE) 
WNV_df_2050s_ssp370 <- as.data.frame(WNV_mean_thr_2050s_ssp370, xy = TRUE) 
WNV_df_2050s_ssp585 <- as.data.frame(WNV_mean_thr_2050s_ssp585, xy = TRUE) 



# Add a column indicating the studied decade + socio-economic scneario, the 
# studied vector/virus, the vector/virus combination (as we plot them in one figure)
I_ricinus_df_2050s_ssp126$time <- "2050s; ssp126"
I_ricinus_df_2050s_ssp126$species <- "Ixodes ricinus"
I_ricinus_df_2050s_ssp126$combination <- "Ixodes ricinus & TBE"

I_ricinus_df_2050s_ssp370$time <- "2050s; ssp370"
I_ricinus_df_2050s_ssp370$species <- "Ixodes ricinus"
I_ricinus_df_2050s_ssp370$combination <- "Ixodes ricinus & TBE"

I_ricinus_df_2050s_ssp585$time <- "2050s; ssp585"
I_ricinus_df_2050s_ssp585$species <- "Ixodes ricinus"
I_ricinus_df_2050s_ssp585$combination <- "Ixodes ricinus & TBE"


C_pipiens_df_2050s_ssp126$time <- "2050s; ssp126"
C_pipiens_df_2050s_ssp126$species <- "Culex pipiens"
C_pipiens_df_2050s_ssp126$combination <- "Culex pipiens & WNV"

C_pipiens_df_2050s_ssp370$time <- "2050s; ssp370"
C_pipiens_df_2050s_ssp370$species <- "Culex pipiens"
C_pipiens_df_2050s_ssp370$combination <- "Culex pipiens & WNV"

C_pipiens_df_2050s_ssp585$time <- "2050s; ssp585"
C_pipiens_df_2050s_ssp585$species <- "Culex pipiens"
C_pipiens_df_2050s_ssp585$combination <- "Culex pipiens & WNV"

TBE_df_2050s_ssp126$time <- "2050s; ssp126"
TBE_df_2050s_ssp126$species <- "TBE"
TBE_df_2050s_ssp126$combination <- "Ixodes ricinus & TBE"

TBE_df_2050s_ssp370$time <- "2050s; ssp370"
TBE_df_2050s_ssp370$species <- "TBE"
TBE_df_2050s_ssp370$combination <- "Ixodes ricinus & TBE"

TBE_df_2050s_ssp585$time <- "2050s; ssp585"
TBE_df_2050s_ssp585$species <- "TBE"
TBE_df_2050s_ssp585$combination <- "Ixodes ricinus & TBE"

WNV_df_2050s_ssp126$time <- "2050s; ssp126"
WNV_df_2050s_ssp126$species <- "WNV"
WNV_df_2050s_ssp126$combination <- "Culex pipiens & WNV"

WNV_df_2050s_ssp370$time <- "2050s; ssp370"
WNV_df_2050s_ssp370$species <- "WNV"
WNV_df_2050s_ssp370$combination <- "Culex pipiens & WNV"

WNV_df_2050s_ssp585$time <- "2050s; ssp585"
WNV_df_2050s_ssp585$species <- "WNV"
WNV_df_2050s_ssp585$combination <- "Culex pipiens & WNV"


# Combine the prepared data frames into one for each socio-economic scenario
aggregated_df_fut_ssp126 <- bind_rows(I_ricinus_df_2050s_ssp126, C_pipiens_df_2050s_ssp126,
                                      TBE_df_2050s_ssp126, WNV_df_2050s_ssp126)

aggregated_df_fut_ssp370 <- bind_rows(I_ricinus_df_2050s_ssp370, C_pipiens_df_2050s_ssp370,
                                      TBE_df_2050s_ssp370, WNV_df_2050s_ssp370)

aggregated_df_fut_ssp585 <- bind_rows(I_ricinus_df_2050s_ssp585, C_pipiens_df_2050s_ssp585,
                                      TBE_df_2050s_ssp585, WNV_df_2050s_ssp585)



# Save the data frames
save(aggregated_df_fut_ssp126, file = "output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp126.RData")
save(aggregated_df_fut_ssp370, file = "output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp370.RData")
save(aggregated_df_fut_ssp585, file = "output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp585.RData")




#-------------------------------------------------------------------------------

# 3. Visualise distribution trends  --------------------------------------------
# for main vectors as well as viruses

# a) Prepare data frame for visualisation --------------------------------------

# Load the needed data
eu_eea_mask <- terra::rast("input_data/spatial_data/eu_eea_mask.tif") # EU/EEA mask at a 0.5° resolution
eu_eea_mask_df <- as.data.frame(eu_eea_mask, xy = TRUE, na.rm = TRUE) # Turn mask into data frame (for plotting)

load("output_data/results/distribution_trends/distribution_trends_past_vector_virus.RData") # Historical target decades
load("output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp126.RData") # Future target decade - ssp126
load("output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp370.RData") # Future target decade - ssp370
load("output_data/results/distribution_trends/distribution_trends_fut_vector_virus_ssp585.RData") # Future target decade - ssp585

# Separate vector from virus data as these require different legends in the plot
distribution_trends_past_vector <- aggregated_df_past[aggregated_df_past$species %in% c("Ixodes ricinus", "Culex pipiens"), ]
distribution_trends_past_virus <- aggregated_df_past[aggregated_df_past$species %in% c("TBE", "WNV"), ]

distribution_trends_fut_vector_ssp126 <- aggregated_df_fut_ssp126[aggregated_df_fut_ssp126$species %in% c("Ixodes ricinus", "Culex pipiens"), ]
distribution_trends_fut_virus_ssp126 <- aggregated_df_fut_ssp126[aggregated_df_fut_ssp126$species %in% c("TBE", "WNV"), ]

distribution_trends_fut_vector_ssp370 <- aggregated_df_fut_ssp370[aggregated_df_fut_ssp370$species %in% c("Ixodes ricinus", "Culex pipiens"), ]
distribution_trends_fut_virus_ssp370 <- aggregated_df_fut_ssp370[aggregated_df_fut_ssp370$species %in% c("TBE", "WNV"), ]

distribution_trends_fut_vector_ssp585 <- aggregated_df_fut_ssp585[aggregated_df_fut_ssp585$species %in% c("Ixodes ricinus", "Culex pipiens"), ]
distribution_trends_fut_virus_ssp585 <- aggregated_df_fut_ssp585[aggregated_df_fut_ssp585$species %in% c("TBE", "WNV"), ]

# Aggregate data frame for plotting
distribution_trends_vector_ssp126 <- bind_rows(distribution_trends_past_vector, distribution_trends_fut_vector_ssp126)
distribution_trends_virus_ssp126 <- bind_rows(distribution_trends_past_virus, distribution_trends_fut_virus_ssp126)

distribution_trends_vector_ssp370 <- bind_rows(distribution_trends_past_vector, distribution_trends_fut_vector_ssp370)
distribution_trends_virus_ssp370 <- bind_rows(distribution_trends_past_virus, distribution_trends_fut_virus_ssp370)

distribution_trends_vector_ssp585 <- bind_rows(distribution_trends_past_vector, distribution_trends_fut_vector_ssp585)
distribution_trends_virus_ssp585 <- bind_rows(distribution_trends_past_virus, distribution_trends_fut_virus_ssp585)

# Make sure the vector/virus combination appears in the correct order
distribution_trends_vector_ssp126$combination <- factor(distribution_trends_vector_ssp126$combination, 
                                                        levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))

distribution_trends_virus_ssp126$combination <- factor(distribution_trends_virus_ssp126$combination, 
                                                       levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))

distribution_trends_vector_ssp370$combination <- factor(distribution_trends_vector_ssp370$combination, 
                                                        levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))

distribution_trends_virus_ssp370$combination <- factor(distribution_trends_virus_ssp370$combination, 
                                                        levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))

distribution_trends_vector_ssp585$combination <- factor(distribution_trends_vector_ssp585$combination, 
                                                        levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))

distribution_trends_virus_ssp585$combination <- factor(distribution_trends_virus_ssp585$combination, 
                                                       levels = c("Ixodes ricinus & TBE", "Culex pipiens & WNV"))



# Visualise vector (greyish scale) and virus (colour scale) distribution 
# on top of each other (cells are shown that were at least one time above
# the maxTSS threshold for the target month within the target decade);
# for three different decades

# b) Visualise historical and future distribution based on ssp126 --------------
ggplot() +
  geom_raster(data = eu_eea_mask_df, aes(x = x, y = y), fill = "skyblue3", alpha = 0.3) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_vector_ssp126, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradient(name = "Vector suitability", low = "grey80", high = "grey10",
                      limits = c(0, 1)) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_virus_ssp126,
              aes(x = x, y = y, fill = mean)) +
  scale_fill_viridis_c(option = "plasma", name = "Virus suitability",
                       limits = c(0, 1)) +
  facet_grid2(time ~ combination,
              strip = strip_themed(background_x = elem_list_rect(fill = c("steelblue3", "lightsteelblue1")),
                                   text_x = elem_list_text(face = c("bold.italic", "bold.italic")))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11.5),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey75", colour = NA),
    strip.text = element_text(size = 15, face = "bold")
  )

# Save the figure
ggsave("output_data/plots/distribution_trends/distribution_trends_Ixodes_TBE_Culex_WNV_ssp126.png", width = 7.5, height = 8.5)



# c) Visualise historical and future distribution based on ssp370 --------------
ggplot() +
  geom_raster(data = eu_eea_mask_df, aes(x = x, y = y), fill = "skyblue3", alpha = 0.3) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_vector_ssp370, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradient(name = "Vector suitability", low = "grey80", high = "grey10",
                      limits = c(0, 1)) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_virus_ssp370,
              aes(x = x, y = y, fill = mean)) +
  scale_fill_viridis_c(option = "plasma", name = "Virus suitability",
                       limits = c(0, 1)) +
  facet_grid2(time ~ combination,
              strip = strip_themed(background_x = elem_list_rect(fill = c("steelblue3", "lightsteelblue1")),
                                   text_x = elem_list_text(face = c("bold.italic", "bold.italic")))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11.5),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey75", colour = NA),
    strip.text = element_text(size = 15, face = "bold")
  )
  

# Save the figure
ggsave("output_data/plots/distribution_trends/distribution_trends_Ixodes_TBE_Culex_WNV_ssp370.png", width = 7.5, height = 8.5)


# d) Visualise historical and future distribution based on ssp585 --------------

ggplot() +
  geom_raster(data = eu_eea_mask_df, aes(x = x, y = y), fill = "skyblue3", alpha = 0.3) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_vector_ssp585, aes(x = x, y = y, fill = mean)) +
  scale_fill_gradient(name = "Vector suitability", low = "grey80", high = "grey10",
                      limits = c(0, 1)) +
  new_scale_fill() +
  geom_raster(data = distribution_trends_virus_ssp585,
              aes(x = x, y = y, fill = mean)) +
  scale_fill_viridis_c(option = "plasma", name = "Virus suitability",
                       limits = c(0, 1)) +
  facet_grid2(time ~ combination,
              strip = strip_themed(background_x = elem_list_rect(fill = c("steelblue3", "lightsteelblue1")),
                                   text_x = elem_list_text(face = c("bold.italic", "bold.italic")))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11.5),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9),
    strip.background = element_rect(fill = "grey75", colour = NA),
    strip.text = element_text(size = 15, face = "bold")
  )

# Save the figure
ggsave("output_data/plots/distribution_trends/distribution_trends_Ixodes_TBE_Culex_WNV_ssp585.png", width = 7.5, height = 8.5)



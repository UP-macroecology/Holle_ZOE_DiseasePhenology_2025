# ZOE project 
# Disease phenology analysis of Culex pipiens in Europe (primary transmitter of WNV)

# ---------------------------------------------------------------------- #
#                          06b. Model prediction                         #
# ---------------------------------------------------------------------- #


# Load needed packages
library(mgcv)
library(maxnet)
library(randomForest)
library(gbm)
library(dismo)
library(terra)
library(tidyverse)


# Load needed data
load("output_data/models/C_pipiens_SDMs.RData") # Load fitted models
load("output_data/validation/C_pipiens_validation.RData") # Load validation results



#-------------------------------------------------------------------------------

# 1. Past monthly predictions from 1970 to 2019 --------------------------------
# under observed climate and land use change
# Based on all four applied algorithms and their ensemble

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
# For all algorithms and their ensemble
r_curr_preds_clim_landuse_ens <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_landuse_ens_bin <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_landuse_glm <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_landuse_gam <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_landuse_rf <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_landuse_brt <- terra::rast(example_data, nlyrs = 600)


for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  # For all algorithms and their ensemble
  r_curr_preds_year_ens <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_ens_bin <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_glm <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_gam <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_rf <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_brt <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Check how many rows the data frame with environmental data would have
    env_df_check <- data.frame(crds(env_data),as.points(env_data))
    
    # Create a matrix to store the predictions of the algorithms
    preds_glm_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_gam_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_rf_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_brt_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    
    
    for (n in 1:length(models_glm)) { # Start of the loop over the number of constructed models with different predictors (using GLM as example)
      
      print(n)
      
      # Extract the predictors within that model (use GLM as example model)
      model_name <- names(models_glm)[n]
      my_preds <- unlist(strsplit(model_name, "\\+"))
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      
      # Insert the predictions in the prepared data frame
      print("GLM")
      preds_glm_month[, n] <- predict(models_glm[[n]], env_df, type='response')
      print("GAM")
      preds_gam_month[, n] <- predict(models_gam[[n]], env_df[,my_preds], type='response')
      print("RF")
      preds_rf_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict(models_rf[[n]][[i]], env_df, type='response')}))
      print("BRT")
      preds_brt_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict.gbm(models_brt[[n]][[i]], env_df, n.trees=models_brt[[n]][[i]]$gbm.call$best.trees, type="response")}))
      
    } # Close the loop over the number of models
    
    # Average the predictions per algorithm and store them with coordinate information
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = rowMeans(preds_glm_month),
                             gam = rowMeans(preds_gam_month),
                             rf = rowMeans(preds_rf_month),
                             brt = rowMeans(preds_brt_month))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Binarise ensemble predictions
    curr_preds$bin_pred = ifelse(curr_preds$mean_prob >= comp_perf[comp_perf$alg == "mean_prob", "thresh"], 1, 0)
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster as well as the rasters based on the different algorithms
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    r_curr_preds_ens_bin <- r_curr_preds[["bin_pred"]]
    r_curr_preds_glm <- r_curr_preds[["glm"]]
    r_curr_preds_gam <- r_curr_preds[["gam"]]
    r_curr_preds_rf <- r_curr_preds[["rf"]]
    r_curr_preds_brt <- r_curr_preds[["brt"]]
    
    # Save the rasters of the different algorithms and their ensemble in the prepared rasterstack
    r_curr_preds_year_ens <- c(r_curr_preds_year_ens, r_curr_preds_ens)
    r_curr_preds_year_ens_bin <- c(r_curr_preds_year_ens_bin, r_curr_preds_ens_bin)
    r_curr_preds_year_glm <- c(r_curr_preds_year_glm, r_curr_preds_glm)
    r_curr_preds_year_gam <- c(r_curr_preds_year_gam, r_curr_preds_gam)
    r_curr_preds_year_rf <- c(r_curr_preds_year_rf, r_curr_preds_rf)
    r_curr_preds_year_brt <- c(r_curr_preds_year_brt, r_curr_preds_brt)
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year_ens) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_ens_bin) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_glm) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_gam) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_rf) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_brt) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_clim_landuse_ens <- c(r_curr_preds_clim_landuse_ens, r_curr_preds_year_ens)
  r_curr_preds_clim_landuse_ens_bin <- c(r_curr_preds_clim_landuse_ens_bin, r_curr_preds_year_ens_bin)
  r_curr_preds_clim_landuse_glm <- c(r_curr_preds_clim_landuse_glm, r_curr_preds_year_glm)
  r_curr_preds_clim_landuse_gam <- c(r_curr_preds_clim_landuse_gam, r_curr_preds_year_gam)
  r_curr_preds_clim_landuse_rf <- c(r_curr_preds_clim_landuse_rf, r_curr_preds_year_rf)
  r_curr_preds_clim_landuse_brt <- c(r_curr_preds_clim_landuse_brt, r_curr_preds_year_brt)
  
  
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_clim_landuse_ens, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_landuse_ens_bin, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_landuse_glm, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_glm_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_landuse_gam, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_gam_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_landuse_rf, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_rf_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_landuse_brt, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_brt_1970_2019.tif", overwrite=T)



#-------------------------------------------------------------------------------

# 2. Past monthly predictions from 1970 to 2019 --------------------------------
# under observed land use change and the counterfactual climate scenario
# (no climate change / detrended climate data)
# Based on all four applied algorithms and their ensemble

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
# For all algorithms and their ensemble
r_curr_preds_noclim_landuse_ens <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_landuse_ens_bin <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_landuse_glm <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_landuse_gam <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_landuse_rf <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_landuse_brt <- terra::rast(example_data, nlyrs = 600)

for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  # For all algorithms and their ensemble
  r_curr_preds_year_ens <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_ens_bin <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_glm <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_gam <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_rf <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_brt <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/CounterClim/processed_data/CounterClim_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Check how many rows the data frame with environmental data would have
    env_df_check <- data.frame(crds(env_data),as.points(env_data))
    
    # Create a matrix to store the predictions of the algorithms
    preds_glm_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_gam_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_rf_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_brt_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    
    
    for (n in 1:length(models_glm)) { # Start of the loop over the number of constructed models with different predictors (using GLM as example)
      
      print(n)
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Extract the predictors within that model (use GLM as example model)
      model_name <- names(models_glm)[n]
      my_preds <- unlist(strsplit(model_name, "\\+"))
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      
      # Insert the predictions in the prepared data frame
      print("GLM")
      preds_glm_month[, n] <- predict(models_glm[[n]], env_df, type='response')
      print("GAM")
      preds_gam_month[, n] <- predict(models_gam[[n]], env_df[,my_preds], type='response')
      print("RF")
      preds_rf_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict(models_rf[[n]][[i]], env_df, type='response')}))
      print("BRT")
      preds_brt_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict.gbm(models_brt[[n]][[i]], env_df, n.trees=models_brt[[n]][[i]]$gbm.call$best.trees, type="response")}))
      
    } # Close the loop over the number of models
    
    
    # Average the predicitons per algorithm and store them with coordinate information
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = rowMeans(preds_glm_month),
                             gam = rowMeans(preds_gam_month),
                             rf = rowMeans(preds_rf_month),
                             brt = rowMeans(preds_brt_month))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Binarise ensemble predictions
    curr_preds$bin_pred = ifelse(curr_preds$mean_prob >= comp_perf[comp_perf$alg == "mean_prob", "thresh"], 1, 0)
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster as well as the rasters based on the different algorithms
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    r_curr_preds_ens_bin <- r_curr_preds[["bin_pred"]]
    r_curr_preds_glm <- r_curr_preds[["glm"]]
    r_curr_preds_gam <- r_curr_preds[["gam"]]
    r_curr_preds_rf <- r_curr_preds[["rf"]]
    r_curr_preds_brt <- r_curr_preds[["brt"]]
    
    # Save the rasters of the different algorithms and their ensemble in the prepared rasterstack
    r_curr_preds_year_ens <- c(r_curr_preds_year_ens, r_curr_preds_ens)
    r_curr_preds_year_ens_bin <- c(r_curr_preds_year_ens_bin, r_curr_preds_ens_bin)
    r_curr_preds_year_glm <- c(r_curr_preds_year_glm, r_curr_preds_glm)
    r_curr_preds_year_gam <- c(r_curr_preds_year_gam, r_curr_preds_gam)
    r_curr_preds_year_rf <- c(r_curr_preds_year_rf, r_curr_preds_rf)
    r_curr_preds_year_brt <- c(r_curr_preds_year_brt, r_curr_preds_brt)
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year_ens) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_ens_bin) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_glm) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_gam) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_rf) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_brt) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_noclim_landuse_ens <- c(r_curr_preds_noclim_landuse_ens, r_curr_preds_year_ens)
  r_curr_preds_noclim_landuse_ens_bin <- c(r_curr_preds_noclim_landuse_ens_bin, r_curr_preds_year_ens_bin)
  r_curr_preds_noclim_landuse_glm <- c(r_curr_preds_noclim_landuse_glm, r_curr_preds_year_glm)
  r_curr_preds_noclim_landuse_gam <- c(r_curr_preds_noclim_landuse_gam, r_curr_preds_year_gam)
  r_curr_preds_noclim_landuse_rf <- c(r_curr_preds_noclim_landuse_rf, r_curr_preds_year_rf)
  r_curr_preds_noclim_landuse_brt <- c(r_curr_preds_noclim_landuse_brt, r_curr_preds_year_brt)
  
  
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_noclim_landuse_ens, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_landuse_ens_bin, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_landuse_glm, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_glm_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_landuse_gam, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_gam_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_landuse_rf, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_rf_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_landuse_brt, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_brt_1970_2019.tif", overwrite=T)




#-------------------------------------------------------------------------------

# 3. Past monthly predictions from 1970 to 2019 --------------------------------
# under counterfactual land use scenario and the factual climate data
# (land use reference year 1901)
# Based on all four applied algorithms and their ensemble

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
# For all algorithms and their ensemble
r_curr_preds_clim_nolanduse_ens <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_nolanduse_ens_bin <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_nolanduse_glm <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_nolanduse_gam <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_nolanduse_rf <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_clim_nolanduse_brt <- terra::rast(example_data, nlyrs = 600)


for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  # For all algorithms and their ensemble
  r_curr_preds_year_ens <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_ens_bin <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_glm <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_gam <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_rf <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_brt <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/CounterLandUse/processed_data/CounterLandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Check how many rows the data frame with environmental data would have
    env_df_check <- data.frame(crds(env_data),as.points(env_data))
    
    # Create a matrix to store the predictions of the algorithms
    preds_glm_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_gam_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_rf_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_brt_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    
    for (n in 1:length(models_glm)) { # Start of the loop over the number of constructed models with different predictors (using GLM as example)
      
      print(n)
      
      # Extract the predictors within that model (use GLM as example model)
      model_name <- names(models_glm)[n]
      my_preds <- unlist(strsplit(model_name, "\\+"))
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      
      # Insert the predictions in the prepared data frame
      print("GLM")
      preds_glm_month[, n] <- predict(models_glm[[n]], env_df, type='response')
      print("GAM")
      preds_gam_month[, n] <- predict(models_gam[[n]], env_df[,my_preds], type='response')
      print("RF")
      preds_rf_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict(models_rf[[n]][[i]], env_df, type='response')}))
      print("BRT")
      preds_brt_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict.gbm(models_brt[[n]][[i]], env_df, n.trees=models_brt[[n]][[i]]$gbm.call$best.trees, type="response")}))
      
    } # Close the loop over the number of models
    
    
    # Average the predicitons per algorithm and store them with coordinate information
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = rowMeans(preds_glm_month),
                             gam = rowMeans(preds_gam_month),
                             rf = rowMeans(preds_rf_month),
                             brt = rowMeans(preds_brt_month))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Binarise ensemble predictions
    curr_preds$bin_pred = ifelse(curr_preds$mean_prob >= comp_perf[comp_perf$alg == "mean_prob", "thresh"], 1, 0)
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster as well as the rasters based on the different algorithms
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    r_curr_preds_ens_bin <- r_curr_preds[["bin_pred"]]
    r_curr_preds_glm <- r_curr_preds[["glm"]]
    r_curr_preds_gam <- r_curr_preds[["gam"]]
    r_curr_preds_rf <- r_curr_preds[["rf"]]
    r_curr_preds_brt <- r_curr_preds[["brt"]]
    
    # Save the rasters of the different algorithms and their ensemble in the prepared rasterstack
    r_curr_preds_year_ens <- c(r_curr_preds_year_ens, r_curr_preds_ens)
    r_curr_preds_year_ens_bin <- c(r_curr_preds_year_ens_bin, r_curr_preds_ens_bin)
    r_curr_preds_year_glm <- c(r_curr_preds_year_glm, r_curr_preds_glm)
    r_curr_preds_year_gam <- c(r_curr_preds_year_gam, r_curr_preds_gam)
    r_curr_preds_year_rf <- c(r_curr_preds_year_rf, r_curr_preds_rf)
    r_curr_preds_year_brt <- c(r_curr_preds_year_brt, r_curr_preds_brt)
    
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year_ens) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_ens_bin) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_glm) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_gam) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_rf) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_brt) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_clim_nolanduse_ens <- c(r_curr_preds_clim_nolanduse_ens, r_curr_preds_year_ens)
  r_curr_preds_clim_nolanduse_ens_bin <- c(r_curr_preds_clim_nolanduse_ens_bin, r_curr_preds_year_ens_bin)
  r_curr_preds_clim_nolanduse_glm <- c(r_curr_preds_clim_nolanduse_glm, r_curr_preds_year_glm)
  r_curr_preds_clim_nolanduse_gam <- c(r_curr_preds_clim_nolanduse_gam, r_curr_preds_year_gam)
  r_curr_preds_clim_nolanduse_rf <- c(r_curr_preds_clim_nolanduse_rf, r_curr_preds_year_rf)
  r_curr_preds_clim_nolanduse_brt <- c(r_curr_preds_clim_nolanduse_brt, r_curr_preds_year_brt)
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_clim_nolanduse_ens, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_nolanduse_ens_bin, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_nolanduse_glm, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_glm_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_nolanduse_gam, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_gam_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_nolanduse_rf, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_rf_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_clim_nolanduse_brt, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_brt_1970_2019.tif", overwrite=T)




#-------------------------------------------------------------------------------

# 4. Past monthly predictions from 1970 to 2019 --------------------------------
# under counterfactual land use scenario and the counterfactual climate scenario
# Based on all four applied algorithms and their ensemble

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
# For all algorithms and their ensemble
r_curr_preds_noclim_nolanduse_ens <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_nolanduse_ens_bin <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_nolanduse_glm <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_nolanduse_gam <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_nolanduse_rf <- terra::rast(example_data, nlyrs = 600)
r_curr_preds_noclim_nolanduse_brt <- terra::rast(example_data, nlyrs = 600)


for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  # For all algorithms and their ensemble
  r_curr_preds_year_ens <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_ens_bin <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_glm <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_gam <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_rf <- terra::rast(example_data, nlyrs = 12)
  r_curr_preds_year_brt <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/CounterClim/processed_data/CounterClim_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/CounterLandUse/processed_data/CounterLandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Check how many rows the data frame with environmental data would have
    env_df_check <- data.frame(crds(env_data),as.points(env_data))
    
    # Create a matrix to store the predictions of the algorithms
    preds_glm_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_gam_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_rf_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    preds_brt_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
    
    for (n in 1:length(models_glm)) { # Start of the loop over the number of constructed models with different predictors (using GLM as example)
      
      print(n)
      
      # Extract the predictors within that model (use GLM as example model)
      model_name <- names(models_glm)[n]
      my_preds <- unlist(strsplit(model_name, "\\+"))
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      
      # Insert the predictions in the prepared data frame
      print("GLM")
      preds_glm_month[, n] <- predict(models_glm[[n]], env_df, type='response')
      print("GAM")
      preds_gam_month[, n] <- predict(models_gam[[n]], env_df[,my_preds], type='response')
      print("RF")
      preds_rf_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict(models_rf[[n]][[i]], env_df, type='response')}))
      print("BRT")
      preds_brt_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict.gbm(models_brt[[n]][[i]], env_df, n.trees=models_brt[[n]][[i]]$gbm.call$best.trees, type="response")}))
      
    } # Close the loop over the number of models
    
    
    # Average the predicitons per algorithm and store them with coordinate information
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = rowMeans(preds_glm_month),
                             gam = rowMeans(preds_gam_month),
                             rf = rowMeans(preds_rf_month),
                             brt = rowMeans(preds_brt_month))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Binarise ensemble predictions
    curr_preds$bin_pred = ifelse(curr_preds$mean_prob >= comp_perf[comp_perf$alg == "mean_prob", "thresh"], 1, 0)
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster as well as the rasters based on the different algorithms
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    r_curr_preds_ens_bin <- r_curr_preds[["bin_pred"]]
    r_curr_preds_glm <- r_curr_preds[["glm"]]
    r_curr_preds_gam <- r_curr_preds[["gam"]]
    r_curr_preds_rf <- r_curr_preds[["rf"]]
    r_curr_preds_brt <- r_curr_preds[["brt"]]
    
    # Save the rasters of the different algorithms and their ensemble in the prepared rasterstack
    r_curr_preds_year_ens <- c(r_curr_preds_year_ens, r_curr_preds_ens)
    r_curr_preds_year_ens_bin <- c(r_curr_preds_year_ens_bin, r_curr_preds_ens_bin)
    r_curr_preds_year_glm <- c(r_curr_preds_year_glm, r_curr_preds_glm)
    r_curr_preds_year_gam <- c(r_curr_preds_year_gam, r_curr_preds_gam)
    r_curr_preds_year_rf <- c(r_curr_preds_year_rf, r_curr_preds_rf)
    r_curr_preds_year_brt <- c(r_curr_preds_year_brt, r_curr_preds_brt)
    
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year_ens) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_ens_bin) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_glm) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_gam) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_rf) <- sprintf("%02d/%s", 1:12, y)
  names(r_curr_preds_year_brt) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_noclim_nolanduse_ens <- c(r_curr_preds_noclim_nolanduse_ens, r_curr_preds_year_ens)
  r_curr_preds_noclim_nolanduse_ens_bin <- c(r_curr_preds_noclim_nolanduse_ens_bin, r_curr_preds_year_ens_bin)
  r_curr_preds_noclim_nolanduse_glm <- c(r_curr_preds_noclim_nolanduse_glm, r_curr_preds_year_glm)
  r_curr_preds_noclim_nolanduse_gam <- c(r_curr_preds_noclim_nolanduse_gam, r_curr_preds_year_gam)
  r_curr_preds_noclim_nolanduse_rf <- c(r_curr_preds_noclim_nolanduse_rf, r_curr_preds_year_rf)
  r_curr_preds_noclim_nolanduse_brt <- c(r_curr_preds_noclim_nolanduse_brt, r_curr_preds_year_brt)
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_noclim_nolanduse_ens, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_nolanduse_ens_bin, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_nolanduse_glm, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_glm_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_nolanduse_gam, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_gam_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_nolanduse_rf, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_rf_1970_2019.tif", overwrite=T)
terra::writeRaster(r_curr_preds_noclim_nolanduse_brt, filename = "output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_brt_1970_2019.tif", overwrite=T)





#-------------------------------------------------------------------------------

# 5. Future monthly predictions from 2030 to 2070 ------------------------------
# under climate change and land use change 
# (3 different scenarios, 5 different climate models)
# Based on all four applied algorithms and their ensemble

# Prepare a vector containing the years for monthly predictions
years <- c(2030:2070)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Create path to directory containing environmental data
datapath_env_fut <- file.path("input_data/environmental_data/ISIMIP3b/")
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Prepare a vector containing the three different forcing scenarios
scenario <- c("ssp126", "ssp370", "ssp585")

# Create a vector containing the different climate models
clim_models <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")


for (s in scenario) { # Start of the loop over the three different forcing scenarios
  
  print(s)
  
  for (l in clim_models) { # Start of the loop over the different climate models
    
    print(l)
    
    # Check if file of ensemble results already exist
    file_exists <- file.exists(paste0("output_data/results/C_pipiens_preds_clim_landuse_ens_2030_2070_",l,"_",s,".tif"))
    
    # If that is the case. skip to the next iteration
    if (file_exists == TRUE) { print("prediction already done")
      next
    }
    
    
    # Load a raster as example template 
    example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
    
    # Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
    # For all algorithms and their ensemble
    r_fut_preds_clim_landuse_ens <- terra::rast(example_data, nlyrs = 492)
    r_fut_preds_clim_landuse_ens_bin <- terra::rast(example_data, nlyrs = 492)
    r_fut_preds_clim_landuse_glm <- terra::rast(example_data, nlyrs = 492)
    r_fut_preds_clim_landuse_gam <- terra::rast(example_data, nlyrs = 492)
    r_fut_preds_clim_landuse_rf <- terra::rast(example_data, nlyrs = 492)
    r_fut_preds_clim_landuse_brt <- terra::rast(example_data, nlyrs = 492)
    
    for (y in years) { # Start of the loop over the prediction years
      
      print(y)
      
      # Load a raster as example template 
      example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
      
      # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
      # For all algorithms and their ensemble
      r_fut_preds_year_ens <- terra::rast(example_data, nlyrs = 12)
      r_fut_preds_year_ens_bin <- terra::rast(example_data, nlyrs = 12)
      r_fut_preds_year_glm <- terra::rast(example_data, nlyrs = 12)
      r_fut_preds_year_gam <- terra::rast(example_data, nlyrs = 12)
      r_fut_preds_year_rf <- terra::rast(example_data, nlyrs = 12)
      r_fut_preds_year_brt <- terra::rast(example_data, nlyrs = 12)
      
      for (m in month) { # Start of the loop over all months of a year
        
        print(m)
        
        # Load the environmental data for the specific year and month
        Climate_data <- terra::rast(paste0(datapath_env_fut, "/Climate/",s,"/processed_data/",l,"/Climate_future_data_",m,"_",y,"_",s,".tif")) # Climate data
        LandUse_data <- terra::rast(paste0(datapath_env_fut, "/LandUse/",s,"/processed_data/LandUse_future_data_",y,"_",s,".tif")) # Land cover data
        
        # Stack the environmental data
        env_data <- c(Climate_data, LandUse_data)
        
        # Check how many rows the data frame with environmental data would have
        env_df_check <- data.frame(crds(env_data),as.points(env_data))
        
        # Create a matrix to store the predictions of the algorithms
        preds_glm_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
        preds_gam_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
        preds_rf_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
        preds_brt_month <- matrix(nrow = nrow(env_df_check), ncol = 3)
        
        
        for (n in 1:length(models_glm)) { # Start of the loop over the number of constructed models with different predictors (using GLM as example)
          
          print(n)
          
          # Extract the predictors within that model (use GLM as example model)
          model_name <- names(models_glm)[n]
          my_preds <- unlist(strsplit(model_name, "\\+"))
          
          # Prepare a data frame with environmental data
          env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
          
          # Make predictions of all models
          print("start of model predictions")
          
          # Insert the predictions in the prepared data frame
          print("GLM")
          preds_glm_month[, n] <- predict(models_glm[[n]], env_df, type='response')
          print("GAM")
          preds_gam_month[, n] <- predict(models_gam[[n]], env_df[,my_preds], type='response')
          print("RF")
          preds_rf_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict(models_rf[[n]][[i]], env_df, type='response')}))
          print("BRT")
          preds_brt_month[, n] <- rowMeans(sapply(1:background_presence_ratio, FUN=function(i){print(i); predict.gbm(models_brt[[n]][[i]], env_df, n.trees=models_brt[[n]][[i]]$gbm.call$best.trees, type="response")}))
          
        } # Close the loop over the number of models
        
        
        # Average the predicitons per algorithm and store them with coordinate information
        fut_preds <- data.frame(env_df[,1:2], 
                                glm = rowMeans(preds_glm_month),
                                gam = rowMeans(preds_gam_month),
                                rf = rowMeans(preds_rf_month),
                                brt = rowMeans(preds_brt_month))
        
        # Make ensemble predictions
        fut_preds$mean_prob = rowMeans(fut_preds[,-c(1:2)])
        
        # Binarise ensemble predictions
        fut_preds$bin_pred = ifelse(fut_preds$mean_prob >= comp_perf[comp_perf$alg == "mean_prob", "thresh"], 1, 0)
        
        # Make Spatrasters from predictions
        r_fut_preds <- terra::rast(fut_preds, crs = crs(env_data))
        
        # Extract the ensemble raster as well as the rasters based on the different algorithms
        r_fut_preds_ens <- r_fut_preds[["mean_prob"]]
        r_fut_preds_ens_bin <- r_fut_preds[["bin_pred"]]
        r_fut_preds_glm <- r_fut_preds[["glm"]]
        r_fut_preds_gam <- r_fut_preds[["gam"]]
        r_fut_preds_rf <- r_fut_preds[["rf"]]
        r_fut_preds_brt <- r_fut_preds[["brt"]]
        
        # Save the rasters of the different algorithms and their ensemble in the prepared rasterstack
        r_fut_preds_year_ens <- c(r_fut_preds_year_ens, r_fut_preds_ens)
        r_fut_preds_year_ens_bin <- c(r_fut_preds_year_ens_bin, r_fut_preds_ens_bin)
        r_fut_preds_year_glm <- c(r_fut_preds_year_glm, r_fut_preds_glm)
        r_fut_preds_year_gam <- c(r_fut_preds_year_gam, r_fut_preds_gam)
        r_fut_preds_year_rf <- c(r_fut_preds_year_rf, r_fut_preds_rf)
        r_fut_preds_year_brt <- c(r_fut_preds_year_brt, r_fut_preds_brt)
        
      } # End of loop over all months
      
      
      # Make sure that raster names are correct
      names(r_fut_preds_year_ens) <- sprintf("%02d/%s", 1:12, y)
      names(r_fut_preds_year_ens_bin) <- sprintf("%02d/%s", 1:12, y)
      names(r_fut_preds_year_glm) <- sprintf("%02d/%s", 1:12, y)
      names(r_fut_preds_year_gam) <- sprintf("%02d/%s", 1:12, y)
      names(r_fut_preds_year_rf) <- sprintf("%02d/%s", 1:12, y)
      names(r_fut_preds_year_brt) <- sprintf("%02d/%s", 1:12, y)
      
      # Stack the all raster for each year
      r_fut_preds_clim_landuse_ens <- c(r_fut_preds_clim_landuse_ens, r_fut_preds_year_ens)
      r_fut_preds_clim_landuse_ens_bin <- c(r_fut_preds_clim_landuse_ens_bin, r_fut_preds_year_ens_bin)
      r_fut_preds_clim_landuse_glm <- c(r_fut_preds_clim_landuse_glm, r_fut_preds_year_glm)
      r_fut_preds_clim_landuse_gam <- c(r_fut_preds_clim_landuse_gam, r_fut_preds_year_gam)
      r_fut_preds_clim_landuse_rf <- c(r_fut_preds_clim_landuse_rf, r_fut_preds_year_rf)
      r_fut_preds_clim_landuse_brt <- c(r_fut_preds_clim_landuse_brt, r_fut_preds_year_brt)
      
      
      
    } # End of the loop over all considered years
    
    
    # Save the raster outputs
    terra::writeRaster(r_fut_preds_clim_landuse_ens, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_",l,"_",s,".tif"), overwrite=T)
    terra::writeRaster(r_fut_preds_clim_landuse_ens_bin, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_",l,"_",s,".tif"), overwrite=T)
    terra::writeRaster(r_fut_preds_clim_landuse_glm, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_glm_2030_2070_",l,"_",s,".tif"), overwrite=T)
    terra::writeRaster(r_fut_preds_clim_landuse_gam, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_gam_2030_2070_",l,"_",s,".tif"), overwrite=T)
    terra::writeRaster(r_fut_preds_clim_landuse_rf, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_rf_2030_2070_",l,"_",s,".tif"), overwrite=T)
    terra::writeRaster(r_fut_preds_clim_landuse_brt, filename = paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_brt_2030_2070_",l,"_",s,".tif"), overwrite=T)
    
    
    
  } # Close the loop over the five climate models
  
  
} # Close the loop over the three forcing scenarios





# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

# ---------------------------------------------------------------------- #
#                          06a. Model prediction                         #
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
load("output_data/models/I_ricinus_SDMs.RData") # Load fitted models
load("output_data/validation/I_ricinus_validation.RData") # Load validation results



#-------------------------------------------------------------------------------

# 1. Past monthly predictions from 1970 to 2019 --------------------------------
# under observed climate and land use change

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
r_curr_preds_clim_landuse <- terra::rast(example_data, nlyrs = 600)


for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  r_curr_preds_year <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Prepare a data frame with environmental data
    env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
    
    # Make predictions of all models
    print("start of model predictions")
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = predict(m_glm, env_df, type='response'),
                             gam = predict(m_gam, env_df[,my_preds], type='response'),
                             rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                             brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    
    # Save ensemble raster in prepared rasterstack
    r_curr_preds_year <- c(r_curr_preds_year, r_curr_preds_ens)
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_clim_landuse <- c(r_curr_preds_clim_landuse, r_curr_preds_year)
  
  
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_clim_landuse, filename = "output_data/results/I_ricinus_preds_clim_landuse_1970_2019.tif", overwrite=T)



#-------------------------------------------------------------------------------

# 2. Past monthly predictions from 1970 to 2019 --------------------------------
# under observed land use change and the counterfactual climate scenario
# (no climate change / detrended climate data)

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
r_curr_preds_noclim_landuse <- terra::rast(example_data, nlyrs = 600)

for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  r_curr_preds_year <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/CounterClim/processed_data/CounterClim_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Prepare a data frame with environmental data
    env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
    
    # Make predictions of all models
    print("start of model predictions")
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = predict(m_glm, env_df, type='response'),
                             gam = predict(m_gam, env_df[,my_preds], type='response'),
                             rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                             brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    
    # Save ensemble raster in prepared rasterstack
    r_curr_preds_year <- c(r_curr_preds_year, r_curr_preds_ens)
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_noclim_landuse <- c(r_curr_preds_noclim_landuse, r_curr_preds_year)
  
  
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_noclim_landuse, filename = "output_data/results/I_ricinus_preds_noclim_landuse_1970_2019.tif", overwrite=T)




#-------------------------------------------------------------------------------

# 3. Past monthly predictions from 1970 to 2019 --------------------------------
# only under climate change (keeping land use constant at their mean across Europe)

# Prepare a vector containing the years for monthly predictions
years <- c(1970:2019)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Prepare path to environmental data
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Load a raster as example template 
example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster

# Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
r_curr_preds_clim <- terra::rast(example_data, nlyrs = 600)


for (y in years) { # Start of the loop over the prediction years
  
  print(y)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land use data raster
  
  # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
  r_curr_preds_year <- terra::rast(example_data, nlyrs = 12)
  
  for (m in month) { # Start of the loop over all months of a year
    
    print(m)
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data)
    
    # Prepare a data frame with environmental data
    env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
    
    # Keep the values of land use variables constant at their mean
    for (pred in c("urban", "primary_forest", "cropland", "secondary_forest", "primary_openland", "pasture", "rangeland")) {
      env_df[[pred]] <- mean(env_df[[pred]], na.rm = TRUE)
    }
    
    # Make predictions of all models
    print("start of model predictions")
    curr_preds <- data.frame(env_df[,1:2], 
                             glm = predict(m_glm, env_df, type='response'),
                             gam = predict(m_gam, env_df[,my_preds], type='response'),
                             rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                             brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
    
    # Make ensemble predictions
    curr_preds$mean_prob = rowMeans(curr_preds[,-c(1:2)])
    
    # Make Spatrasters from predictions
    r_curr_preds <- terra::rast(curr_preds, crs = crs(env_data))
    
    # Extract the ensemble raster
    r_curr_preds_ens <- r_curr_preds[["mean_prob"]]
    
    # Save ensemble raster in prepared rasterstack
    r_curr_preds_year <- c(r_curr_preds_year, r_curr_preds_ens)
    
    
  } # End of loop over all months
  
  
  # Make sure that raster names are correct
  names(r_curr_preds_year) <- sprintf("%02d/%s", 1:12, y)
  
  # Stack the all raster for each year
  r_curr_preds_clim <- c(r_curr_preds_clim, r_curr_preds_year)
  
  
  
} # End of the loop over all considered years


# Save the raster outputs
terra::writeRaster(r_curr_preds_clim, filename = "output_data/results/I_ricinus_preds_clim_1970_2019.tif", overwrite=T)




#-------------------------------------------------------------------------------

# 4. Future monthly predictions from 2030 to 2070 ------------------------------
# under climate change (3 different scenarios), keeping land use steady at 2019


# Prepare a vector containing the years for monthly predictions
years <- c(2030:2070)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Create path to directory containing environmental data
datapath_env_fut <- file.path("input_data/environmental_data/ISIMIP3b/")
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Prepare a vector containing the three different climate forcing scenarios
scenario <- c("ssp126", "ssp370", "ssp585")

for (s in scenario) { # Start of the loop over the three different climate forcing scenarios
  
  print(s)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
  
  # Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
  r_fut_preds_clim_nolanduse <- terra::rast(example_data, nlyrs = 492)
  
  for (y in years) { # Start of the loop over the prediction years
    
    print(y)
    
    # Load a raster as example template 
    example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
    
    # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
    r_fut_preds_year <- terra::rast(example_data, nlyrs = 12)
    
    for (m in month) { # Start of the loop over all months of a year
      
      print(m)
      
      # Load the environmental data for the specific year and month
      Climate_data <- terra::rast(paste0(datapath_env_fut, "/Climate/",s,"/processed_data/Climate_future_data_",m,"_",y,"_",s,".tif")) # Climate data
      LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land cover data
      
      # Stack the environmental data
      env_data <- c(Climate_data, LandUse_data)
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      fut_preds <- data.frame(env_df[,1:2], 
                              glm = predict(m_glm, env_df, type='response'),
                              gam = predict(m_gam, env_df[,my_preds], type='response'),
                              rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                              brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
      
      # Make ensemble predictions
      fut_preds$mean_prob = rowMeans(fut_preds[,-c(1:2)])
      
      # Make Spatrasters from predictions
      r_fut_preds <- terra::rast(fut_preds, crs = crs(env_data))
      
      # Extract the ensemble raster
      r_fut_preds_ens <- r_fut_preds[["mean_prob"]]
      
      # Save ensemble raster in prepared rasterstack
      r_fut_preds_year <- c(r_fut_preds_year, r_fut_preds_ens)
      
      
    } # End of loop over all months
    
    
    # Make sure that raster names are correct
    names(r_fut_preds_year) <- sprintf("%02d/%s", 1:12, y)
    
    # Stack the all raster for each year
    r_fut_preds_clim_nolanduse <- c(r_fut_preds_clim_nolanduse, r_fut_preds_year)
    
    
    
  } # End of the loop over all considered years
  
  
  # Save the raster outputs
  terra::writeRaster(r_fut_preds_clim_nolanduse, filename = paste0("output_data/results/I_ricinus_preds_clim_nolanduse_2030_2070_",s,".tif"), overwrite=T)
  
  
} # Close the loop over the three climate forcing scenarios




#-------------------------------------------------------------------------------

# 5. Future monthly predictions from 2030 to 2070 ------------------------------
# under climate change and land use change (3 different scenarios)

# Prepare a vector containing the years for monthly predictions
years <- c(2030:2070)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Create path to directory containing environmental data
datapath_env_fut <- file.path("input_data/environmental_data/ISIMIP3b/")
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Prepare a vector containing the three different forcing scenarios
scenario <- c("ssp126", "ssp370", "ssp585")

for (s in scenario) { # Start of the loop over the three different forcing scenarios
  
  print(s)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
  
  # Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
  r_fut_preds_clim_landuse <- terra::rast(example_data, nlyrs = 492)
  
  for (y in years) { # Start of the loop over the prediction years
    
    print(y)
    
    # Load a raster as example template 
    example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
    
    # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
    r_fut_preds_year <- terra::rast(example_data, nlyrs = 12)
    
    for (m in month) { # Start of the loop over all months of a year
      
      print(m)
      
      # Load the environmental data for the specific year and month
      Climate_data <- terra::rast(paste0(datapath_env_fut, "/Climate/",s,"/processed_data/Climate_future_data_",m,"_",y,"_",s,".tif")) # Climate data
      LandUse_data <- terra::rast(paste0(datapath_env_fut, "/LandUse/",s,"/processed_data/LandUse_future_data_",y,"_",s,".tif")) # Land cover data
      
      # Stack the environmental data
      env_data <- c(Climate_data, LandUse_data)
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Make predictions of all models
      print("start of model predictions")
      fut_preds <- data.frame(env_df[,1:2], 
                              glm = predict(m_glm, env_df, type='response'),
                              gam = predict(m_gam, env_df[,my_preds], type='response'),
                              rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                              brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
      
      # Make ensemble predictions
      fut_preds$mean_prob = rowMeans(fut_preds[,-c(1:2)])
      
      # Make Spatrasters from predictions
      r_fut_preds <- terra::rast(fut_preds, crs = crs(env_data))
      
      # Extract the ensemble raster
      r_fut_preds_ens <- r_fut_preds[["mean_prob"]]
      
      # Save ensemble raster in prepared rasterstack
      r_fut_preds_year <- c(r_fut_preds_year, r_fut_preds_ens)
      
      
    } # End of loop over all months
    
    
    # Make sure that raster names are correct
    names(r_fut_preds_year) <- sprintf("%02d/%s", 1:12, y)
    
    # Stack the all raster for each year
    r_fut_preds_clim_landuse <- c(r_fut_preds_clim_landuse, r_fut_preds_year)
    
    
    
  } # End of the loop over all considered years
  
  
  # Save the raster outputs
  terra::writeRaster(r_fut_preds_clim_landuse, filename = paste0("output_data/results/I_ricinus_preds_clim_landuse_2030_2070_",s,".tif"), overwrite=T)
  
  
} # Close the loop over the three forcing scenarios



#-------------------------------------------------------------------------------

# 6. Future monthly predictions from 2030 to 2070 ------------------------------
# only under climate change (keeping land use constant at their mean across Europe)

# Prepare a vector containing the years for monthly predictions
years <- c(2030:2070)

# Prepare vector containing the months of prediction
month <- str_pad(1:12, width = 2, pad = "0")

# Create path to directory containing environmental data
datapath_env_fut <- file.path("input_data/environmental_data/ISIMIP3b/")
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Prepare a vector containing the three different forcing scenarios
scenario <- c("ssp126", "ssp370", "ssp585")

for (s in scenario) { # Start of the loop over the three different forcing scenarios
  
  print(s)
  
  # Load a raster as example template 
  example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
  
  # Prepare a Spatraster to store all prediction rasters of all prediction years, using the env_data as template
  r_fut_preds_clim <- terra::rast(example_data, nlyrs = 492)
  
  for (y in years) { # Start of the loop over the prediction years
    
    print(y)
    
    # Load a raster as example template 
    example_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_2019.tif")) # Land use data raster
    
    # Prepare a Spatraster to store the prediction rasters of one prediction year, using the env_data as template
    r_fut_preds_year <- terra::rast(example_data, nlyrs = 12)
    
    for (m in month) { # Start of the loop over all months of a year
      
      print(m)
      
      # Load the environmental data for the specific year and month
      Climate_data <- terra::rast(paste0(datapath_env_fut, "/Climate/",s,"/processed_data/Climate_future_data_",m,"_",y,"_",s,".tif")) # Climate data
      LandUse_data <- terra::rast(paste0(datapath_env_fut, "/LandUse/",s,"/processed_data/LandUse_future_data_",y,"_",s,".tif")) # Land cover data
      
      # Stack the environmental data
      env_data <- c(Climate_data, LandUse_data)
      
      # Prepare a data frame with environmental data
      env_df <- data.frame(crds(env_data[[my_preds]]),as.points(env_data[[my_preds]]))
      
      # Keep the values of land use variables constant at their mean
      for (pred in c("urban", "primary_forest", "cropland", "secondary_forest", "primary_openland", "pasture", "rangeland")) {
        env_df[[pred]] <- mean(env_df[[pred]], na.rm = TRUE)
      }
      
      # Make predictions of all models
      print("start of model predictions")
      fut_preds <- data.frame(env_df[,1:2], 
                              glm = predict(m_glm, env_df, type='response'),
                              gam = predict(m_gam, env_df[,my_preds], type='response'),
                              rf = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict(m_rf[[i]], env_df, type='response')})),
                              brt = rowMeans(sapply(1:ratio_presence_background, FUN=function(i){print(i); predict.gbm(m_brt[[i]], env_df, n.trees=m_brt[[i]]$gbm.call$best.trees, type="response")})))
      
      # Make ensemble predictions
      fut_preds$mean_prob = rowMeans(fut_preds[,-c(1:2)])
      
      # Make Spatrasters from predictions
      r_fut_preds <- terra::rast(fut_preds, crs = crs(env_data))
      
      # Extract the ensemble raster
      r_fut_preds_ens <- r_fut_preds[["mean_prob"]]
      
      # Save ensemble raster in prepared rasterstack
      r_fut_preds_year <- c(r_fut_preds_year, r_fut_preds_ens)
      
      
    } # End of loop over all months
    
    
    # Make sure that raster names are correct
    names(r_fut_preds_year) <- sprintf("%02d/%s", 1:12, y)
    
    # Stack the all raster for each year
    r_fut_preds_clim <- c(r_fut_preds_clim, r_fut_preds_year)
    
    
    
  } # End of the loop over all considered years
  
  
  # Save the raster outputs
  terra::writeRaster(r_fut_preds_clim, filename = paste0("output_data/results/I_ricinus_preds_clim_2030_2070_",s,".tif"), overwrite=T)
  
  
} # Close the loop over the three forcing scenarios


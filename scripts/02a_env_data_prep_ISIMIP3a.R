# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#             02a. Historical environmental data preparation             #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We process historical ISIMIP3a (Inter-Sectoral Impact Model Intercomparison 
# Project Phase 3a) climate and land-use datasets, including both factual and 
# counterfactual scenarios. The climate data are provided as daily outputs at 
# a 0.5° spatial resolution. To match our target temporal resolution, we aggregate 
# daily temperature and humidity values into monthly means and daily precipitation 
# into monthly totals. In addition, the annual land-use data at 0.5° resolution 
# are aggregated into eight distinct land-use categories for each year.

# Load needed packages
library(terra) # terra_1.7-55
library(tidyverse) # tidyverse_2.0.0



#-------------------------------------------------------------------------------

# 1. Create mask of Europe (0.5° - WGS84) --------------------------------------

# Load needed objects
# Shapefile of Europe (downloaded from the ArcGIS Hub: https://hub.arcgis.com/datasets/bdcb40c0b6124f6d99f10b9b23647712/explore)
europe_shp <- terra::vect("input_data/spatial_data/Europe_NUTS_0_Demographics_and_Boundaries.shp")

# Define the extent of Europe with an approximate bounding box
europe_extent <- ext(-31, 40, 34, 72)

# Create an empty raster with a resolution of 0.5°
europe_raster <- terra::rast(europe_extent, resolution = c(0.5, 0.5))

# Rasterise the Europe shapefile to create a mask
europe_mask <- terra::rasterize(europe_shp, europe_raster, background = NA)

# Save the resulting raster
terra::writeRaster(europe_mask, filename = "input_data/spatial_data/europe_mask.tif", overwrite = TRUE)





#-------------------------------------------------------------------------------

# 2. Prepare observed climate data ---------------------------------------------

# CHELSA GSWP3-W5E5 (recommended by Dirk Karger)
# https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/ -> configure download
# Bounding box: South: 34 North: 75 West: -31 East: 60
# Extracted files from zip document
# 1970-2019


# Load in all data of the climate variables (pr, tas, tasmax, tasmin, hurs) 
# and calculate monthly averages (temperature, relative humidity) 
# or monthly sums (precipitation)
# for each year

# Prepare path to data folder
datapath_Climate_data <- file.path("input_data/environmental_data/ISIMIP3a/Climate")

# Create a vector containing the months of a year
months <- str_pad(1:12, width = 2, pad = "0")

# Create a vector containing the years of interest
years <- 1970:2019

# Create a vector containing the five climatic target variables
clim_variables <- c("pr", "tas", "tasmax", "tasmin", "hurs")

# Loop through the different year-month combinations for each variable,
# process climate data and write a raster for each year-month combination 
# containing all five climate variables

for (y in years) { # Start of the loop over the different years
  
  print(y)
  
  for (m in months) { # Start of the loop over the different months
    
    print(m)
    
    # Check if climate file was already processed for respective month and year
    file_exists <- file.exists(paste0(datapath_Climate_data, "/processed_data/Climate_data_",m,"_",y,".tif"))
    
    if (file_exists == FALSE) { # If file does not exist, start processing
      
      # Write a vector with the year-month combination
      year_month_select <- paste(y, m, sep = "-")
      
      for (c in clim_variables) { # Start of the loop over the five different climate variables
        
        print(c)
        
        # Load climate files that contain data for the years from 1970 to 2019
        # for the respective variable
        Climate_data_rasters_1 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1961_1970.nc"))
        Climate_data_rasters_2 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1971_1980.nc"))
        Climate_data_rasters_3 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1981_1990.nc"))
        Climate_data_rasters_4 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1991_2000.nc"))
        Climate_data_rasters_5 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2001_2010.nc"))
        Climate_data_rasters_6 <- terra::rast(paste0(datapath_Climate_data, "/raw_data/gswp3-w5e5_obsclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2011_2019.nc"))
        
        # Stack all rasters
        Climate_data_rasters <- c(Climate_data_rasters_1, Climate_data_rasters_2, Climate_data_rasters_3,
                                  Climate_data_rasters_4, Climate_data_rasters_5, Climate_data_rasters_6)
        
        # Extract the time information from the rasters, starting from the year 1970
        Climate_data_rasters_dates <- as.Date(time(Climate_data_rasters), origin = "1970-01-01")
        
        # Format these dates to only contain month and year
        Climate_data_rasters_dates <- format(Climate_data_rasters_dates, "%Y-%m")
        
        # Extract daily rasters that correspond to that year-month combination
        year_month_rasters <- which(Climate_data_rasters_dates == year_month_select)
        year_month_rasters <- Climate_data_rasters[[year_month_rasters]]

        
        # Calculate mean values for corresponding month for temperature variables 
        # and relative humidity variable, monthly sums for precipitation data 
        # and adapt the unit of variable
        if (c == "pr") {
        
        # Precipitation is given in the unit kg m-2 s-1, calculate to kg m-2 month-1
        days_month <- nlyr(year_month_rasters) # Extract the days of month by number of raster layers
        seconds_day <- 24 * 3600 # Extract the number of seconds for the day
        year_month_rasters_proc <- sum(year_month_rasters * seconds_day) # Calculate new unit: Multiply each daily flux by seconds in a day, then sum over all days of a month
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc) # Assign raster an adjusted name
        
        
        } else if (c %in% c("tas", "tasmax", "tasmin")) { year_month_rasters_proc <- mean(year_month_rasters) # Average Temperature values (given in K) over each month
        
        # Temperature is given in the unit K, calculate to °C
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc - 273.15) # Subtract 273.15 to get °C values and assign raster an adjusted name
        
        
        } else if (c == "hurs") { year_month_rasters_proc <- mean(year_month_rasters) # Average relative humidity values (given in %) over each month
        
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc) # Assign raster an adjusted name
          
        
        } # End of if-condition
        
        
      } # Close the loop over all five climate variables
      
      # Stack the climatic variables of the same year-month combination
      Climate_rasters_processed <- c(year_month_rasters_processed_pr, year_month_rasters_processed_tas,
                                     year_month_rasters_processed_tasmax, year_month_rasters_processed_tasmin,
                                     year_month_rasters_processed_hurs)
      
      # Crop the extent to match the bounding box of the mask and mask the values 
      # outside of the terrestrial continent of Europe (e.g. ocean area)
      Climate_rasters_processed <- terra::crop(Climate_rasters_processed, europe_mask)
      Climate_rasters_processed <- terra::mask(Climate_rasters_processed, europe_mask)
      
      # Add names to the raster layers
      names(Climate_rasters_processed) <- c("pr", "tas", "tasmax", "tasmin", "hurs")
      
      # Save the processed raster as tif file
      terra::writeRaster(Climate_rasters_processed, filename = paste0(datapath_Climate_data, "/processed_data/Climate_data_",m,"_",y,".tif"), overwrite = TRUE)
      
      
      
    } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next month
    } # Close if condition
    
    
    
  } # Close the loop over all months of a year
  
} # Close the loop over all years




#-------------------------------------------------------------------------------

# 3. Prepare detrended climate data --------------------------------------------

# Scenario of no climate change (counterclim)
# CHELSA GSWP3-W5E5 (recommended by Dirk Karger)
# https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/atmosphere/gswp3-w5e5/ -> configure download
# Bounding box: South: 34 North: 75 West: -31 East: 60
# Extracted files from zip document
# 1970-2019


# Load in all data of the climate variables (pr, tas, tasmax, tasmin, hurs) 
# and calculate monthly averages (temperature, relative humidity) 
# or monthly sums (precipitation)
# for each year

# Prepare path to data folder
datapath_CounterClim_data <- file.path("input_data/environmental_data/ISIMIP3a/CounterClim")

# Create a vector containing the months of a year
months <- str_pad(1:12, width = 2, pad = "0")

# Create a vector containing the years of interest
years <- 1970:2019

# Create a vector containing the four climatic target variables
clim_variables <- c("pr", "tas", "tasmax", "tasmin", "hurs")

# Loop through the different year-month combinations for each variable,
# process climate data and write a raster for each year-month combination 
# containing all five climate variables

for (y in years) { # Start of the loop over the different years
  
  print(y)
  
  for (m in months) { # Start of the loop over the different months
    
    print(m)
    
    # Check if climate file was already processed for respective month and year
    file_exists <- file.exists(paste0(datapath_CounterClim_data, "/processed_data/CounterClim_data_",m,"_",y,".tif"))
    
    if (file_exists == FALSE) { # If file does not exist, start processing
      
      # Write a vector with the year-month combination
      year_month_select <- paste(y, m, sep = "-")
      
      for (c in clim_variables) { # Start of the loop over the five different climate variables
        
        print(c)
        
        # Load climate files that contain data for the years from 1970 to 2019
        # for the respective variable
        CounterClim_data_rasters_1 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1961_1970.nc"))
        CounterClim_data_rasters_2 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1971_1980.nc"))
        CounterClim_data_rasters_3 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1981_1990.nc"))
        CounterClim_data_rasters_4 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_1991_2000.nc"))
        CounterClim_data_rasters_5 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2001_2010.nc"))
        CounterClim_data_rasters_6 <- terra::rast(paste0(datapath_CounterClim_data, "/raw_data/gswp3-w5e5_counterclim_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2011_2019.nc"))
        
        # Stack all rasters
        CounterClim_data_rasters <- c(CounterClim_data_rasters_1, CounterClim_data_rasters_2, CounterClim_data_rasters_3,
                                      CounterClim_data_rasters_4, CounterClim_data_rasters_5, CounterClim_data_rasters_6)
        
        # Extract the time information from the rasters, starting from the year 1970
        CounterClim_data_rasters_dates <- as.Date(time(CounterClim_data_rasters), origin = "1970-01-01")
        
        # Format these dates to only contain month and year
        CounterClim_data_rasters_dates <- format(CounterClim_data_rasters_dates, "%Y-%m")
        
        # Extract daily rasters that correspond to that year-month combination
        year_month_rasters <- which(CounterClim_data_rasters_dates == year_month_select)
        year_month_rasters <- CounterClim_data_rasters[[year_month_rasters]]
        

        
        # Calculate mean values for corresponding month for temperature variables 
        # and relative humidity variable, monthly sums for precipitation data 
        # and adapt the unit of variable
        if (c == "pr") {
        
        # Precipitation is given in the unit kg m-2 s-1, calculate to kg m-2 month-1
        days_month <- nlyr(year_month_rasters) # Extract the days of month by number of raster layers
        seconds_day <- 24 * 3600 # Extract the number of seconds for the day
        year_month_rasters_proc <- sum(year_month_rasters * seconds_day) # Calculate new unit: Multiply each daily flux by seconds in a day, then sum over all days of a month
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc) # Assign raster an adjusted name
        
        
        } else if (c %in% c("tas", "tasmax", "tasmin")) { year_month_rasters_proc <- mean(year_month_rasters) # Average Temperature values (given in K) over each month
        
        # Temperature is given in the unit K, calculate to °C
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc - 273.15) # Subtract 273.15 to get °C values and assign raster an adjusted name
        
        
        } else if (c == "hurs") { year_month_rasters_proc <- mean(year_month_rasters) # Average relative humidity values (given in %) over each month
        
        assign(paste0("year_month_rasters_processed_", c), year_month_rasters_proc) # Assign raster an adjusted name# Save raster with adjusted name
        
        
        } # End of if-condition
        
        
      } # Close the loop over all five climate variables
      
      # Stack the climatic variables of the same year-month combination
      CounterClim_rasters_processed <- c(year_month_rasters_processed_pr, year_month_rasters_processed_tas,
                                         year_month_rasters_processed_tasmax, year_month_rasters_processed_tasmin,
                                         year_month_rasters_processed_hurs)
      
      # Crop the extent to match the bounding box of the mask and mask the values 
      # outside of the terrestrial continent of Europe (e.g. ocean area)
      CounterClim_rasters_processed <- terra::crop(CounterClim_rasters_processed, europe_mask)
      CounterClim_rasters_processed <- terra::mask(CounterClim_rasters_processed, europe_mask)
      
      # Add names to the raster layers
      names(CounterClim_rasters_processed) <- c("pr", "tas", "tasmax", "tasmin", "hurs")
      
      # Save the processed raster as tif file
      terra::writeRaster(CounterClim_rasters_processed, filename = paste0(datapath_CounterClim_data, "/processed_data/CounterClim_data_",m,"_",y,".tif"), overwrite = TRUE)
      
      
      
    } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next month
    } # Close if condition
    
    
    
  } # Close the loop over all months of a year
  
} # Close the loop over all years




#-------------------------------------------------------------------------------

# 4. Prepare historical land use data ------------------------------------------

# Four different datasets downloaded (Downscaling to 5 crops, Forests and natural
# vegetation, Managed pastures, and Urban areas) from:
# https://data.isimip.org/search/tree/ISIMIP3a/InputData/socioeconomic/landuse/histsoc
  

# Prepare path to data folder
datapath_LandUse_data <- file.path("input_data/environmental_data/ISIMIP3a/LandUse")

# Create a vector containing the years of interest
years <- 1970:2019

# Load in yearly land use data
LandUse_data_rasters_crops <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-5crops_histsoc_annual_1901_2021.nc"))
LandUse_data_rasters_fornatveg <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-forests-and-natural-vegetation_histsoc_annual_1901_2021.nc"))
LandUse_data_rasters_pastures <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-pastures_histsoc_annual_1901_2021.nc"))
LandUse_data_rasters_urban <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-urbanareas_histsoc_annual_1901_2021.nc"))

# Extract the rasters for each variable of interest
print("forested primary land")
LandUse_data_rasters_primf <- LandUse_data_rasters_fornatveg[[grep("primary_forests", names(LandUse_data_rasters_fornatveg))]]

print("non-forested primary land")
LandUse_data_rasters_primn <- LandUse_data_rasters_fornatveg[[grep("primary_nonforests", names(LandUse_data_rasters_fornatveg))]]

print("potentially forested secondary land")
LandUse_data_rasters_secdf <- LandUse_data_rasters_fornatveg[[grep("secondary_forests", names(LandUse_data_rasters_fornatveg))]]

print("potentially non-forested secondary land")
LandUse_data_rasters_secdn <- LandUse_data_rasters_fornatveg[[grep("secondary_nonforests", names(LandUse_data_rasters_fornatveg))]]

print("managed pasture")
LandUse_data_rasters_pastr <- LandUse_data_rasters_pastures[[grep("managed_pastures", names(LandUse_data_rasters_pastures))]]

print("rangeland")
LandUse_data_rasters_range <- LandUse_data_rasters_pastures[[grep("rangeland", names(LandUse_data_rasters_pastures))]]

print("C3 annual crops")
LandUse_data_rasters_crop1 <- LandUse_data_rasters_crops[[grep("c3ann_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C4 annual crops")
LandUse_data_rasters_crop2 <- LandUse_data_rasters_crops[[grep("c4ann_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C3 perennial crops")
LandUse_data_rasters_crop3 <- LandUse_data_rasters_crops[[grep("c3per_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C4 perennial crops")
LandUse_data_rasters_crop4 <- LandUse_data_rasters_crops[[grep("c4per_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C3 nitrogen-fixing crops")
LandUse_data_rasters_crop5 <- LandUse_data_rasters_crops[[grep("c3nfx_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]


for (y in years) { # Start of the loop over the years of interest
  
  print(y)
  
  # Check if land use file was already processed for respective year
  file_exists <- file.exists(paste0(datapath_LandUse_data, "/processed_data/LandUse_data_",y,".tif"))
  
  if (file_exists == FALSE) { # If file does not exist, start processing
    
    # Extract the raster of the respective year for each land use variable,
    # crop it to the European extent and mask it to the terrestrial European
    # area
    # (in the case of cropland grids, rainfed and irrigated croplands are added 
    # up for each cropland type before they are summarised into a general 
    # cropland layer)
    print("forested primary land")
    LandUse_data_raster_primf_year <- LandUse_data_rasters_primf[[grep(y, time(LandUse_data_rasters_primf))]]
    LandUse_data_raster_primf_year <- terra::crop(LandUse_data_raster_primf_year, europe_mask)
    LandUse_data_raster_primf_year <- terra::mask(LandUse_data_raster_primf_year, europe_mask)
    
    print("non-forested primary land")
    LandUse_data_raster_primn_year <- LandUse_data_rasters_primn[[grep(y, time(LandUse_data_rasters_primn))]]
    LandUse_data_raster_primn_year <- terra::crop(LandUse_data_raster_primn_year, europe_mask)
    LandUse_data_raster_primn_year <- terra::mask(LandUse_data_raster_primn_year, europe_mask)
    
    print("potentially forested secondary land")
    LandUse_data_raster_secdf_year <- LandUse_data_rasters_secdf[[grep(y, time(LandUse_data_rasters_secdf))]]
    LandUse_data_raster_secdf_year <- terra::crop(LandUse_data_raster_secdf_year, europe_mask)
    LandUse_data_raster_secdf_year <- terra::mask(LandUse_data_raster_secdf_year, europe_mask)
    
    print("potentially non-forested secondary land")
    LandUse_data_raster_secdn_year <- LandUse_data_rasters_secdn[[grep(y, time(LandUse_data_rasters_secdn))]]
    LandUse_data_raster_secdn_year <- terra::crop(LandUse_data_raster_secdn_year, europe_mask)
    LandUse_data_raster_secdn_year <- terra::mask(LandUse_data_raster_secdn_year, europe_mask)
    
    print("managed pasture")
    LandUse_data_raster_pastr_year <- LandUse_data_rasters_pastr[[grep(y, time(LandUse_data_rasters_pastr))]]
    LandUse_data_raster_pastr_year <- terra::crop(LandUse_data_raster_pastr_year, europe_mask)
    LandUse_data_raster_pastr_year <- terra::mask(LandUse_data_raster_pastr_year, europe_mask)
    
    print("rangeland")
    LandUse_data_raster_range_year <- LandUse_data_rasters_range[[grep(y, time(LandUse_data_rasters_range))]]
    LandUse_data_raster_range_year <- terra::crop(LandUse_data_raster_range_year, europe_mask)
    LandUse_data_raster_range_year <- terra::mask(LandUse_data_raster_range_year, europe_mask)
    
    print("C3 annual crops")
    LandUse_data_raster_crop1_year <- LandUse_data_rasters_crop1[[grep(y, time(LandUse_data_rasters_crop1))]]
    LandUse_data_raster_crop1_year <- terra::crop(LandUse_data_raster_crop1_year, europe_mask)
    LandUse_data_raster_crop1_year <- terra::mask(LandUse_data_raster_crop1_year, europe_mask)
    LandUse_data_raster_crop1_year <- LandUse_data_raster_crop1_year[[1]] + LandUse_data_raster_crop1_year[[2]]
    
    print("C4 annual crops")
    LandUse_data_raster_crop2_year <- LandUse_data_rasters_crop2[[grep(y, time(LandUse_data_rasters_crop2))]]
    LandUse_data_raster_crop2_year <- terra::crop(LandUse_data_raster_crop2_year, europe_mask)
    LandUse_data_raster_crop2_year <- terra::mask(LandUse_data_raster_crop2_year, europe_mask)
    LandUse_data_raster_crop2_year <- LandUse_data_raster_crop2_year[[1]] + LandUse_data_raster_crop2_year[[2]]
    
    print("C3 perennial crops")
    LandUse_data_raster_crop3_year <- LandUse_data_rasters_crop3[[grep(y, time(LandUse_data_rasters_crop3))]]
    LandUse_data_raster_crop3_year <- terra::crop(LandUse_data_raster_crop3_year, europe_mask)
    LandUse_data_raster_crop3_year <- terra::mask(LandUse_data_raster_crop3_year, europe_mask)
    LandUse_data_raster_crop3_year <- LandUse_data_raster_crop3_year[[1]] + LandUse_data_raster_crop3_year[[2]]
    
    print("C4 perennial crops")
    LandUse_data_raster_crop4_year <- LandUse_data_rasters_crop4[[grep(y, time(LandUse_data_rasters_crop4))]]
    LandUse_data_raster_crop4_year <- terra::crop(LandUse_data_raster_crop4_year, europe_mask)
    LandUse_data_raster_crop4_year <- terra::mask(LandUse_data_raster_crop4_year, europe_mask)
    LandUse_data_raster_crop4_year <- LandUse_data_raster_crop4_year[[1]] + LandUse_data_raster_crop4_year[[2]]
    
    print("C3 nitrogen-fixing crops")
    LandUse_data_raster_crop5_year <- LandUse_data_rasters_crop5[[grep(y, time(LandUse_data_rasters_crop5))]]
    LandUse_data_raster_crop5_year <- terra::crop(LandUse_data_raster_crop5_year, europe_mask)
    LandUse_data_raster_crop5_year <- terra::mask(LandUse_data_raster_crop5_year, europe_mask)
    LandUse_data_raster_crop5_year <- LandUse_data_raster_crop5_year[[1]] + LandUse_data_raster_crop5_year[[2]]
    
    print("urban land")
    LandUse_data_raster_urban_year <- LandUse_data_rasters_urban[[grep(y, time(LandUse_data_rasters_urban))]]
    LandUse_data_raster_urban_year <- terra::crop(LandUse_data_raster_urban_year, europe_mask)
    LandUse_data_raster_urban_year <- terra::mask(LandUse_data_raster_urban_year, europe_mask)
    
    
    # Summarise variable types (cropland)
    LandUse_data_cropland <- sum(LandUse_data_raster_crop1_year, LandUse_data_raster_crop2_year, 
                                 LandUse_data_raster_crop3_year, LandUse_data_raster_crop4_year,
                                 LandUse_data_raster_crop5_year)
    
    # Stack the raster files for the respective year
    LandUse_rasters_processed <- c(LandUse_data_raster_primf_year, LandUse_data_raster_primn_year,
                                   LandUse_data_raster_secdf_year, LandUse_data_raster_secdn_year,
                                   LandUse_data_raster_pastr_year, LandUse_data_raster_range_year,
                                   LandUse_data_cropland, LandUse_data_raster_urban_year)
    
    # Add names to the raster layers
    names(LandUse_rasters_processed) <- c("primary_forest", "primary_openland", "secondary_forest", "secondary_openland",
                                          "pasture", "rangeland", "cropland", "urban")
    
    # Save the processed raster files
    terra::writeRaster(LandUse_rasters_processed, filename = paste0(datapath_LandUse_data, "/processed_data/LandUse_data_",y,".tif"), overwrite = TRUE)
    
    
  } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next year
  } # Close if condition
  
} # End of the loop over all years of interest




#-------------------------------------------------------------------------------

# 5. Prepare counterfactual land use data --------------------------------------

# Four different datasets downloaded (Downscaling to 5 crops, Forests and natural
# vegetation, Managed pastures, and Urban areas) from:
# https://data.isimip.org/search/tree/ISIMIP3a/InputData/socioeconomic/landuse/1901soc

# Prepare path to data folder
datapath_LandUse_data <- file.path("input_data/environmental_data/ISIMIP3a/CounterLandUse")

# Create a vector containing the years of interest
years <- 1970:2019

# Load in yearly land use data
LandUse_data_rasters_crops <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-5crops_1901soc_annual_1901_2021.nc"))
LandUse_data_rasters_fornatveg <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-forests-and-natural-vegetation_1901soc_annual_1901_2021.nc"))
LandUse_data_rasters_pastures <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-pastures_1901soc_annual_1901_2021.nc"))
LandUse_data_rasters_urban <- terra::rast(paste0(datapath_LandUse_data, "/raw_data/landuse-urbanareas_1901soc_annual_1901_2021.nc"))

# Extract the rasters for each variable of interest
print("forested primary land")
LandUse_data_rasters_primf <- LandUse_data_rasters_fornatveg[[grep("primary_forests", names(LandUse_data_rasters_fornatveg))]]

print("non-forested primary land")
LandUse_data_rasters_primn <- LandUse_data_rasters_fornatveg[[grep("primary_nonforests", names(LandUse_data_rasters_fornatveg))]]

print("potentially forested secondary land")
LandUse_data_rasters_secdf <- LandUse_data_rasters_fornatveg[[grep("secondary_forests", names(LandUse_data_rasters_fornatveg))]]

print("potentially non-forested secondary land")
LandUse_data_rasters_secdn <- LandUse_data_rasters_fornatveg[[grep("secondary_nonforests", names(LandUse_data_rasters_fornatveg))]]

print("managed pasture")
LandUse_data_rasters_pastr <- LandUse_data_rasters_pastures[[grep("managed_pastures", names(LandUse_data_rasters_pastures))]]

print("rangeland")
LandUse_data_rasters_range <- LandUse_data_rasters_pastures[[grep("rangeland", names(LandUse_data_rasters_pastures))]]

print("C3 annual crops")
LandUse_data_rasters_crop1 <- LandUse_data_rasters_crops[[grep("c3ann_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C4 annual crops")
LandUse_data_rasters_crop2 <- LandUse_data_rasters_crops[[grep("c4ann_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C3 perennial crops")
LandUse_data_rasters_crop3 <- LandUse_data_rasters_crops[[grep("c3per_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C4 perennial crops")
LandUse_data_rasters_crop4 <- LandUse_data_rasters_crops[[grep("c4per_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]

print("C3 nitrogen-fixing crops")
LandUse_data_rasters_crop5 <- LandUse_data_rasters_crops[[grep("c3nfx_(irrigated|rainfed)", names(LandUse_data_rasters_crops))]]


for (y in years) { # Start of the loop over the years of interest
  
  print(y)
  
  # Check if land use file was already processed for respective year
  file_exists <- file.exists(paste0(datapath_LandUse_data, "/processed_data/LandUse_data_",y,".tif"))
  
  if (file_exists == FALSE) { # If file does not exist, start processing
    
    # Extract the raster of the respective year for each land use variable,
    # crop it to the European extent and mask it to the terrestrial European
    # area
    # (in the case of cropland grids, rainfed and irrigated croplands are added 
    # up for each cropland type before they are summarised into a general 
    # cropland layer)
    print("forested primary land")
    LandUse_data_raster_primf_year <- LandUse_data_rasters_primf[[grep(y, time(LandUse_data_rasters_primf))]]
    LandUse_data_raster_primf_year <- terra::crop(LandUse_data_raster_primf_year, europe_mask)
    LandUse_data_raster_primf_year <- terra::mask(LandUse_data_raster_primf_year, europe_mask)
    
    print("non-forested primary land")
    LandUse_data_raster_primn_year <- LandUse_data_rasters_primn[[grep(y, time(LandUse_data_rasters_primn))]]
    LandUse_data_raster_primn_year <- terra::crop(LandUse_data_raster_primn_year, europe_mask)
    LandUse_data_raster_primn_year <- terra::mask(LandUse_data_raster_primn_year, europe_mask)
    
    print("potentially forested secondary land")
    LandUse_data_raster_secdf_year <- LandUse_data_rasters_secdf[[grep(y, time(LandUse_data_rasters_secdf))]]
    LandUse_data_raster_secdf_year <- terra::crop(LandUse_data_raster_secdf_year, europe_mask)
    LandUse_data_raster_secdf_year <- terra::mask(LandUse_data_raster_secdf_year, europe_mask)
    
    print("potentially non-forested secondary land")
    LandUse_data_raster_secdn_year <- LandUse_data_rasters_secdn[[grep(y, time(LandUse_data_rasters_secdn))]]
    LandUse_data_raster_secdn_year <- terra::crop(LandUse_data_raster_secdn_year, europe_mask)
    LandUse_data_raster_secdn_year <- terra::mask(LandUse_data_raster_secdn_year, europe_mask)
    
    print("managed pasture")
    LandUse_data_raster_pastr_year <- LandUse_data_rasters_pastr[[grep(y, time(LandUse_data_rasters_pastr))]]
    LandUse_data_raster_pastr_year <- terra::crop(LandUse_data_raster_pastr_year, europe_mask)
    LandUse_data_raster_pastr_year <- terra::mask(LandUse_data_raster_pastr_year, europe_mask)
    
    print("rangeland")
    LandUse_data_raster_range_year <- LandUse_data_rasters_range[[grep(y, time(LandUse_data_rasters_range))]]
    LandUse_data_raster_range_year <- terra::crop(LandUse_data_raster_range_year, europe_mask)
    LandUse_data_raster_range_year <- terra::mask(LandUse_data_raster_range_year, europe_mask)
    
    print("C3 annual crops")
    LandUse_data_raster_crop1_year <- LandUse_data_rasters_crop1[[grep(y, time(LandUse_data_rasters_crop1))]]
    LandUse_data_raster_crop1_year <- terra::crop(LandUse_data_raster_crop1_year, europe_mask)
    LandUse_data_raster_crop1_year <- terra::mask(LandUse_data_raster_crop1_year, europe_mask)
    LandUse_data_raster_crop1_year <- LandUse_data_raster_crop1_year[[1]] + LandUse_data_raster_crop1_year[[2]]
    
    print("C4 annual crops")
    LandUse_data_raster_crop2_year <- LandUse_data_rasters_crop2[[grep(y, time(LandUse_data_rasters_crop2))]]
    LandUse_data_raster_crop2_year <- terra::crop(LandUse_data_raster_crop2_year, europe_mask)
    LandUse_data_raster_crop2_year <- terra::mask(LandUse_data_raster_crop2_year, europe_mask)
    LandUse_data_raster_crop2_year <- LandUse_data_raster_crop2_year[[1]] + LandUse_data_raster_crop2_year[[2]]
    
    print("C3 perennial crops")
    LandUse_data_raster_crop3_year <- LandUse_data_rasters_crop3[[grep(y, time(LandUse_data_rasters_crop3))]]
    LandUse_data_raster_crop3_year <- terra::crop(LandUse_data_raster_crop3_year, europe_mask)
    LandUse_data_raster_crop3_year <- terra::mask(LandUse_data_raster_crop3_year, europe_mask)
    LandUse_data_raster_crop3_year <- LandUse_data_raster_crop3_year[[1]] + LandUse_data_raster_crop3_year[[2]]
    
    print("C4 perennial crops")
    LandUse_data_raster_crop4_year <- LandUse_data_rasters_crop4[[grep(y, time(LandUse_data_rasters_crop4))]]
    LandUse_data_raster_crop4_year <- terra::crop(LandUse_data_raster_crop4_year, europe_mask)
    LandUse_data_raster_crop4_year <- terra::mask(LandUse_data_raster_crop4_year, europe_mask)
    LandUse_data_raster_crop4_year <- LandUse_data_raster_crop4_year[[1]] + LandUse_data_raster_crop4_year[[2]]
    
    print("C3 nitrogen-fixing crops")
    LandUse_data_raster_crop5_year <- LandUse_data_rasters_crop5[[grep(y, time(LandUse_data_rasters_crop5))]]
    LandUse_data_raster_crop5_year <- terra::crop(LandUse_data_raster_crop5_year, europe_mask)
    LandUse_data_raster_crop5_year <- terra::mask(LandUse_data_raster_crop5_year, europe_mask)
    LandUse_data_raster_crop5_year <- LandUse_data_raster_crop5_year[[1]] + LandUse_data_raster_crop5_year[[2]]
    
    print("urban land")
    LandUse_data_raster_urban_year <- LandUse_data_rasters_urban[[grep(y, time(LandUse_data_rasters_urban))]]
    LandUse_data_raster_urban_year <- terra::crop(LandUse_data_raster_urban_year, europe_mask)
    LandUse_data_raster_urban_year <- terra::mask(LandUse_data_raster_urban_year, europe_mask)
    
    # Summarize variable types (cropland)
    LandUse_data_cropland <- sum(LandUse_data_raster_crop1_year, LandUse_data_raster_crop2_year, 
                                 LandUse_data_raster_crop3_year, LandUse_data_raster_crop4_year,
                                 LandUse_data_raster_crop5_year)
    
    # Stack the raster files for the respective year
    LandUse_rasters_processed <- c(LandUse_data_raster_primf_year, LandUse_data_raster_primn_year,
                                   LandUse_data_raster_secdf_year, LandUse_data_raster_secdn_year,
                                   LandUse_data_raster_pastr_year, LandUse_data_raster_range_year,
                                   LandUse_data_cropland, LandUse_data_raster_urban_year) 
    
    # Add names to the raster layers
    names(LandUse_rasters_processed) <- c("primary_forest", "primary_openland", "secondary_forest", "secondary_openland",
                                          "pasture", "rangeland", "cropland", "urban")
    
    # Save the processed raster files
    terra::writeRaster(LandUse_rasters_processed, filename = paste0(datapath_LandUse_data, "/processed_data/CounterLandUse_data_",y,".tif"), overwrite = TRUE)
    
    
  } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next year
  } # Close if condition
  
} # End of the loop over all years of interest



# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#              02b. Future environmental data preparation                #
# ---------------------------------------------------------------------- #


# What is done within this script:

# We process future climate datasets from ISIMIP3b (Inter-Sectoral Impact Model 
# Intercomparison Project Phase 3b) for five different climate models
# and three different scenarios (ssp126, ssp370, ssp585). The climate data are 
# provided as daily outputs at a 0.5° spatial resolution. To match our target 
# temporal resolution, we aggregate daily temperature and humidity values into 
# monthly means and daily precipitation into monthly totals. In addition, we 
# process future LUH2 (Land Use Harmonization) land-use datasets for the same
# three scenarios. Provided at a 0.25° resolution, the rasters are aggregated
# by a factor of 2 to achieve a 0.5° spatial resolution and summarise them into 
# eight distinct land-use categories for each year.


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

# 2. Prepare future climate data -----------------------------------------------

# CHELSA MPI-ESM1-2-HR, GFDL-ESM4, IPSL-CM6A-LR, MRI-ESM2-0, UKESM1-0-LL
# https://data.isimip.org/search/tree/ISIMIP3b/InputData/climate/atmosphere/mpi-esm1-2-hr/ -> configure download
# Bounding box: South: 34 North: 75 West: -31 East: 60
# Extracted files from zip document
# 2020-2059


# Load in all data of the climate variables (pr, tas, tasmax, tasmin, hurs) 
# and calculate monthly averages (temperature, relative humidity) or monthly sums 
# (precipitation) for each future year and for there different climate forcing 
# scenarios per climate model

# Prepare path to data folder
datapath_Climate_data <- file.path("input_data/environmental_data/ISIMIP3b/Climate")

# Create a vector containing the months of a year
months <- str_pad(1:12, width = 2, pad = "0")

# Create a vector containing the years of interest
years <- 2020:2059

# Create a vector containing the five climatic target variables
clim_variables <- c("pr", "tas", "tasmax", "tasmin", "hurs")

# Create a vector containing the three different climate forcing scenarios
clim_scenario <- c("ssp126", "ssp370", "ssp585")

# Create a vector containing the different climate models
clim_models <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")

# Loop through the different year-month combinations for each variable,
# process climate data and write a raster for each year-month combination 
# containing all five climate variables (for 3 different climate scenarios
# based on 5 different climate models)

for (s in clim_scenario) { # Start of the loop over the three different climate forcing scenarios
  
  print(s)
  
  for (l in clim_models) { # Start the loop over the 5 different climate models
    
    print(l)
    
    for (y in years) { # Start of the loop over the different years
      
      print(y)
      
      for (m in months) { # Start of the loop over the different months
        
        print(m)
        
        # Check if climate file was already processed for respective month and year
        file_exists <- file.exists(paste0(datapath_Climate_data, "/",s,"/processed_data/",l,"/Climate_future_data_",m,"_",y,"_",s,".tif"))
        if (file_exists == FALSE) { # If file does not exist, start processing
          
          # Write a vector with the year-month combination
          year_month_select <- paste(y, m, sep = "-")
          
          for (c in clim_variables) { # Start of the loop over the five different climate variables
            
            print(c)
            
            # Adapt the file name depending on the model
            if (l %in% c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0")) { v <- "r1i1p1f1"
            } else if (l == "ukesm1-0-ll") { v <- "r1i1p1f2"
            }
            
            # Load climate files that contain data for the years from 2030 to 2070
            # for the respective variable
            Climate_data_rasters_1 <- terra::rast(paste0(datapath_Climate_data, "/",s,"/raw_data/",l,"/",l,"_",v,"_w5e5_",s,"_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2015_2020.nc"))
            Climate_data_rasters_2 <- terra::rast(paste0(datapath_Climate_data, "/",s,"/raw_data/",l,"/",l,"_",v,"_w5e5_",s,"_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2021_2030.nc"))
            Climate_data_rasters_3 <- terra::rast(paste0(datapath_Climate_data, "/",s,"/raw_data/",l,"/",l,"_",v,"_w5e5_",s,"_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2031_2040.nc"))
            Climate_data_rasters_4 <- terra::rast(paste0(datapath_Climate_data, "/",s,"/raw_data/",l,"/",l,"_",v,"_w5e5_",s,"_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2041_2050.nc"))
            Climate_data_rasters_5 <- terra::rast(paste0(datapath_Climate_data, "/",s,"/raw_data/",l,"/",l,"_",v,"_w5e5_",s,"_",c,"_lon-31.0to60.0lat34.0to75.0_daily_2051_2060.nc"))
            
            # Stack all rasters
            Climate_data_rasters <- c(Climate_data_rasters_1, Climate_data_rasters_2, Climate_data_rasters_3,
                                      Climate_data_rasters_4, Climate_data_rasters_5)
            
            # Extract the time information from the rasters, starting from the year 2020
            Climate_data_rasters_dates <- as.Date(time(Climate_data_rasters), origin = "2020-01-01")
            
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
          terra::writeRaster(Climate_rasters_processed, filename = paste0(datapath_Climate_data, "/",s,"/processed_data/",l,"/Climate_future_data_",m,"_",y,"_",s,".tif"), overwrite = TRUE)
          
          
          
        } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next month
        } # Close if condition
        
        
        
      } # Close the loop over all months of a year
      
    } # Close the loop over all years
    
  } # Close the loop over the 5 different climate models
  
  
} # Close of the loop over the climate forcing scenarios




#-------------------------------------------------------------------------------

# 3. Prepare future land use data ----------------------------------------------

# Download future land use scenarios (states at resolution of 0.25°) from LUH2 
# (https://luh.umd.edu/data.shtml) - v2f Release
# Three different future scenarios: RCP2.6 SSP1, RCP7.0 SSP3, RCP8.5 SSP5

# Prepare path to data folder
datapath_LandUse_data <- file.path("input_data/environmental_data/ISIMIP3b/LandUse")

# Create a vector containing the years of interest
years <- 2020:2059

# Create a vector containing the three different land use forcing scenarios
landuse_scenario <- c("ssp126", "ssp370", "ssp585")

# Create a vector containing the models used to generate forcing scenario
model <- c("IMAGE", "AIM", "MAGPIE")

for (s in landuse_scenario) {
  
  print(s)
  
  print("Load data")
  
  # Assign the models that were used to generate the respective forcing scenario
  if (s == "ssp126") { m <- "IMAGE"
  } else if (s == "ssp370") { m <- "AIM"
  } else if (s == "ssp585") { m <- "MAGPIE"
  }
  
  # Load in yearly land use data
  LandUse_data_rasters_future <- terra::rast(paste0(datapath_LandUse_data, "/",s,"/raw_data/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-",m,"-",s,"-2-1-f_gn_2015-2100.nc"))
  
  
  # Extract the rasters for each variable of interest
  print("forested primary land")
  LandUse_data_rasters_primf <- LandUse_data_rasters_future[[grep("primf", names(LandUse_data_rasters_future))]]
  
  print("non-forested primary land")
  LandUse_data_rasters_primn <- LandUse_data_rasters_future[[grep("primn", names(LandUse_data_rasters_future))]]
  
  print("potentially forested secondary land")
  LandUse_data_rasters_secdf <- LandUse_data_rasters_future[[grep("secdf", names(LandUse_data_rasters_future))]]
  
  print("potentially non-forested secondary land")
  LandUse_data_rasters_secdn <- LandUse_data_rasters_future[[grep("secdn", names(LandUse_data_rasters_future))]]
  
  print("managed pasture")
  LandUse_data_rasters_pastr <- LandUse_data_rasters_future[[grep("pastr", names(LandUse_data_rasters_future))]]
  
  print("rangeland")
  LandUse_data_rasters_range <- LandUse_data_rasters_future[[grep("range", names(LandUse_data_rasters_future))]]
  
  print("C3 annual crops")
  LandUse_data_rasters_crop1 <- LandUse_data_rasters_future[[grep("c3ann", names(LandUse_data_rasters_future))]]
  
  print("C4 annual crops")
  LandUse_data_rasters_crop2 <- LandUse_data_rasters_future[[grep("c4ann", names(LandUse_data_rasters_future))]]
  
  print("C3 perennial crops")
  LandUse_data_rasters_crop3 <- LandUse_data_rasters_future[[grep("c3per", names(LandUse_data_rasters_future))]]
  
  print("C4 perennial crops")
  LandUse_data_rasters_crop4 <- LandUse_data_rasters_future[[grep("c4per", names(LandUse_data_rasters_future))]]
  
  print("C3 nitrogen-fixing crops")
  LandUse_data_rasters_crop5 <- LandUse_data_rasters_future[[grep("c3nfx", names(LandUse_data_rasters_future))]]
  
  print("urban")
  LandUse_data_rasters_urban <- LandUse_data_rasters_future[[grep("urban", names(LandUse_data_rasters_future))]]
  
  
  for (y in years) { # Start of the loop over the years of interest
    
    print(y)
    
    # Check if land use file was already processed for respective year
    file_exists <- file.exists(paste0(datapath_LandUse_data, "/",s,"/processed_data/LandUse_future_data_",y,"_",s,".tif"))
    if (file_exists == FALSE) { # If file does not exist, start processing
      
      # Extract the raster of the respective year for each land use variable,
      # aggregate the rasters by a factor of 2 to change the resolution to 0.5°,
      # crop it to the European extent and mask it to the terrestrial European
      # area
      print("forested primary land")
      LandUse_data_raster_primf_year <- LandUse_data_rasters_primf[[grep(y, time(LandUse_data_rasters_primf))]]
      LandUse_data_raster_primf_year <- terra::aggregate(LandUse_data_raster_primf_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_primf_year <- terra::crop(LandUse_data_raster_primf_year, europe_mask)
      LandUse_data_raster_primf_year <- terra::mask(LandUse_data_raster_primf_year, europe_mask)
      
      print("non-forested primary land")
      LandUse_data_raster_primn_year <- LandUse_data_rasters_primn[[grep(y, time(LandUse_data_rasters_primn))]]
      LandUse_data_raster_primn_year <- terra::aggregate(LandUse_data_raster_primn_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_primn_year <- terra::crop(LandUse_data_raster_primn_year, europe_mask)
      LandUse_data_raster_primn_year <- terra::mask(LandUse_data_raster_primn_year, europe_mask)
      
      print("potentially forested secondary land")
      LandUse_data_raster_secdf_year <- LandUse_data_rasters_secdf[[grep(y, time(LandUse_data_rasters_secdf))]]
      LandUse_data_raster_secdf_year <- terra::aggregate(LandUse_data_raster_secdf_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_secdf_year <- terra::crop(LandUse_data_raster_secdf_year, europe_mask)
      LandUse_data_raster_secdf_year <- terra::mask(LandUse_data_raster_secdf_year, europe_mask)
      
      print("potentially non-forested secondary land")
      LandUse_data_raster_secdn_year <- LandUse_data_rasters_secdn[[grep(y, time(LandUse_data_rasters_secdn))]]
      LandUse_data_raster_secdn_year <- terra::aggregate(LandUse_data_raster_secdn_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_secdn_year <- terra::crop(LandUse_data_raster_secdn_year, europe_mask)
      LandUse_data_raster_secdn_year <- terra::mask(LandUse_data_raster_secdn_year, europe_mask)
      
      print("managed pasture")
      LandUse_data_raster_pastr_year <- LandUse_data_rasters_pastr[[grep(y, time(LandUse_data_rasters_pastr))]]
      LandUse_data_raster_pastr_year <- terra::aggregate(LandUse_data_raster_pastr_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_pastr_year <- terra::crop(LandUse_data_raster_pastr_year, europe_mask)
      LandUse_data_raster_pastr_year <- terra::mask(LandUse_data_raster_pastr_year, europe_mask)
      
      print("rangeland")
      LandUse_data_raster_range_year <- LandUse_data_rasters_range[[grep(y, time(LandUse_data_rasters_range))]]
      LandUse_data_raster_range_year <- terra::aggregate(LandUse_data_raster_range_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_range_year <- terra::crop(LandUse_data_raster_range_year, europe_mask)
      LandUse_data_raster_range_year <- terra::mask(LandUse_data_raster_range_year, europe_mask)
      
      print("C3 annual crops")
      LandUse_data_raster_crop1_year <- LandUse_data_rasters_crop1[[grep(y, time(LandUse_data_rasters_crop1))]]
      LandUse_data_raster_crop1_year <- terra::aggregate(LandUse_data_raster_crop1_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_crop1_year <- terra::crop(LandUse_data_raster_crop1_year, europe_mask)
      LandUse_data_raster_crop1_year <- terra::mask(LandUse_data_raster_crop1_year, europe_mask)
      
      print("C4 annual crops")
      LandUse_data_raster_crop2_year <- LandUse_data_rasters_crop2[[grep(y, time(LandUse_data_rasters_crop2))]]
      LandUse_data_raster_crop2_year <- terra::aggregate(LandUse_data_raster_crop2_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_crop2_year <- terra::crop(LandUse_data_raster_crop2_year, europe_mask)
      LandUse_data_raster_crop2_year <- terra::mask(LandUse_data_raster_crop2_year, europe_mask)
      
      print("C3 perennial crops")
      LandUse_data_raster_crop3_year <- LandUse_data_rasters_crop3[[grep(y, time(LandUse_data_rasters_crop3))]]
      LandUse_data_raster_crop3_year <- terra::aggregate(LandUse_data_raster_crop3_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_crop3_year <- terra::crop(LandUse_data_raster_crop3_year, europe_mask)
      LandUse_data_raster_crop3_year <- terra::mask(LandUse_data_raster_crop3_year, europe_mask)
      
      print("C4 perennial crops")
      LandUse_data_raster_crop4_year <- LandUse_data_rasters_crop4[[grep(y, time(LandUse_data_rasters_crop4))]]
      LandUse_data_raster_crop4_year <- terra::aggregate(LandUse_data_raster_crop4_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_crop4_year <- terra::crop(LandUse_data_raster_crop4_year, europe_mask)
      LandUse_data_raster_crop4_year <- terra::mask(LandUse_data_raster_crop4_year, europe_mask)
      
      print("C3 nitrogen-fixing crops")
      LandUse_data_raster_crop5_year <- LandUse_data_rasters_crop5[[grep(y, time(LandUse_data_rasters_crop5))]]
      LandUse_data_raster_crop5_year <- terra::aggregate(LandUse_data_raster_crop5_year, fact = 2, fun = mean, na.rm = TRUE)
      LandUse_data_raster_crop5_year <- terra::crop(LandUse_data_raster_crop5_year, europe_mask)
      LandUse_data_raster_crop5_year <- terra::mask(LandUse_data_raster_crop5_year, europe_mask)
      
      print("urban land")
      LandUse_data_raster_urban_year <- LandUse_data_rasters_urban[[grep(y, time(LandUse_data_rasters_urban))]]
      LandUse_data_raster_urban_year <- terra::aggregate(LandUse_data_raster_urban_year, fact = 2, fun = mean, na.rm = TRUE)
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
      terra::writeRaster(LandUse_rasters_processed, filename = paste0(datapath_LandUse_data, "/",s,"/processed_data/LandUse_future_data_",y,"_",s,".tif"), overwrite = TRUE)
      
      
    } else if (file_exists == TRUE) { print("already done") # If file already exists, start with next year
    } # Close if condition
    
  } # End of the loop over all years of interest
  
  
} # End of the loop over the three land use forcing scenarios




#-------------------------------------------------------------------------------

# 4. Prepare quantification of climate and land-use changes --------------------

# Quantify the projected changes from the baseline (2010s) to the future (2050s) 
# of different climate and land-use variables per SSP scenario
# This is done to evaluate the implications for climate and land-use change 
# in Europe according to the three different SSP pathways

# (a) Climate - temperature and precipitation ----------------------------------

# Prepare path to data folders
datapath_Climate_data_future <- file.path("input_data/environmental_data/ISIMIP3b/Climate")
datapath_Climate_data_baseline <- file.path("input_data/environmental_data/ISIMIP3a/Climate")

# Create a vector containing the three different climate forcing scenarios
clim_scenario <- c("ssp126", "ssp370", "ssp585")

# Create a vector containing the different climate models
clim_models <- c("gfdl-esm4", "ipsl-cm6a-lr", "mpi-esm1-2-hr", "mri-esm2-0", "ukesm1-0-ll")

# Prepare a data frame to store the quantified climate changes per SSP and
# climate variable
climate_change_quantification <- data.frame()

# Prepare baseline (2010s) 
baseline_files <- list.files(path = paste0(datapath_Climate_data_baseline, "/processed_data/"), pattern = "201", full.names = TRUE) # List files
baseline_stack <- terra::rast(baseline_files) # Create raster stacks

baseline_stack_pr <- baseline_stack[[grep("^pr$", names(baseline_stack))]] # Precipitation stack
baseline_stack_tas <- baseline_stack[[grep("^tas$", names(baseline_stack))]] # Mean temperature stack

baseline_tas_mean <- mean(values(baseline_stack_tas), na.rm = TRUE) # Compute the long-term mean across all months in the 2010s

 # Compute decadal mean annual precipitation
years <- as.numeric(str_extract(basename(baseline_files), "\\d{4}"))
baseline_pr_annual <- tapp(baseline_stack_pr, years, sum)
baseline_pr_mean <- mean(values(baseline_pr_annual), na.rm = TRUE)

# Prepare future (2050s)
# Loop through the different SSP pathways and the 5 different climate models
# Calculate the decadal temperature mean and the decadal mean annual precipitation
# for each climate model. These model-specific decadal values are then compared
# to a common baseline period (2010s) to compute climate change signals.
# Finally, the resulting temperature and precipitation changes are averaged
# across the five climate models to obtain a ensemble mean for each SSP scenario.

for (s in clim_scenario) { # Start of the loop over the three different forcing scenarios
  
  print(s)
  
  tas_models <- c()
  pr_models <- c()
  
  for (l in clim_models) { # Start the loop over the 5 different climate models
    
    print(l)
    
    # List files
    future_files <- list.files(path = paste0(datapath_Climate_data_future, "/",s,"/processed_data/",l,"/"), pattern = "205", full.names = TRUE)
    
    # Create raster stacks
    future_stack <- terra::rast(future_files)
    
    # Create different raster stacks for precipitation and mean temperature
    future_stack_pr <- future_stack[[grep("^pr$", names(future_stack))]]
    future_stack_tas <- future_stack[[grep("^tas$", names(future_stack))]]
    
    # Compute the long-term mean across all months in the 2050s
    future_tas_mean <- mean(values(future_stack_tas), na.rm = TRUE)
    
    # Compute decadal mean annual precipitation
    years <- as.numeric(str_extract(basename(future_files), "\\d{4}"))
    future_pr_annual <- tapp(future_stack_pr, years, sum)
    future_pr_mean <- mean(values(future_pr_annual), na.rm = TRUE)
    
    # Quantify projected changes from baseline to future
    tas_change_m <- future_tas_mean - baseline_tas_mean
    pr_change_m <- future_pr_mean - baseline_pr_mean
    
    # Store results per climate model
    tas_models <- c(tas_models, tas_change_m)
    pr_models <- c(pr_models, pr_change_m)
    
  } # Close loop over the 5 climate models
  
  # Create an ensemble over the different climate models under each SSP
  tas_change <- mean(tas_models, na.rm = TRUE)
  pr_change  <- mean(pr_models, na.rm = TRUE)
  pr_change_pct <- (pr_change / baseline_pr_mean) * 100
  
  # Store SSP results
  climate_change_quantification <- rbind(climate_change_quantification, 
                                         data.frame(SSP = s, tas_change = tas_change,
                                                    pr_change = pr_change,
                                                    pr_change_pct = pr_change_pct))

} # Close loop over the 3 SSP pathways

print(climate_change_quantification)
  
    
    
    
# (b) Land use -----------------------------------------------------------------

# Prepare path to data folder
datapath_LandUse_data_future <- file.path("input_data/environmental_data/ISIMIP3b/LandUse")
datapath_LandUse_data_baseline <- file.path("input_data/environmental_data/ISIMIP3a/LandUse")

# Create a vector containing the three different environmental forcing scenarios
landuse_scenario <- c("ssp126", "ssp370", "ssp585")

# Prepare a data frame to store the quantified land-use changes per SSP and
# land-use variable
landuse_change_quantification <- data.frame()

# Prepare baseline (2010s)
baseline_files <- list.files(path = paste0(datapath_LandUse_data_baseline, "/processed_data/"), pattern = "201.*\\.tif$", full.names = TRUE) # List files
baseline_stack <- terra::rast(baseline_files) # Create raster stacks

baseline_stack_primary_forest <- baseline_stack[[grep("^primary_forest$", names(baseline_stack))]] # Primary forest stack
baseline_stack_primary_openland <- baseline_stack[[grep("^primary_openland$", names(baseline_stack))]] # Secondary open land stack
baseline_stack_secondary_forest <- baseline_stack[[grep("^secondary_forest$", names(baseline_stack))]] # Secondary forest stack
baseline_stack_secondary_openland <- baseline_stack[[grep("^secondary_openland$", names(baseline_stack))]] # Secondary open land stack
baseline_stack_pasture <- baseline_stack[[grep("^pasture$", names(baseline_stack))]] # Pasture stack
baseline_stack_rangeland <- baseline_stack[[grep("^rangeland$", names(baseline_stack))]] # Rangeland stack
baseline_stack_cropland <- baseline_stack[[grep("^cropland$", names(baseline_stack))]] # Cropland stack
baseline_stack_urban <- baseline_stack[[grep("^urban$", names(baseline_stack))]] # Urban stack

# Compute the long-term mean across all years in the 2010s
baseline_primary_forest_mean <- mean(values(baseline_stack_primary_forest), na.rm = TRUE)
baseline_primary_openland_mean <- mean(values(baseline_stack_primary_openland), na.rm = TRUE)
baseline_secondary_forest_mean <- mean(values(baseline_stack_secondary_forest), na.rm = TRUE)
baseline_secondary_openland_mean <- mean(values(baseline_stack_secondary_openland), na.rm = TRUE)
baseline_pasture_mean <- mean(values(baseline_stack_pasture), na.rm = TRUE)
baseline_rangeland_mean <- mean(values(baseline_stack_rangeland), na.rm = TRUE)
baseline_cropland_mean <- mean(values(baseline_stack_cropland), na.rm = TRUE)
baseline_urban_mean <- mean(values(baseline_stack_urban), na.rm = TRUE)


# Prepare future (2050s)
# Changes in land-use composition between the baseline period (2010s) and future
# projections (2050s), under three different SSP scenarios are quantified.
# For each land-use category, the difference between calculated baseline values 
# as the mean across the baseline period and future values as the mean across 
# the 2050s are extracted, resulting in scenario-specific shift in land-use
# fractions.
for (s in clim_scenario) { # Start of the loop over the three different forcing scenarios
  
  print(s)
  
  # List files
  future_files <- list.files(path = paste0(datapath_LandUse_data_future, "/",s,"/processed_data/"), pattern = "205.*\\.tif$", full.names = TRUE)
  
  # Create raster stacks
  future_stack <- terra::rast(future_files)
  
  # Create different raster stacks for the different land-use categories
  future_stack_primary_forest <- future_stack[[grep("^primary_forest$", names(future_stack))]] # Primary forest stack
  future_stack_primary_openland <- future_stack[[grep("^primary_openland$", names(future_stack))]] # Secondary open land stack
  future_stack_secondary_forest <- future_stack[[grep("^secondary_forest$", names(future_stack))]] # Secondary forest stack
  future_stack_secondary_openland <- future_stack[[grep("^secondary_openland$", names(future_stack))]] # Secondary open land stack
  future_stack_pasture <- future_stack[[grep("^pasture$", names(future_stack))]] # Pasture stack
  future_stack_rangeland <- future_stack[[grep("^rangeland$", names(future_stack))]] # Rangeland stack
  future_stack_cropland <- future_stack[[grep("^cropland$", names(future_stack))]] # Cropland stack
  future_stack_urban <- future_stack[[grep("^urban$", names(future_stack))]] # Urban stack
  
  # Compute the long-term mean across all years in the 2050s
  future_primary_forest_mean <- mean(values(future_stack_primary_forest), na.rm = TRUE)
  future_primary_openland_mean <- mean(values(future_stack_primary_openland), na.rm = TRUE)
  future_secondary_forest_mean <- mean(values(future_stack_secondary_forest), na.rm = TRUE)
  future_secondary_openland_mean <- mean(values(future_stack_secondary_openland), na.rm = TRUE)
  future_pasture_mean <- mean(values(future_stack_pasture), na.rm = TRUE)
  future_rangeland_mean <- mean(values(future_stack_rangeland), na.rm = TRUE)
  future_cropland_mean <- mean(values(future_stack_cropland), na.rm = TRUE)
  future_urban_mean <- mean(values(future_stack_urban), na.rm = TRUE)
  
  # Quantify change from baseline period to future period
  primary_forest_change <- future_primary_forest_mean - baseline_primary_forest_mean
  primary_openland_change <- future_primary_openland_mean - baseline_primary_openland_mean
  secondary_forest_change <- future_secondary_forest_mean - baseline_secondary_forest_mean
  secondary_openland_change <- future_secondary_openland_mean - baseline_secondary_openland_mean
  pasture_change <- future_pasture_mean - baseline_pasture_mean
  rangeland_change <- future_rangeland_mean - baseline_rangeland_mean
  cropland_change <- future_cropland_mean - baseline_cropland_mean
  urban_change <- future_urban_mean - baseline_urban_mean
  
  # Store SSP results
  landuse_change_quantification <- rbind(landuse_change_quantification, 
                                         data.frame(SSP = s, 
                                                    primary_forest_change = primary_forest_change,
                                                    primary_openland_change = primary_openland_change,
                                                    secondary_forest_change = secondary_forest_change,
                                                    secondary_openland_change = secondary_openland_change,
                                                    pasture_change = pasture_change,
                                                    rangeland_change = rangeland_change,
                                                    cropland_change = cropland_change,
                                                    urban_change = urban_change))
  
} # Close the loop over the different SSP pathways
  
print(landuse_change_quantification)

# save the data frame containing the results
save(climate_change_quantification, landuse_change_quantification, file = "input_data/environmental_data/climate_landuse_change_quantification.RData")

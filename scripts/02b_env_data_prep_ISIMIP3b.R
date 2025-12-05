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




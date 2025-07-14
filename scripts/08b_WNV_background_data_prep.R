# ZOE project 
# Disease phenology analysis of West Virus Fever in Europe

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                   08b. WNV infection data preparation                  #
# ---------------------------------------------------------------------- #




# Load needed packages
library(terra)
library(sf)
library(sfheaders)
library(tidyverse)

# Load needed objects
source("scripts/00_functions.R") # Get the thin function
nuts_3_raster_mask <- terra::rast("input_data/spatial_data/nuts_3_raster_mask.tif") # Background mask of EU/EEA countries in a 50 km resolution
load("output_data/data/WNV_occurrences_cleaned.RData") # Infection occurrence data
r_curr_preds_clim_landuse <- terra::rast("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif") # Ensemble predictions of occurrence probability (under factual climate and land use) of the main vector species Culex pipiens
load("output_data/data/nuts_3_WNV.RData") # Data frame containing all infection points, also of municipalities consisting of > 5 cells


# Prepare path to data folder
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Create a sequence of dates with monthly steps, that are temporally
# covered by the environmental data and infection data
start_year <- min(WNV_occurrences_cleaned$year) # Find the year of earliest observation
start_date <- as.Date(paste0(start_year,"-01-01")) # Define start date
end_date <- as.Date("2019-12-01") # Define end date
date_sequence_month <- seq.Date(from = start_date, to = end_date, by = "month") # create a monthly sequence
date_sequence <- format(date_sequence_month, "%m/%Y") # Extract year and month from dates

# Format the months within the data frames of cleaned infection data
WNV_occurrences_cleaned$month <- sprintf("%02d", as.numeric(WNV_occurrences_cleaned$month))
nuts_3_WNV$DateUsedForStatisticsMonth <- sprintf("%02d", as.numeric(nuts_3_WNV$DateUsedForStatisticsMonth))


# Prepare the EU/EEA mask to only have values of 1 or NA
eu_eea_mask <- nuts_3_raster_mask
eu_eea_mask[] <- ifelse(!is.na(eu_eea_mask[]), 1, NA)

# Create a vector that contains the two-letter iso codes of the EU/EEA member states
# that send data from their surveillance systems to ECDC (If not identical with the
# NUTS3 country code, also add that one; Examples: Greece (GR, EL), United Kingdom (GB, UK))
ECDC_reporting_countries <- c("AT", "BE", "BG", "HR", "CY", "CZ", "DK", "EE", "FI",
                              "FR", "DE", "GB", "UK", "EL", "GR", "HU", "IE", "IT", "LV", "LT", 
                              "LU", "MT", "NL", "PL", "PT", "RO", "SK", "SI", "LI", "ES", 
                              "SE", "IS", "LI", "NO")

# Select only the relevant columns, order them, and rename the id column
WNV_occurrences_cleaned <- WNV_occurrences_cleaned[c("lon", "lat", "occ", "year", "month")]

# Create an empty data frame to store the thinned presences and absence data
WNV_occ_env <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(WNV_occ_env) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "C_pipiens")



# Loop over all dates in the sequence, 
# remove duplicates in 50 km² cells, generate background data,
# and match the presence and background data with the time-specific environmental data.
for (d in date_sequence) { # Start of the loop over all dates
  
  # Extract month and year from date
  m <- substr(d, 1, 2)
  print(m)
  y <- substr(d, 4, 7)
  print(y)
  
  # Subset the infection occurrence data frame by each date (month and year)
  subset_year_month <- subset(WNV_occurrences_cleaned, WNV_occurrences_cleaned$year == y & WNV_occurrences_cleaned$month == m)
  
  
  if (nrow(subset_year_month) > 0) { # Just continue with preparation process if occurrences are available for the date
    
    print("start with data prep")
    
    
#-------------------------------------------------------------------------------
    
# 1. Remove duplicates in cells ------------------------------------------------
    
    # Select coordinates of the occurrences
    occ_coords <- as.data.frame(subset_year_month[, c("lon", "lat")])
    
    # Extract an ID per cell
    cellnumbers <- terra::extract(eu_eea_mask, occ_coords, cells = TRUE)
    
    # Only keep cells that are not duplicated
    occ_coords <- occ_coords[!duplicated(cellnumbers[,"cell"]),]
    
    
#-------------------------------------------------------------------------------
    
# 2. Generation of background data ---------------------------------------------
    
    print("start with background generation")
    
    # Retrieve all NUTS3 municipalities that reported infections for that month-year
    # combination (also NUTS3 municipalities consisting of > 5 cells)
    subset_year_month_all <- nuts_3_WNV[nuts_3_WNV$DateUsedForStatisticsYear == y & 
                                           nuts_3_WNV$DateUsedForStatisticsMonth == m, ]
    
    subset_year_month_all_nuts <- unique(subset_year_month_all$new_municipality)
    
    
    # Create absence points in all cells of NUTS3 municipalities where no infection
    # was reported even though the country is part of the EU/EEA mandatory surveillance
    # data reporting system
    # First, create points in all cells corresponding to a NUTS3 municipality
    points_all_cells <- as.points(nuts_3_raster_mask, values = TRUE)
    
    # Extract NUTS3 values of the raster
    points_all_cells$NUTS_3 <- values(nuts_3_raster_mask)
    
    # Identify NUTS3 codes that do not have recorded infections for that month-year
    # combination
    all_nuts3 <- unique(points_all_cells$NUTS_ID)
    non_reporting_nuts3 <- setdiff(all_nuts3, subset_year_month_all_nuts)
    
    # Extract first two letters (country codes) of NUTS3 regions
    points_all_cells$country_code <- substr(points_all_cells$NUTS_ID, 1, 2)
    
    # Only keep points where NUTS3 is non-reporting and country is in the EU/EEA
    # mandatory surveillance system
    abs_points <- points_all_cells[points_all_cells$NUTS_ID %in% non_reporting_nuts3 &
                                     points_all_cells$country_code %in% ECDC_reporting_countries, ]
    
    # Extract the coordinates of the presence points Spatvector
    abs_coords <- terra::geom(abs_points)
    abs_coords <- data.frame(lon = abs_coords[, "x"], lat = abs_coords[, "y"])
    
    
    # Add information on presence and absence
    occ_coords$occ <- 1
    abs_coords$occ <- 0
    
    
    
#-------------------------------------------------------------------------------
    
# 3. Spatial thinning of the presence and background data ----------------------
    
    print("spatial thinning")
    
    # Thinning of presences
    # Transform data frame into sf object to use in thin function
    occ_coords_sf <- st_as_sf(occ_coords, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
    
    # Using the thin function with a thinning distance of 50 km (checkerboard pattern)
    occ_coords_thinned <- thin(occ_coords_sf, thin_dist = 50000, runs = 1, ncores = 1)
    
    # Merge data frames to only retained thinned presences and background
    occ_coords_thinned <- merge(occ_coords_thinned, occ_coords, by = c("lon", "lat"))
    
    
    # Thinning of background data
    # Transform data frame into sf object to use in thin function
    abs_coords_sf <- st_as_sf(abs_coords, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
    
    # Using the thin function with a thinning distance of 50 km (checkerboard pattern)
    abs_coords_thinned <- thin(abs_coords_sf, thin_dist = 50000, runs = 1, ncores = 1)
    
    # Merge data frames to only retained thinned presences and background
    abs_coords_thinned <- merge(abs_coords_thinned, abs_coords, by = c("lon", "lat"))
    
    
    # Join presence and background data
    WNV_occ_thinned <- rbind(occ_coords_thinned, abs_coords_thinned)
    
    # Add the year and month as information in columns
    WNV_occ_thinned$year <- y
    WNV_occ_thinned$month <- m
    
    
    
#-------------------------------------------------------------------------------
    
# 4. Join with environmental data  ---------------------------------------------
    
    print("matching env. data")
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    Species_data <- r_curr_preds_clim_landuse[[paste0(m,"/",y)]] # Occurrence probabilities of the vector species Culex pipiens (under factual climate and land use)
    
    # Make sure the extents of climate matches the
    # species prediction data
    Climate_data <- terra::crop(Climate_data, Species_data)
    
    # Make sure the species data has the correct name
    names(Species_data) <- "C_pipiens"
    
    # Stack the environmental data
    env_data <- c(Climate_data, Species_data)
    
    # Extract the environmental values per occurrence cell
    WNV_occ_env_date <- cbind(WNV_occ_thinned, terra::extract(x = env_data, y = WNV_occ_thinned[,c('lon','lat')]))
    
    # Drop NA for the environmental variables 
    WNV_occ_env_date <- WNV_occ_env_date %>% drop_na()
    
    # Check for duplicates
    duplicated(WNV_occ_env_date$ID)
    
    # Only retain non-duplicated cells
    WNV_occ_env_date <- WNV_occ_env_date[!duplicated(WNV_occ_env_date$ID),]
    
    # Add the thinned presences and background data belonging to the specific year and month to the 
    # prepared results data frame
    WNV_occ_env <- rbind(WNV_occ_env, WNV_occ_env_date)
    
    
  } else if (nrow(subset_year_month) == 0) { print("no data available for the year-month combination")
  } # End of if-condition
  
} # End of loop over dates



# Get a summary of presence and background data numbers
table(WNV_occ_env$occ) # 1: 629, 0: 67646
print(table(WNV_occ_env$month[WNV_occ_env$occ == 1]))

# Save the resulting data frame, containing thinned presence and background data,
# joined with the respective environmental data of year and month
save(WNV_occ_env, file = "output_data/data/WNV_occ_env.RData")



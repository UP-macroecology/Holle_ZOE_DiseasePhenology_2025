# ZOE project 
# Disease phenology analysis of TBEV in Europe 

# ---------------------------------------------------------------------- #
#                    07a. Background data preparation                    #
# ---------------------------------------------------------------------- #


# Load needed packages
library(terra)
library(sf)
library(sfheaders)
library(tidyverse)
library(countrycode)

# Load needed objects
source("scripts/00_functions.R") # Get the thin function
nuts_3_raster_mask <- terra::rast("input_data/spatial_data/nuts_3_raster_mask.tif") # Background mask of EU/EEA countries in a 50 km resolution
load("output_data/data/TBEV_occurrences_cleaned.RData") # Infection occurrence data
load("output_data/data/I_ricinus_occurrences_cleaned.RData") # The presence information of tick species Ixodes ricinus
r_curr_preds_clim_landuse <- terra::rast("output_data/results/I_ricinus_preds_clim_landuse_ens_1970_2019.tif") # Ensemble predictions of occurrence probability of Ixodes ricinus
load("output_data/validation/I_ricinus_validation.RData") # Performance measures with presence-absence threshold (MaxSSS)

# Prepare path to data folder
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Create a sequence of dates with monthly steps, that are temporally
# covered by the environmental data and infection data
start_year <- min(TBEV_occurrences_cleaned$year) # Find the year of earliest observation
start_date <- as.Date(paste0(start_year,"-01-01")) # Define start date
end_date <- as.Date("2019-12-01") # Define end date
date_sequence_month <- seq.Date(from = start_date, to = end_date, by = "month") # create a monthly sequence
date_sequence <- format(date_sequence_month, "%m/%Y") # Extract year and month from dates

# Format the months within the data frames of cleaned infection and species occurrences
TBEV_occurrences_cleaned$month <- sprintf("%02d", as.numeric(TBEV_occurrences_cleaned$month))
I_ricinus_occurrences_cleaned$month <- sprintf("%02d", as.numeric(I_ricinus_occurrences_cleaned$month))

# Convert the country column into a two-letter ISO code for all entries from VectorMap
I_ricinus_occurrences_cleaned$country[I_ricinus_occurrences_cleaned$database == "VectorMap"] <- countrycode(I_ricinus_occurrences_cleaned$country[I_ricinus_occurrences_cleaned$database == "VectorMap"], 
                                                                                                            origin = "country.name", destination = "iso2c")

# Find the threshold that maximises the sum of sensitivity and specificity
threshold <- comp_perf["mean_prob", "thresh"]

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
TBEV_occurrences_cleaned <- TBEV_occurrences_cleaned[c("lon", "lat", "occ", "year", "month")]

# Create an empty data frame to store the thinned presences and background data
TBEV_occ_env <- data.frame(matrix(ncol = 20, nrow = 0))
colnames(TBEV_occ_env) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "primary_forest", "primary_openland", "secondary_forest", 
                            "secondary_openland", "pasture", "rangeland", "cropland", "urban", "I_ricinus")



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
  subset_year_month <- subset(TBEV_occurrences_cleaned, TBEV_occurrences_cleaned$year == y & TBEV_occurrences_cleaned$month == m)
  
  
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
    
    # Subset the presence information of associated vector species Ixodes ricinus
    # for the respective year and month
    subset_year_month_species <- subset(I_ricinus_occurrences_cleaned, I_ricinus_occurrences_cleaned$year == y & I_ricinus_occurrences_cleaned$month == m)
    
    # Remove occurrences of Ixodes ricinus that do not occur in countries with
    # surveillance data reporting to ECDC
    subset_year_month_species <- subset_year_month_species[subset_year_month_species$country %in% ECDC_reporting_countries, ]
    
    # Remove duplicates within a cell for vector species
    # Select coordinates of the occurrences
    occ_coords_species <- as.data.frame(subset_year_month_species[, c("lon", "lat")])
    
    # Extract an ID per cell
    cellnumbers_species <- terra::extract(eu_eea_mask, occ_coords_species, cells = TRUE)
    
    # Only keep cells that are not duplicated
    occ_coords_species <- occ_coords_species[!duplicated(cellnumbers_species[,"cell"]),]
    
    # Vectorise species and disease occurrences
    occ_coords_TBEV_vect <- terra::vect(occ_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
    occ_coords_species_vect <- terra::vect(occ_coords_species, geom = c("lon", "lat"), crs = "EPSG:4326")
    
    # Rasterise vector species and infection occurrence points
    occ_coords_TBEV_r <- terra::rasterize(occ_coords_TBEV_vect, eu_eea_mask, field = 1, background = NA)
    occ_coords_species_r <- terra::rasterize(occ_coords_species_vect, eu_eea_mask, field = 1, background = NA)
    
    # Identify cells with species presences but no infection occurrence
    abs_raster <- occ_coords_species_r == 1 & is.na(occ_coords_TBEV_r)
    
    # Convert raster cells back to points
    abs_points <- as.points(abs_raster, values = TRUE, na.rm = TRUE)
    
    # Extract the coordinates of the presence points Spatvector
    abs_coords <- terra::geom(abs_points)
    abs_coords <- data.frame(lon = abs_coords[, "x"], lat = abs_coords[, "y"])
    
    # Add information on presence
    occ_coords$occ <- 1
    
    
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
    
    # Add the year and month as information in columns
    occ_coords_thinned$year <- y
    occ_coords_thinned$month <- m
    
    
    
    if (nrow(abs_coords) >= 1) { # If there is more than one generated absence, start the thinning for background data
      
      # Thinning of background data
      # Add information on background
      abs_coords$occ <- 0
      
      # Transform data frame into sf object to use in thin function
      abs_coords_sf <- st_as_sf(abs_coords, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
      
      # Using the thin function with a thinning distance of 50 km (checkerboard pattern)
      abs_coords_thinned <- thin(abs_coords_sf, thin_dist = 50000, runs = 1, ncores = 1)
      
      # Merge data frames to only retained thinned presences and background
      abs_coords_thinned <- merge(abs_coords_thinned, abs_coords, by = c("lon", "lat"))
      
      # Add the year and month as information in columns
      abs_coords_thinned$year <- y
      abs_coords_thinned$month <- m
      
      
    } else if(nrow(abs_coords) == 0) { print("no background available for thinning")
      
      # Rename abs_coords
      abs_coords_thinned <- abs_coords
      
    } # End of if-condition
    
    
    # Join presence and background data (if available)
    TBEV_occ_thinned <- rbind(occ_coords_thinned, abs_coords_thinned)
    
    
    
    
#-------------------------------------------------------------------------------
    
# 4. Join with environmental data  ---------------------------------------------
    
    print("matching env. data")
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    Species_data <- r_curr_preds_clim_landuse[[paste0(m,"/",y)]] # Occurrence probabilities of the vector species Ixodes ricinus
    
    # Make sure the extents of climate and land use data matches the
    # species prediction data
    Climate_data <- terra::crop(Climate_data, Species_data)
    LandUse_data <- terra::crop(LandUse_data, Species_data)
    
    # Make sure the species data has the correct name
    names(Species_data) <- "I_ricinus"
    
    # Stack the environmental data
    env_data <- c(Climate_data, LandUse_data, Species_data)
    
    # Extract the environmental values per occurrence cell
    TBEV_occ_env_date <- cbind(TBEV_occ_thinned, terra::extract(x = env_data, y = TBEV_occ_thinned[,c('lon','lat')]))
    
    # Drop NA for the environmental variables 
    TBEV_occ_env_date <- TBEV_occ_env_date %>% drop_na()
    
    # Check for duplicates
    duplicated(TBEV_occ_env_date$ID)
    
    # Only retain non-duplicated cells
    TBEV_occ_env_date <- TBEV_occ_env_date[!duplicated(TBEV_occ_env_date$ID),]
    
    # Set all occurrence probability values below a threshold of the maxSSS for the 
    # generated background data points 
    TBEV_occ_env_date$I_ricinus[TBEV_occ_env_date$occ == 0 & TBEV_occ_env_date$I_ricinus < threshold] <- 0
    
    # Add the thinned presences and background data belonging to the specific year and month to the 
    # prepared results data frame
    TBEV_occ_env <- rbind(TBEV_occ_env, TBEV_occ_env_date)
    
    
  } else if (nrow(subset_year_month) == 0) { print("no data available for the year-month combination")
  } # End of if-condition
  
} # End of loop over dates





# Get a summary of presence and background data numbers
table(TBEV_occ_env$occ) # 1:2514, 0: 1491
print(table(TBEV_occ_env$month[TBEV_occ_env$occ == 1]))

# Save the resulting data frame, containing thinned presence and background data,
# joined with the respective environmental data of year and month
save(TBEV_occ_env, file = "output_data/data/TBEV_occ_env_vector.RData")


# Map the thinned presences and background data
# library(maps)
maps::map('world',xlim=c(-31,40), ylim=c(34,72))
points(TBEV_occ_env$lon[TBEV_occ_env$occ == 0], TBEV_occ_env$lat[TBEV_occ_env$occ == 0], col='steelblue4',  pch=19, cex = 0.5)
points(TBEV_occ_env$lon[TBEV_occ_env$occ == 1], TBEV_occ_env$lat[TBEV_occ_env$occ == 1], col='goldenrod',  pch=19, cex = 0.5)




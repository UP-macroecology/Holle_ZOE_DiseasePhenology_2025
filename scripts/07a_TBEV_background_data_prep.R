# ZOE project 
# Disease phenology analysis of TBEV in Europe 

# ---------------------------------------------------------------------- #
#                    03c. Background data preparation                    #
# ---------------------------------------------------------------------- #


# Load needed packages
library(terra)
library(sf)
library(sfheaders)
library(tidyverse)

# Load needed objects
source("scripts/00_functions.R") # Get the thin function
europe_mask <- terra::rast(paste0("input_data/spatial_data/europe_mask_50km.tif")) # Background mask of Europe in a 50 km resolution
load("output_data/data/TBEV_occurrences.RData") # Infection occurrence data
load("output_data/data/I_ricinus_occurrences_cleaned.RData") # The presence information of tick species Ixodes ricinus
r_curr_preds_clim_landuse <- terra::rast("output_data/results/I_ricinus_preds_clim_landuse_1970_2019.tif") # Predictions of occurrence probability of Ixodes ricinus

# Prepare path to data folder
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Create a sequence of dates with monthly steps, that are temporally
# covered by the environmental data and infection data
start_date <- as.Date("2008-01-01") # Define start date
end_date <- as.Date("2019-12-01") # Define end date
date_sequence_month <- seq.Date(from = start_date, to = end_date, by = "month") # create a monthly sequence
date_sequence <- format(date_sequence_month, "%m/%Y") # Extract year and month from dates

# Format the months within the data frame of cleaned occurrences
TBEV_occurrences$month <- sprintf("%02d", TBEV_occurrences$month)

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
  subset_year_month <- subset(TBEV_occurrences, TBEV_occurrences$year == y & TBEV_occurrences$month == m)
  
  if (nrow(subset_year_month) > 0) { # Just continue with preparation process if occurrences are available for the date
    
    print("start with data prep")
    
    
#-------------------------------------------------------------------------------
    
# 1. Remove duplicates in cells ------------------------------------------------
    
    # Select coordinates of the occurrences
    occ_coords <- as.data.frame(subset_year_month[, c("lon", "lat")])
    
    # Extract an ID per cell
    cellnumbers <- terra::extract(europe_mask, occ_coords, cells = TRUE)
    
    # Only keep cells that are not duplicated
    occ_coords <- occ_coords[!duplicated(cellnumbers[,"cell"]),]
    
    
#-------------------------------------------------------------------------------
    
# 2. Generation of background data ---------------------------------------------
    
    print("start with background generation")
    
    # Subset the presence information of associated vector species Ixodes ricinus
    # for the respective year and month
    subset_year_month_species <- subset(I_ricinus_occurrences_cleaned, I_ricinus_occurrences_cleaned$year == y & I_ricinus_occurrences_cleaned$month == m)
    
    # Remove duplicates within a cell for vector species
    # Select coordinates of the occurrences
    occ_coords_species <- as.data.frame(subset_year_month_species[, c("lon", "lat")])
    
    # Extract an ID per cell
    cellnumbers_species <- terra::extract(europe_mask, occ_coords_species, cells = TRUE)
    
    # Only keep cells that are not duplicated
    occ_coords_species <- occ_coords_species[!duplicated(cellnumbers_species[,"cell"]),]
    
    # Vectorise species and disease occurrences
    occ_coords_TBEV_vect <- terra::vect(occ_coords, geom = c("lon", "lat"), crs = "EPSG:4326")
    occ_coords_species_vect <- terra::vect(occ_coords_species, geom = c("lon", "lat"), crs = "EPSG:4326")
    
    # Rasterise vector species and infection occurrence points
    occ_coords_TBEV_r <- terra::rasterize(occ_coords_TBEV_vect, europe_mask, field = 1, background = 0)
    occ_coords_species_r <- terra::rasterize(occ_coords_species_vect, europe_mask, field = 1, background = 0)
    
    # Identify cells with species presences but no infection occurrence
    abs_raster <- occ_coords_species_r == 1 & is.na(occ_coords_TBEV_r)
    
    # Extract cells where the species is present and the disease absent
    abs_raster <- which(values(abs_raster) == 1)
    
    # Convert raster cells back to points
    abs_points <- as.points(abs_raster, values = TRUE, na.rm = TRUE)
    
    # Extract the coordinates of the presence points Spatvector
    abs_coords <- terra::geom(abs_points)
    abs_coords <- data.frame(lon = abs_coords[, "x"], lat = abs_coords[, "y"])
    
    # Add information on presence and background
    occ_coords$occ <- 1
    abs_coords$occ <- 0
    
    
#-------------------------------------------------------------------------------
    
# 3. Spatial thinning of the presence and background data ----------------------
    
    print("spatial thinning")
    
    # Thinning of presences
    # Transform data frame into sf object to use in thin function
    occ_coords_sf <- st_as_sf(occ_coords, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
    
    # Using the thin function with a thinning distance of 100 km (2 cells)
    occ_coords_thinned <- thin(occ_coords_sf, thin_dist = 100000, runs = 1, ncores = 1)
    
    # Merge data frames to only retained thinned presences and background
    occ_coords_thinned <- merge(occ_coords_thinned, occ_coords, by = c("lon", "lat"))
    
    
    # Thinning of background data
    # Transform data frame into sf object to use in thin function
    abs_coords_sf <- st_as_sf(abs_coords_150, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
    
    # Using the thin function with a thinning distance of 100 km (2 cells)
    abs_coords_thinned <- thin(abs_coords_sf, thin_dist = 100000, runs = 1, ncores = 1)
    
    # Merge data frames to only retained thinned presences and background
    abs_coords_thinned <- merge(abs_coords_thinned, abs_coords_150, by = c("lon", "lat"))
    
    
    # Join presence and background data
    TBEV_occ_thinned <- rbind(occ_coords_thinned, abs_coords_thinned)
    
    # Add the year and month as information in columns
    TBEV_occ_thinned$year <- y
    TBEV_occ_thinned$month <- m
    
    
#-------------------------------------------------------------------------------
    
# 4. Join with environmental data  ---------------------------------------------
    
    print("matching env. data")
    
    # Load the environmental data for the specific year and month
    Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
    LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
    Species_data <- r_curr_preds_clim_landuse[[paste0(m,"/",y)]] # Occurrence probabilities of the vector species Ixodes ricinus
    
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
    
    
    
    # Add the thinned presences and background data belonging to the specific year and month to the 
    # prepared results data frame
    TBEV_occ_env <- rbind(TBEV_occ_env, TBEV_occ_env_date)
    
    
    
  } else if (nrow(subset_year_month) == 0) { print("no data available")
  } # End of if-condition
  
} # End of loop over dates
    
    
# Get a summary of presence and background data numbers
table(TBEV_occ_env$occ) # 
print(table(TBEV_occ_env$month[I_ricinus_occ_env$occ == 1]))

# Save the resulting data frame, containing thinned presence and background data,
# joined with the respective environmental data of year and month
save(TBEV_occ_env, file = "output_data/data/TBEV_occ_env.RData")
    


    
    
    


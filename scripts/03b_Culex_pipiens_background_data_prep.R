# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#           03b. Background data preparation - Culex pipiens             #
# ---------------------------------------------------------------------- #


# What is done within this script:

# We generate background data for our vector species Culex pipiens by randomly 
# selecting locations within a specified buffer distance of 100 km from the 
# presence points, aiming for a presence-background ratio of 1:10, excluding 
# cells containing the actual presence locations. This process is conducted 
# separately based on occurrences within the same month of a given year, 
# resulting in temporally matched background data. To avoid spatial 
# autocorrelation, we thin both the monthly presence and background data of the 
# species using a 50 km threshold. Finally, we match the species data - 
# comprising both presence and background data - with the month- and year-specific 
# climate predictors, as well as year-specific land-use predictors.

# Load needed packages
library(terra) # terra_1.7-55
library(sf) # sf_1.0-16
library(sfheaders) # sfheaders_0.4.3
library(tidyverse) # tidyverse_2.0.0


# Load needed objects
source("scripts/00_functions.R") # Get the thin function

# Prepare path to data folder
datapath_env <- file.path("input_data/environmental_data/ISIMIP3a/")

# Create a sequence of dates with monthly steps, that are temporally
# covered by the environmental data
start_date <- as.Date("1970-01-01") # Define start date
end_date <- as.Date("2019-12-01") # Define end date (end of monthly climate data)
date_sequence_month <- seq.Date(from = start_date, to = end_date, by = "month") # Create a monthly sequence
date_sequence <- format(date_sequence_month, "%m/%Y") # Extract year and month from dates

# Read in the background mask of Europe (0.5° resolution)
europe_mask <- terra::rast(paste0("input_data/spatial_data/europe_mask.tif"))

# Load cleaned occurrence data
load("output_data/data/C_pipiens_occurrences_cleaned.RData")

# Format the months within the data frame of cleaned occurrences
C_pipiens_occurrences_cleaned$month <- sprintf("%02d", C_pipiens_occurrences_cleaned$month)

# Create an empty data frame to store the thinned presences
C_pipiens_occ_env <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(C_pipiens_occ_env) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "primary_forest", "primary_openland", "secondary_forest", 
                                 "secondary_openland", "pasture", "rangeland", "cropland", "urban")




# Loop over all dates in the sequence, 
# remove duplicates in 0.5° cells, generate background data,
# and match the presence and background data with the time-specific environmental data.

for (d in date_sequence) { # Start of the loop over all dates
  
  # Extract month and year from date
  m <- substr(d, 1, 2)
  print(m)
  y <- substr(d, 4, 7)
  print(y)
  
  # Subset the cleaned occurrence data frame by each date (month and year)
  subset_year_month <- C_pipiens_occurrences_cleaned[
    C_pipiens_occurrences_cleaned$year  == y & 
      C_pipiens_occurrences_cleaned$month == m, 
  ]
  
  if (nrow(subset_year_month) > 0) { # Just continue with preparation process if occurrences are available for the date
    
    print("start with data prep")
    
    
#-------------------------------------------------------------------------------
    
# 1. Remove duplicates in cells ------------------------------------------------
    
    # Select coordinates of the occurrences
    occ_coords <- as.data.frame(subset_year_month[, c("lon", "lat")])
    
    # Extract an ID per cell
    cellnumbers <- terra::extract(europe_mask, occ_coords, cells = TRUE)
    
    # Remove duplicates based on raster cell (if multiple points fall in the 
    # same raster cell, only one point is kept)
    occ_coords <- occ_coords[!duplicated(cellnumbers[,"cell"]),]
    
    
    
#-------------------------------------------------------------------------------
    
# 2. Generation of background data ---------------------------------------------
    
    # Convert presence points into vectors
    presences <- vect(occ_coords, crs = "+proj=longlat +datum=WGS84")
    
    # Check if presences fall within the European continent and remove those that fall outside
    values_ext <- terra::extract(europe_mask, presences)
    presences_europe <- presences[!is.na(values_ext[, 2])]
    
    # Extract the coordinates of the presence points Spatvector
    occ_coords <- terra::geom(presences_europe)
    occ_coords <- data.frame(lon = occ_coords[, "x"], lat = occ_coords[, "y"])
    
    
    
    if (nrow(occ_coords) > 0) { # Just continue if at least one presence of the respective month and year are within the continent of Europe
      
      # Place a buffer of 100 km around presence locations
      buf_100 <- buffer(presences_europe, width = 100000)
      
      # use mask_buf to rasterize buf_100 (which has been a vector so far; !raster required for later steps)
      buf_100 <- rasterize(buf_100, europe_mask)
      
      # set raster cells outside the buffer to NA
      buf_100 <- terra::mask(europe_mask, buf_100, overwrite = TRUE)
      
      # randomly select background data within the buffer, excluding presence locations (aiming to sample 10x as many background points as presences)
      occ_cells_100 <- terra::extract(buf_100, occ_coords, cells = TRUE)[,"cell"]
      buf_cells_100 <- terra::extract(buf_100, crds(buf_100), cells = TRUE)[,"cell"]
      diff_cells_100 <- setdiff(buf_cells_100, occ_cells_100)
      
      abs_indices_100 <- sample(diff_cells_100, ifelse(length(diff_cells_100) < nrow(occ_coords)*10, length(diff_cells_100), nrow(occ_coords)*10))
      abs_coords_100 <- as.data.frame(xyFromCell(buf_100, abs_indices_100))
      colnames(abs_coords_100) = c("lon", "lat")
      
      # Add information on presence and background
      occ_coords$occ <- 1
      abs_coords_100$occ <- 0
      
      
      
      
#-------------------------------------------------------------------------------
      
# 3. Spatial thinning of the presence and background data ----------------------
      
      print("spatial thinning")
      
      # Thinning of presences
      # Transform data frame into sf object to use in thin function
      occ_coords_sf <- st_as_sf(occ_coords, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
      
      # Using the thin function with a thinning distance of 50 km (aiming for a checkerboard-like pattern)
      occ_coords_thinned <- thin(occ_coords_sf, thin_dist = 50000, runs = 1, ncores = 1)
      
      # Merge data frames to only retain thinned presences
      occ_coords_thinned <- merge(occ_coords_thinned, occ_coords, by = c("lon", "lat"))
      
      
      
      # Thinning of background data
      # Transform data frame into sf object to use in thin function
      abs_coords_sf <- st_as_sf(abs_coords_100, coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")
      
      # Using the thin function with a thinning distance of 50 km (aiming for a checkerboard-like pattern)
      abs_coords_thinned <- thin(abs_coords_sf, thin_dist = 50000, runs = 1, ncores = 1)
      
      # Merge data frames to only retain thinned background data
      abs_coords_thinned <- merge(abs_coords_thinned, abs_coords_100, by = c("lon", "lat"))
      
      
      
      # Join presence and background data
      C_pipiens_occ_thinned <- rbind(occ_coords_thinned, abs_coords_thinned)
      
      # Add the year and month as information in columns
      C_pipiens_occ_thinned$year <- y
      C_pipiens_occ_thinned$month <- m
      
      
      
#-------------------------------------------------------------------------------
      
# 4. Join with environmental data  ---------------------------------------------
      
      print("matching env. data")
      
      # Load the environmental data for the specific year and month
      Climate_data <- terra::rast(paste0(datapath_env, "/Climate/processed_data/Climate_data_",m,"_",y,".tif")) # Climate data
      LandUse_data <- terra::rast(paste0(datapath_env, "/LandUse/processed_data/LandUse_data_",y,".tif")) # Land cover data
      
      # Stack the environmental data
      env_data <- c(Climate_data, LandUse_data)
      
      # Extract the environmental values per occurrence cell
      C_pipiens_occ_env_date <- cbind(C_pipiens_occ_thinned, terra::extract(x = env_data, y = C_pipiens_occ_thinned[,c("lon", "lat")]))
      
      # Drop NA for the environmental variables 
      C_pipiens_occ_env_date <- C_pipiens_occ_env_date %>% drop_na()
      
      # Check for duplicates
      duplicated(C_pipiens_occ_env_date$ID)
      
      # Only retain non-duplicated cells
      C_pipiens_occ_env_date <- C_pipiens_occ_env_date[!duplicated(C_pipiens_occ_env_date$ID),]
      
      
      
      # Add the thinned presences and background data belonging to the specific year and month to the 
      # prepared results data frame
      C_pipiens_occ_env <- rbind(C_pipiens_occ_env, C_pipiens_occ_env_date)
      
      
    } else if (nrow(occ_coords) == 0) { print("no data available for European continent")
    } # End of if-condition
    
    
  } else if (nrow(subset_year_month) == 0) { print("no data available for month-year combination")
  } # End of if-condition
  
} # End of loop over dates



# Get a summary of presence and background data numbers
print(table(C_pipiens_occ_env$occ))
print(table(C_pipiens_occ_env$month[C_pipiens_occ_env$occ == 1]))

# Save the resulting data frame, containing thinned presence and background data,
# joined with the respective environmental data of year and month
save(C_pipiens_occ_env, file = "output_data/data/C_pipiens_occ_env.RData")




#-------------------------------------------------------------------------------

# 5. Occurrence mapping --------------------------------------------------------

# Map the thinned presences and background data
png("output_data/plots/presence_background/C_pipiens_presence_pseudoabsence.png", width = 2000, height = 2000, res = 300)


maps::map('world',xlim=c(-31,40), ylim=c(34,72), 
          col = "gray97",
          fill = TRUE,
          border = "gray30")

maps::map.axes(cex.axis = 0.75)


points(C_pipiens_occ_env$lon[C_pipiens_occ_env$occ == 0], C_pipiens_occ_env$lat[C_pipiens_occ_env$occ == 0], col='steelblue4',  pch=19, cex = 0.5)
points(C_pipiens_occ_env$lon[C_pipiens_occ_env$occ == 1], C_pipiens_occ_env$lat[C_pipiens_occ_env$occ == 1], col='goldenrod',  pch=19, cex = 0.5)

legend(title = "Culex pipiens:", x = -28, y = 50, legend = c("Pseudoabsence", "Presence"), col = c("steelblue4", "goldenrod"), pch = 19, pt.cex = 1, bty = "n")

dev.off()

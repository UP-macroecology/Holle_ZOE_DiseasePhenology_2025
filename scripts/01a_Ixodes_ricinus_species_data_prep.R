# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                    01a. Species data preparation                       #
# ---------------------------------------------------------------------- #


# Load needed packages
library(rgbif)
library(CoordinateCleaner)
library(dplyr)
library(tibble)
library(maps)
library(lubridate)


# Load needed objects
# Load downloaded occurrences that were available from VectorMap
I_ricinus_vectormap <- read.csv("input_data/raw_species_data/Ixodes_ricinus_occurrences_VectorMap_020824.csv") 


#-------------------------------------------------------------------------------

# 1. Download GBIF data --------------------------------------------------------

# Check for synonyms
name_suggest(q = "Ixodes ricinus", rank = "species")

# Check how many records with coordinate information are available for the extent of Europe
# within more recent years
occ_count(scientificName = "Ixodes ricinus", hasCoordinate = TRUE,
          decimalLongitude = "-31,40", decimalLatitude = "34,72",
          year = "1970,2019") # 4447

# Download these occurrences
I_ricinus_gbif_list <- occ_search(scientificName = "Ixodes ricinus", hasCoordinate = TRUE,
                                  decimalLongitude = "-31,40", decimalLatitude = "34,72",
                                  year = "1970,2019", limit = 10000)

# Extract data frame containing occurrence records
I_ricinus_gbif <- I_ricinus_gbif_list$data


# Subset data frame to only contain relevant columns
I_ricinus_gbif <- I_ricinus_gbif[, c("species", "decimalLatitude","decimalLongitude",
                                     "country", "year", "month", "institutionCode", 
                                     "coordinateUncertaintyInMeters")]

# Change the column names
colnames(I_ricinus_gbif) <- c("species", "lat", "lon", "country", "year", "month", "datasource", "coordinate_uncertainty")

# Add a column name that indicates the used database
I_ricinus_gbif$database <- "GBIF"



#-------------------------------------------------------------------------------

# 2. Downloaded VectorMap data -------------------------------------------------

# Subset the data frame to only contain relevant columns
I_ricinus_vectormap <- I_ricinus_vectormap[, c("ScientificName", "DecimalLatitude", "DecimalLongitude",
                                               "Country", "EarliestYearCollected", "EarliestDateCollected", "InstitutionCode", 
                                               "CoordinateUncertaintyInMeters")]

# Extract the month of collection from Date format
earliestdatecollected <- I_ricinus_vectormap$EarliestDateCollected
I_ricinus_vectormap$EarliestDateCollected <- sapply(earliestdatecollected, function(x) month(mdy_hms(x)))

# Change the column names
colnames(I_ricinus_vectormap) <- c("species", "lat", "lon", "country", "year", "month", "datasource", "coordinate_uncertainty")

# Only retain occurrences on the European continent
I_ricinus_vectormap <- I_ricinus_vectormap[I_ricinus_vectormap$lon >= -31 & I_ricinus_vectormap$lon <= 40 & 
                                             I_ricinus_vectormap$lat >= 34 & I_ricinus_vectormap$lat <= 72, ]


# Add a column name that indicates the used database
I_ricinus_vectormap$database <- "VectorMap"


#-------------------------------------------------------------------------------

# 3. Occurrence cleaning -------------------------------------------------------

  
# Bind the occurrences from the two sources in one data frame
I_ricinus_occurrences <- rbind(I_ricinus_gbif, I_ricinus_vectormap)

# Clean coordinates for the monthly SDMs (keep occurrences with the same coordinates
# if recorded in different months)
I_ricinus_occurrences_cleaned <- I_ricinus_occurrences %>%
  mutate_at(vars(lon, lat), round, 4) %>%  
  dplyr::filter(!(is.na(lat) | is.na(lon)), # only records with coordinates
                !(lat == lon | lat == 0 | lon == 0), # coordinates should not be equal or zero
                !is.na(month), # only records with specified month
                !(is.na(coordinate_uncertainty) | coordinate_uncertainty > 25000), # coordinate precision < 25000m  to still fall within cell (50km resolution)
                !(year < 1970 | year > 2019)) %>% # recent years
  distinct(lon, lat, month, year, .keep_all = TRUE) %>% # remove entries with duplicate coordinate records within the same month of the same year
  clean_coordinates(lon = "lon", lat = "lat", species = "species", countries = "country", 
                    tests = c("centroids", "outliers", "institutions"))

# Extract the cleaned occurrence locations
I_ricinus_occurrences_cleaned <- I_ricinus_occurrences_cleaned %>%  
  dplyr::filter(.summary == TRUE) %>% # remove occurrences that were flagged by coordinateCleaner
  rowid_to_column(var = "occ_id") %>% # create unique identifier for each occurrence
  dplyr::select(occ_id, species, lon, lat, year, month, country, datasource, coordinate_uncertainty, database) # select only relevant columns



#-------------------------------------------------------------------------------

# 4. Map occurrence points -----------------------------------------------------

# Plot a map with the extent of Europe and plot all remaining points after cleaning
maps::map('world',xlim=c(-31,40), ylim=c(34,72))
points(I_ricinus_occurrences_cleaned$lon, I_ricinus_occurrences_cleaned$lat, col='goldenrod',  pch=19, cex = 0.5)



#-------------------------------------------------------------------------------

# 5. Save resulting data -------------------------------------------------------
save(I_ricinus_occurrences_cleaned, file = "output_data/data/I_ricinus_occurrences_cleaned.RData")
  


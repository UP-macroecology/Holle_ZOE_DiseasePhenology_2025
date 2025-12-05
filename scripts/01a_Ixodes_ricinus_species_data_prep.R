# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#            01a. Species data preparation - Ixodes ricinus              #
# ---------------------------------------------------------------------- #


# What is done within this script:

# We process raw Ixodes ricinus occurrence data (GBIF, VectorMap) to retain only
# the relevant columns; clean the coordinates by removing records with erroneous 
# timestamps and coordinates; and remove duplicate occurrences reported within
# the same calendar month of the same year to achieve a target temporal resolution 
# of one month.



# Load needed packages
library(CoordinateCleaner) # CoordinateCleaner_3.0.1
library(tidyverse) # tidyverse_2.0.0
library(tibble) # tibble_3.2.1
library(lubridate) # lubridate_1.9.3
library(sf) # sf_1.0-16
library(maps) # maps_3.4.1


# Load needed objects
# Load downloaded occurrences that were available from VectorMap and GBIF
I_ricinus_vectormap <- read.csv("input_data/raw_species_data/Ixodes_ricinus_occurrences_VectorMap_020824.csv")
I_ricinus_gbif <- read.delim("input_data/raw_species_data/Ixodes_ricinus_occurrences_GBIF_280125.csv") 



#-------------------------------------------------------------------------------

# 1. Process raw GBIF data -----------------------------------------------------

# GBIF data was downloaded with the following filters: Geometry: POLYGON((-31 34,40 34,40 72,-31 72,-31 34)),
# Has Coordinate: TRUE, Scientific name: Ixodes ricinus, Year: between start of 1970 and end of 2019
# GBIF.org (28 January 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.c2krqj
# DOI: 10.15468/dl.c2krqj
# 4483 occurrences included in download

# {
#"and" : [
#  "Geometry POLYGON((-31 34,40 34,40 72,-31 72,-31 34))",
#  "HasCoordinate is true",
#  "HasGeospatialIssue is false",
#  "TaxonKey is Ixodes ricinus Linnaeus, 1758",
#  "Year 1970-2019"
#]
#}

# Subset data frame to only contain relevant columns
I_ricinus_gbif <- I_ricinus_gbif[, c("species", "decimalLatitude","decimalLongitude",
                                     "countryCode", "year", "month", "institutionCode", 
                                     "coordinateUncertaintyInMeters")]

# Change the column names
colnames(I_ricinus_gbif) <- c("species", "lat", "lon", "country", "year", "month", "datasource", "coordinate_uncertainty")

# Add a column name that indicates the used database
I_ricinus_gbif$database <- "GBIF"



#-------------------------------------------------------------------------------

# 2. Process raw VectorMap data ------------------------------------------------

# Remove rows that have no entry/NAs in their date columns
I_ricinus_vectormap <- I_ricinus_vectormap %>%
  filter(!is.na(EarliestDateCollected) & EarliestDateCollected != "",
         !is.na(LatestDateCollected) & LatestDateCollected != "")

# Only keep occurrences where the columns "EarliestDateCollected" and
# "LatestDateCollected" coincide in their year and month
I_ricinus_vectormap <- I_ricinus_vectormap %>% 
  mutate(EarliestDateCollected = mdy_hms(EarliestDateCollected),
         LatestDateCollected = mdy_hms(LatestDateCollected),
         EarliestYearMonth = format(EarliestDateCollected, "%Y-%m"),
         LatestYearMonth = format(LatestDateCollected, "%Y-%m")
  ) %>%
  filter(EarliestYearMonth == LatestYearMonth) %>%
  select(-EarliestYearMonth, -LatestYearMonth)

# Subset the data frame to only contain relevant columns
I_ricinus_vectormap <- I_ricinus_vectormap[, c("ScientificName", "DecimalLatitude", "DecimalLongitude",
                                               "Country", "EarliestYearCollected", "EarliestDateCollected", "InstitutionCode", 
                                               "CoordinateUncertaintyInMeters")]

# Extract the month of collection from Date format
earliestdatecollected <- I_ricinus_vectormap$EarliestDateCollected
I_ricinus_vectormap$EarliestDateCollected <- sapply(earliestdatecollected, function(x) month(ymd(x)))

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

# Clean coordinates in preparation for the spatiotemporal (monthly) SDMs
# (keep occurrences with the same coordinates if recorded in different months)
I_ricinus_occurrences_cleaned <- I_ricinus_occurrences %>%
  mutate_at(vars(lon, lat), round, 4) %>%  
  dplyr::filter(!(is.na(lat) | is.na(lon)), # Only records with coordinates
                !(lat == lon | lat == 0 | lon == 0), # Coordinates should not be equal or zero
                !is.na(month), # Only records with specified month
                !(is.na(coordinate_uncertainty) | coordinate_uncertainty > 25000), # Coordinate precision < 25000m  to still fall within cell (0.5°)
                !(year < 1970 | year > 2019)) %>% # Recent years
  distinct(lon, lat, month, year, .keep_all = TRUE) %>% # Remove entries with duplicate coordinate records within the same month of the same year
  clean_coordinates(lon = "lon", lat = "lat", species = "species", countries = "country", 
                    tests = c("centroids", "outliers", "institutions"))

# Extract the cleaned occurrence locations
I_ricinus_occurrences_cleaned <- I_ricinus_occurrences_cleaned %>%  
  dplyr::filter(.summary == TRUE) %>% # Remove occurrences that were flagged by coordinateCleaner
  rowid_to_column(var = "occ_id") %>% # Create unique identifier for each occurrence
  dplyr::select(occ_id, species, lon, lat, year, month, country, datasource, coordinate_uncertainty, database) # Select only relevant columns


                                             

#-------------------------------------------------------------------------------

# 4. Map occurrence points -----------------------------------------------------

# Plot a map with the extent of Europe and plot all remaining points after cleaning
# for visual inspection
maps::map('world',xlim=c(-31,40), ylim=c(34,72), 
          col = "gray97",
          fill = TRUE,
          border = "gray30")

maps::map.axes()

points(I_ricinus_occurrences_cleaned$lon, I_ricinus_occurrences_cleaned$lat, col='goldenrod',  pch=19, cex = 0.5)




#-------------------------------------------------------------------------------

# 5. Save resulting data -------------------------------------------------------
save(I_ricinus_occurrences_cleaned, file = "output_data/data/I_ricinus_occurrences_cleaned.RData")
  


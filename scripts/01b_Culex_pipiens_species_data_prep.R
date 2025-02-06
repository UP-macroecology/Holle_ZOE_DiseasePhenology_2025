# ZOE project 
# Disease phenology analysis of Culex pipiens in Europe (primary transmitter of WNV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                     01b. Species data preparation                      #
# ---------------------------------------------------------------------- #


# Load needed packages
library(CoordinateCleaner)
library(tibble)
library(maps)
library(lubridate)
library(tidyverse)


# Load needed objects
# Load downloaded occurrences that were available from VectorMap and GBIF
C_pipiens_vectormap <- read.csv("input_data/raw_species_data/Culex_pipiens_occurrences_VectorMap_200924.csv") 
C_pipiens_gbif <- read.delim("input_data/raw_species_data/Culex_pipiens_occurrences_GBIF_280125.csv") 


#-------------------------------------------------------------------------------

# 1. Download GBIF data --------------------------------------------------------

# GBIF data was downloaded with the following filters: Geometry: POLYGON((-31 34,40 34,40 72,-31 72,-31 34)),
# Has Coordinate: TRUE, Scientific name: Culex pipiens, Year: between start of 1970 and end of 2019
# GBIF.org (28 January 2025) GBIF Occurrence Download https://doi.org/10.15468/dl.2q36w8
# DOI: 10.15468/dl.2q36w8
# 7257 occurrences included in download

# {
#"and" : [
#  "Geometry POLYGON((-31 34,40 34,40 72,-31 72,-31 34))",
#  "HasCoordinate is true",
#  "HasGeospatialIssue is false",
#  "TaxonKey is Culex pipiens Linnaeus, 1758",
#  "Year 1970-2019"
#]
#}


# Subset data frame to only contain relevant columns
C_pipiens_gbif <- C_pipiens_gbif[, c("species", "decimalLatitude","decimalLongitude",
                                     "countryCode", "year", "month", "institutionCode", 
                                     "coordinateUncertaintyInMeters")]

# Change the column names
colnames(C_pipiens_gbif) <- c("species", "lat", "lon", "country", "year", "month", "datasource", "coordinate_uncertainty")

# Add a column name that indicates the used database
C_pipiens_gbif$database <- "GBIF"


#-------------------------------------------------------------------------------

# 2. Downloaded VectorMap data -------------------------------------------------

# Only keep occurrences where the columns "EarliestDateCollected" and
# "LatestDateCollected" coincide in their year and month
C_pipiens_vectormap <- C_pipiens_vectormap %>% 
  mutate(EarliestDateCollected = ymd_hm(EarliestDateCollected),
         LatestDateCollected = ymd_hm(LatestDateCollected),
         EarliestYearMonth = format(EarliestDateCollected, "%Y-%m"),
         LatestYearMonth = format(LatestDateCollected, "%Y-%m")
  ) %>%
  filter(EarliestYearMonth == LatestYearMonth) %>%
  select(-EarliestYearMonth, -LatestYearMonth)

# Subset the data frame to only contain relevant columns
C_pipiens_vectormap <- C_pipiens_vectormap[, c("ScientificName", "DecimalLatitude", "DecimalLongitude",
                                               "Country", "EarliestYearCollected", "EarliestDateCollected", "InstitutionCode", 
                                               "CoordinateUncertaintyInMeters")]

# Extract the month of collection from Date format
earliestdatecollected <- C_pipiens_vectormap$EarliestDateCollected
C_pipiens_vectormap$EarliestDateCollected <- sapply(earliestdatecollected, function(x) month(mdy_hms(x)))

# Change the column names
colnames(C_pipiens_vectormap) <- c("species", "lat", "lon", "country", "year", "month", "datasource", "coordinate_uncertainty")

# Only retain occurrences on the European continent
C_pipiens_vectormap <- C_pipiens_vectormap[C_pipiens_vectormap$lon >= -31 & C_pipiens_vectormap$lon <= 40 & 
                                            C_pipiens_vectormap$lat >= 34 & C_pipiens_vectormap$lat <= 72, ]


# Add a column name that indicates the used database
C_pipiens_vectormap$database <- "VectorMap"


#-------------------------------------------------------------------------------

# 3. Occurrence cleaning -------------------------------------------------------


# Bind the occurrences from the two sources in one data frame
C_pipiens_occurrences <- rbind(C_pipiens_gbif, C_pipiens_vectormap)

# Clean coordinates for the monthly SDMs (keep occurrences with the same coordinates
# if recorded in different months)
C_pipiens_occurrences_cleaned <- C_pipiens_occurrences %>%
  mutate_at(vars(lon, lat), round, 4) %>%  
  dplyr::filter(!(is.na(lat) | is.na(lon)), # only records with coordinates
                !(lat == lon | lat == 0 | lon == 0), # coordinates should not be equal or zero
                !is.na(month), # only records with specified month
                !(is.na(coordinate_uncertainty) | coordinate_uncertainty > 25000), # coordinate precision < 25000m  to still fall within cell (50km2)
                !(year < 1970 | year > 2019)) %>% # recent years
  distinct(lon, lat, month, year, .keep_all = TRUE) %>% # remove entries with duplicate coordinate records within the same month of the same year
  clean_coordinates(lon = "lon", lat = "lat", species = "species", countries = "country", 
                    tests = c("centroids", "outliers", "institutions"))

# Extract the cleaned occurrence locations
C_pipiens_occurrences_cleaned <- C_pipiens_occurrences_cleaned %>%  
  dplyr::filter(.summary == TRUE) %>% # remove occurrences that were flagged by coordinateCleaner
  rowid_to_column(var = "occ_id") %>% # create unique identifier for each occurrence
  dplyr::select(occ_id, species, lon, lat, year, month, country, datasource, coordinate_uncertainty, database) # select only relevant columns



#-------------------------------------------------------------------------------

# 4. Map occurrence points -----------------------------------------------------

# Plot a map with the extent of Europe and plot all remaining points after cleaning
maps::map('world',xlim=c(-31,40), ylim=c(34,72))
points(C_pipiens_occurrences_cleaned$lon, C_pipiens_occurrences_cleaned$lat, col='goldenrod',  pch=19, cex = 0.5)



#-------------------------------------------------------------------------------

# 5. Save resulting data -------------------------------------------------------
save(C_pipiens_occurrences_cleaned, file = "output_data/data/C_pipiens_occurrences_cleaned.RData")
  
  



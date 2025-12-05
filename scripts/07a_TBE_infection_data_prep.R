# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                   07a. Infection data preparation - TBE                #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We process locally acquired, confirmed TBE cases in Europe (provided by TESSy/ECDC)
# to generate spatially explicit infection data. This is done by rasterising
# the European NUTS3 municipalities to the target spatial resolution of 0.5° and
# extracting the central coordinates of municipalities with observed infections. 
# For NUTS3 municipalities that consist of only one cell after rasterisation, 
# these coordinates are used as the final location for the infection data point.
# In cases where NUTS3 municipalities are too small to be represented by a 0.5°
# cell in the rasterisation process but had reported TBE infections, their 
# central coordinates are likewise considered as the infection data point.
# Infection records stemming from large NUTS3 municipalities spanning more than 
# one cell are excluded to minimise the spatial uncertainty of infection locations.
# Finally, we examine which countries provided information on TBE infections
# and identify those reporting at the NUTS3 level. This is important for generating
# the background dataset in script 08a.


# Load needed packages
library(giscoR) # giscoR_0.6.0
library(terra) # terra_1.7-55
library(ggplot2) # ggplot2_4.0.0
library(tidyverse) # tidyverse_2.0.0
library(sf) # sf_1.0-16
library(countrycode) # countrycode_1.6.0

# Load needed data
europe_mask <- terra::rast("input_data/spatial_data/europe_mask.tif") # 0.5° raster template of Europe
TBE_infection_data <- read.csv("input_data/raw_infection_data/TBE.csv") # TBE infection data provided by ECDC/TESSy



#-------------------------------------------------------------------------------

# 1. Extracting the uncertainty of ECDC disease data ---------------------------

# ECDC human case infection data is provided at the NUTS3 level,
# Retain a map showing the European municipalities on NUTS 3 level (year 2021)
# (As we use data until 2019; UK was still reporting surveillance data to ECDC)
# "4326": WGS84
nuts_3 <- gisco_get_nuts(
  year = "2021",
  epsg = "4326",
  cache = TRUE,
  update_cache = FALSE,
  cache_dir = NULL,
  verbose = FALSE,
  resolution = "10",
  spatialtype = "RG",
  country = NULL,
  nuts_id = NULL,
  nuts_level = "3"
)

# Rasterise the NUTS3 multipolygon and mask the values that do not belong to the European continent
nuts_3_raster <- terra::rasterize(nuts_3, europe_mask, field = "NUTS_ID", touches = TRUE)
nuts_3_raster_mask <- terra::mask(nuts_3_raster, europe_mask)

# Store the raster map of nuts 3 municipalities for later usage (during background data generation)
writeRaster(nuts_3_raster_mask, "input_data/spatial_data/nuts_3_raster_mask.tif", overwrite = TRUE)

# Count the number of raster cells for each nuts3 municipality
cell_counts_nuts_3 <- terra::freq(nuts_3_raster_mask) %>%
  as.data.frame() %>%
  rename(NUTS_ID = value, num_cells = count)

# Join the cell counts back to the nuts_3 multipolygon data
nuts_3 <- nuts_3 %>%
  left_join(cell_counts_nuts_3, by = "NUTS_ID")

# Rasterise the cell count into the 0.5° raster and mask values based on raster template
nuts_3_cell_count <- terra::rasterize(nuts_3, europe_mask, field = "num_cells", touches = TRUE)
nuts_3_cell_count <- terra::mask(nuts_3_cell_count, europe_mask)

# Convert raster to data frame for ggplot2 visualization
nuts_3_cell_count_df <- as.data.frame(nuts_3_cell_count, xy = TRUE, na.rm = TRUE)

# Plot the raster showing the number of cells within a NUTS3 municipality (uncertainty)
ggplot(nuts_3_cell_count_df, aes(x = x, y = y, fill = num_cells)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cell count per NUTS3\nmunicipality (0.5° resolution)", option = "viridis") +
  theme_minimal() +
  labs(
    title = "Uncertainty of ECDC infection data on NUTS3 level",
    x = "Longitude", y = "Latitude"
  ) + theme(
    plot.title = element_text(size = 18, face = "bold")
  )

# Save the plot
ggsave("output_data/plots/maps/NUTS3_uncertainty.png", width = 8, height = 5)




#-------------------------------------------------------------------------------

# 2. Process TBE infection data ------------------------------------------------

# Provided by TESSy and ECDC (human case infection data aggregated by NUTS3 level 
# (place of infection, year of infection, month of infection)

# There are two columns referring to the place of infection (PlaceOfInfection,
# PlaceOfInfectionEVD), create a column that brings the info of the two columns
# together
TBE_infection_data <- TBE_infection_data %>%
  mutate(NUTS_ID = case_when(
    (PlaceofInfection %in% c("NULL", "UNK") & PlaceOfInfectionEVD %in% c("NULL", "UNK")) ~ "NULL", # Both columns are "NULL" or "UNK"
    (PlaceofInfection %in% c("NULL", "UNK") & !PlaceOfInfectionEVD %in% c("NULL", "UNK")) ~ PlaceOfInfectionEVD, # PlaceOfInfection is "NULL" or "UNK", use PlaceOfInfectionEVD
    (!PlaceofInfection %in% c("NULL", "UNK") & PlaceOfInfectionEVD %in% c("NULL", "UNK")) ~ PlaceofInfection, # PlaceOfInfectionEVD is "NULL" or "UNK", use PlaceOfInfection
    (nchar(PlaceofInfection) == 5 & nchar(PlaceOfInfectionEVD) != 5) ~ PlaceofInfection, # If both columns have NUTS entries, consider the one with 5 strings as these refer to NUTS3 categorisations
    (nchar(PlaceofInfection) != 5 & nchar(PlaceOfInfectionEVD) == 5) ~ PlaceOfInfectionEVD,
    (nchar(PlaceofInfection) == 5 & nchar(PlaceOfInfectionEVD) == 5) ~ PlaceOfInfectionEVD,
    TRUE ~ NA_character_
  ))


# Only keep rows of confirmed cases that were not imported and were reported until
# the year 2020 (as these are not covered by environmental data)
TBE_infection_data <- subset(TBE_infection_data, 
                             Classification == "CONF" & 
                               DateUsedForStatisticsYear <= 2019 & 
                               Imported == "N")


# Join the infection data with NUTS3 geographic information
nuts_3_TBE <- nuts_3 %>%
  left_join(TBE_infection_data, by = "NUTS_ID")

# Remove rows of municipalities that do not have any observed infection data 
nuts_3_TBE <- nuts_3_TBE[!is.na(nuts_3_TBE$NumberOfCases), ]

# Extract the centroid information of each NUTS3 municipality and append the 
# info to the data frame
nuts_3_TBE$centroid <- st_centroid(nuts_3_TBE$geometry)
nuts_3_TBE$x <- st_coordinates(nuts_3_TBE$centroid)[, 1]
nuts_3_TBE$y <- st_coordinates(nuts_3_TBE$centroid)[, 2]

# Remove infection entries within the same NUTS3 municipality if occurring in the same 
# month of a certain year
nuts_3_TBE <- nuts_3_TBE %>%
  distinct(NUTS_ID, DateUsedForStatisticsMonth, DateUsedForStatisticsYear, .keep_all = TRUE)

# For each infection entry, check the number of 0.5° cells within the reporting
# municipality. If the NUTS 3 municipality consists of more than one 0.5° cell, 
# randomly select one cell within the municipality, extract their central 
# x and y coordinate and fill this information into the data frame
for (i in 1:nrow(nuts_3_TBE)) { # Start of the loop over all rows
  
  num_cells <- nuts_3_TBE$num_cells[i] # Extract information of cell number
  
  if (!is.na(num_cells) && num_cells > 1) {
    
    # Extract the raster cells of the respective municipality
    nuts_ID <- nuts_3_TBE$NUTS_ID[i]
    r_nuts_ID <- ifel(nuts_3_raster_mask == nuts_ID, 1, NA)
    
    # Count the number of cells and sample one random cell within the municipality
    count_cells <- sum(values(r_nuts_ID) == 1, na.rm = TRUE)
    random_cell_index <- sample(1:count_cells, 1)
    cell_number <- which(values(r_nuts_ID) == 1)[random_cell_index]
    
    # Extract the central coordinates of that cell and add them to the infection data frame
    sampled_coordinates <- xyFromCell(r_nuts_ID, cell_number)
    nuts_3_TBE$x[i] <- sampled_coordinates[1]
    nuts_3_TBE$y[i] <- sampled_coordinates[2]
    
  }
} # Close the loop over all rows

# As some NUTS3 municipalities are not depicted as cells due to their small size,
# we assign a column indicating the new NUTS3 municipality code they are belonging to now
# by finding the NUTS3 IDs for each entry infection point
# Remove the multipolygon column
nuts_3_TBE <- st_drop_geometry(nuts_3_TBE)
infection_points <- terra::vect(nuts_3_TBE, geom = c("x", "y"), crs = crs(nuts_3_raster_mask))
NUTS3_values <- terra::extract(nuts_3_raster_mask, infection_points)
nuts_3_TBE$new_municipality <- NUTS3_values[, 2] 

# Save the data frame stating all infection occurrences and their respective
# location (for later usage in absence generation)
save(nuts_3_TBE, file = "output_data/data/nuts_3_TBE.RData")

# Remove entries that stem from municipalities that consist of more than one cell
# as this increases the uncertainty of the reported location (keep entries with
# NA values as these belong to NUTS3 municipalities that were too small to be
# rasterised with a 0.5° resolution)
nuts_3_TBE_filtered <- nuts_3_TBE[nuts_3_TBE$num_cells <= 1 | is.na(nuts_3_TBE$num_cells), ]

# Select only relevant columns of the data frame
TBE_occurrences <- nuts_3_TBE_filtered %>%
  rowid_to_column(var = "occ_id") %>% # create unique identifier for each occurrence
  dplyr::select(occ_id, x, y, DateUsedForStatisticsYear, DateUsedForStatisticsMonth, NUTS_ID)

# Remove the multipolygon column
TBE_occurrences <- st_drop_geometry(TBE_occurrences)

# Change the name of the columns
colnames(TBE_occurrences) <- c("occ_id", "lon", "lat", "year", "month", "NUTS_ID")

# Add a column indicating the presence of infection
TBE_occurrences$occ <- 1

# Convert the occurrence data frame to a spatial object
TBE_occurrences_sp <- st_as_sf(TBE_occurrences, coords = c("lon", "lat"), crs = st_crs(europe_mask))

# Extract raster values at the coordinate locations
occurrences_values <- terra::extract(europe_mask, TBE_occurrences_sp)

# Only keep the occurrences where the values is 1 (meaning that it lays on the European continent)
TBE_occurrences_cleaned <- TBE_occurrences[!is.na(occurrences_values[,"layer"]) & occurrences_values[,"layer"] == 1, ]



#-------------------------------------------------------------------------------

# 3. Visualisation of processed TBE infection data -----------------------------

# Plot the extracted infection occurrences
ggplot(nuts_3_cell_count_df, aes(x = x, y = y, fill = num_cells)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cell count\nper NUTS3 municipality", option = "viridis") +
  geom_point(data = TBE_occurrences_cleaned, aes(x = lon, y = lat, color = as.factor(occ)), size = 2, inherit.aes = FALSE) +
  scale_color_manual(
    name = "TBE infection",
    values = c("1" = "red")
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of ECDC infection data on NUTS3 level & TBE infection data",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold") 
  )




#-------------------------------------------------------------------------------

# 4. Save processed TBE infection data -----------------------------------------

# Save the data frame with occurrence points
save(TBE_occurrences_cleaned, file = "output_data/data/TBE_occurrences_cleaned.RData")



#-------------------------------------------------------------------------------

# 5. Check countries that provided infection data ------------------------------

# Create a vector containing all reporting countries
reporting_countries <- unique(TBE_infection_data$ReportingCountry)

# Get the country names based on ISO 2-Letter Code
reporting_countries_full <- countrycode(reporting_countries, origin = "iso2c", destination = "country.name")
print(reporting_countries_full)



#-------------------------------------------------------------------------------

# 6. Check NUTS3-level reportings ----------------------------------------------

# Add a column for code length and inferred NUTS level to original data frame
# with reported infections
TBE_infection_data_NUTS_check <- TBE_infection_data %>%
  mutate(
    nuts_length = nchar(NUTS_ID),
    country_code = substr(NUTS_ID, 1, 2),
    level = case_when(
      nuts_length == 2 ~ "Country",
      nuts_length == 3 ~ "NUTS1",
      nuts_length == 4 ~ "NUTS2",
      nuts_length == 5 ~ "NUTS3",
      TRUE ~ "Other"
    )
  )

# Summarise which countries did not report at NUTS3
countries_not_nuts3 <- TBE_infection_data_NUTS_check %>%
  group_by(country_code) %>%
  summarise(
    reported_levels = paste(unique(level), collapse = ", "),
    reported_nuts3 = "NUTS3" %in% level
  ) %>%
  filter(!reported_nuts3)

print(countries_not_nuts3)


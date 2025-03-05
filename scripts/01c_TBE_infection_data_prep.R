# ZOE project 
# Disease phenology analysis of TBEV in Europe 

# ---------------------------------------------------------------------- #
#                 01c. TBEV infection data preparation                   #
# ---------------------------------------------------------------------- #



# Load needed packages
library(giscoR)
library(terra)
library(ggplot2)
library(tidyverse)
library(sf)

# Load needed data
europe_mask_50km <- terra::rast("input_data/spatial_data/europe_mask_50km.tif") # 50 km raster template of Europe
TBEV_infection_data <- read.csv("input_data/raw_infection_data/TBE.csv") # TBE infection data provided by ECDC/TESSy



#-------------------------------------------------------------------------------

# 1. Extracting the uncertainty of ECDC disease data ---------------------------

# ECDC human case infection data is provided at the NUTS3 level,
# Retain a map showing the European municipalities on NUTS 3 level (year 2021)
# (As we use data until 2019; UK was still reporting surveillance data to ECDC)
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
nuts_3_raster <- terra::rasterize(nuts_3, europe_mask_50km, field = "NUTS_ID", touches = TRUE)
nuts_3_raster_mask <- terra::mask(nuts_3_raster, europe_mask_50km)

# Store the raster map of nuts 3 municipalities for later usage during background data generation
writeRaster(nuts_3_raster_mask, "input_data/spatial_data/nuts_3_raster_mask.tif", overwrite = TRUE)

# Count the number of raster cells for each polygon
cell_counts_nuts_3 <- terra::freq(nuts_3_raster_mask) %>%
  as.data.frame() %>%
  rename(NUTS_ID = value, num_cells = count)

# Join the cell counts back to the nuts_3 multipolygon data
nuts_3 <- nuts_3 %>%
  left_join(cell_counts_nuts_3, by = "NUTS_ID")

# Rasterise the cell count into the 50km raster and mask values based on raster template
nuts_3_cell_count <- terra::rasterize(nuts_3, europe_mask_50km, field = "num_cells", touches = TRUE)
nuts_3_cell_count <- terra::mask(nuts_3_cell_count, europe_mask_50km)

# Convert raster to data frame for ggplot2 visualization
nuts_3_cell_count_df <- as.data.frame(nuts_3_cell_count, xy = TRUE, na.rm = TRUE)

# Plot the raster showing the number of cells within a NUTS3 municipality (uncertainty)
ggplot(nuts_3_cell_count_df, aes(x = x, y = y, fill = num_cells)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cell count per NUTS3\nmunicipality (50km resolution)", option = "viridis") +
  theme_minimal() +
  labs(
    title = "Uncertainty of ECDC infection data on NUTS3 level",
    x = "Longitude", y = "Latitude"
  ) + theme(
    plot.title = element_text(size = 18, face = "bold")
  )


ggsave("output_data/plots/maps/NUTS3_uncertainty.png", width = 8, height = 5)




#-------------------------------------------------------------------------------

# 2. TBEV infection data -------------------------------------------------------

# Provided by TESSy (human case infection data aggregated by NUTS3 level (place of infection),
# year of infection, month of infection)

# There are two columns referring to the place of infection (PlaceOfInfection,
# PlaceOfInfectionEVD), create a column that brings the info of the two columns
# together
TBEV_infection_data <- TBEV_infection_data %>%
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
TBEV_infection_data <- subset(TBEV_infection_data, 
                              Classification == "CONF" & 
                                DateUsedForStatisticsYear <= 2019 & 
                                Imported == "N")


# Join the infection data with NUTS3 geographic information
nuts_3_TBEV <- nuts_3 %>%
  left_join(TBEV_infection_data, by = "NUTS_ID")

# Remove rows of municipalities that do not have any observed infection data 
nuts_3_TBEV <- nuts_3_TBEV[!is.na(nuts_3_TBEV$NumberOfCases), ]

# Extract the centroid information of each NUTS3 municipality and append the 
# info to the data frame
nuts_3_TBEV$centroid <- st_centroid(nuts_3_TBEV$geometry)
nuts_3_TBEV$x <- st_coordinates(nuts_3_TBEV$centroid)[, 1]
nuts_3_TBEV$y <- st_coordinates(nuts_3_TBEV$centroid)[, 2]

# Remove infection entries within the same cell if occurring in the same 
# month of a certain year
nuts_3_TBEV <- nuts_3_TBEV %>%
  distinct(NUTS_ID, DateUsedForStatisticsMonth, DateUsedForStatisticsYear, .keep_all = TRUE)

# For each infection entry, check the number of 50 km cells within the reporting
# municipality. If the NUTS 3 municipality consists of more than one 50 km cell, 
# randomly select one cell within the municipality, extract their central 
# x and y coordinate and fill this information into the data frame
for (i in 1:nrow(nuts_3_TBEV)) { # Start of the loop over all rows
  
  num_cells <- nuts_3_TBEV$num_cells[i] # Extract information of cell number
  
  if (!is.na(num_cells) && num_cells > 1) {
    
    nuts_ID <- nuts_3_TBEV$NUTS_ID[i]
    nuts_ID_subset <- subset(nuts_3_TBEV, nuts_3_TBEV$NUTS_ID == nuts_ID)
    nuts_ID_r <- terra::rasterize(nuts_ID_subset, europe_mask_50km, field = "num_cells", touches = TRUE)
    
    count_cells <- sum(values(nuts_ID_r) == num_cells, na.rm = TRUE)
    random_cell_index <- sample(1:count_cells, 1)
    
    cell_number <- which(values(nuts_ID_r) == num_cells)[random_cell_index]
    
    sampled_coordinates <- xyFromCell(nuts_ID_r, cell_number)
    
    nuts_3_TBEV$x[i] <- sampled_coordinates[1]
    nuts_3_TBEV$y[i] <- sampled_coordinates[2]
    
  }
} # Close the loop over all rows

# Save the data frame stating all infection occurrences and their respective
# location (for later usage in absence generation)
save(nuts_3_TBEV, file = "output_data/data/nuts_3_TBEV.RData")

# Remove entries that stem from municipalities that consist of more than 5 cells
# as this increases the uncertainty of the reported location
nuts_3_TBEV_filtered <- nuts_3_TBEV %>%
  anti_join(
    nuts_3_cell_count_df %>% filter(num_cells > 5), 
    by = c("x", "y")
  )

# Select only relevant columns of the data frame
TBEV_occurrences <- nuts_3_TBEV_filtered %>%
  rowid_to_column(var = "occ_id") %>% # create unique identifier for each occurrence
  dplyr::select(occ_id, x, y, DateUsedForStatisticsYear, DateUsedForStatisticsMonth, NUTS_ID)

# Remove the multipolygon column
TBEV_occurrences <- st_drop_geometry(TBEV_occurrences)

# Change the name of the columns
colnames(TBEV_occurrences) <- c("occ_id", "lon", "lat", "year", "month", "NUTS_ID")

# Add a column indicating the presence of infection
TBEV_occurrences$occ <- 1

# Convert the occurrence data frame to a spatial object
TBEV_occurrences_sp <- st_as_sf(TBEV_occurrences, coords = c("lon", "lat"), crs = st_crs(europe_mask_50km))

# Extract raster values at the coordinate locations
occurrences_values <- terra::extract(europe_mask_50km, TBEV_occurrences_sp)

# Only keep the occurrences where the values is 1 (meaning that it lays on the European continent)
TBEV_occurrences_cleaned <- TBEV_occurrences[!is.na(occurrences_values[,"layer"]) & occurrences_values[,"layer"] == 1, ]


# Plot the extracted infection occurrences
ggplot(nuts_3_cell_count_df, aes(x = x, y = y, fill = num_cells)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cell count\nper NUTS3 municipality", option = "viridis") +
  geom_point(data = TBEV_occurrences_cleaned, aes(x = lon, y = lat, color = as.factor(occ)), size = 2, inherit.aes = FALSE) +
  scale_color_manual(
    name = "TBEV infection",
    values = c("1" = "red")
  ) +
  theme_minimal() +
  labs(
    title = "Uncertainty of ECDC infection data on NUTS3 level & TBEV infection data",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold") 
  )




# Save the data frame with occurrence points
save(TBEV_occurrences_cleaned, file = "output_data/data/TBEV_occurrences_cleaned.RData")



#-------------------------------------------------------------------------------

# 3. Check countries that provided infection data ------------------------------

# Load the needed package
library(countrycode)

# Create a vector containing all reporting countries
reporting_countries <- unique(TBEV_infection_data$ReportingCountry)

# Get the country names based on ISO 2-Letter Code
reporting_countries_full <- countrycode(reporting_countries, origin = "iso2c", destination = "country.name")

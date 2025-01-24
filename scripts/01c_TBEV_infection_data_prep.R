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



#-------------------------------------------------------------------------------

# 1. Extracting the uncertainty of ECDC disease data ---------------------------

# ECDC human case infection data is provided at the NUTS3 level,
# Retain a map showing the European municipalities on NUTS 3 level (year 2024)
nuts_3 <- gisco_get_nuts(
  year = "2024",
  epsg = "4326",
  cache = TRUE,
  update_cache = FALSE,
  cache_dir = NULL,
  verbose = FALSE,
  resolution = "20",
  spatialtype = "RG",
  country = NULL,
  nuts_id = NULL,
  nuts_level = "3"
)

# Rasterise the NUTS3 multipolygon and mask the values that do not belong to the European continent
nuts_3_raster <- terra::rasterize(nuts_3, europe_mask_50km, field = "NUTS_ID", touches = TRUE)
nuts_3_raster_mask <- terra::mask(nuts_3_raster, europe_mask_50km)

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

# Load the infection data
TBEV_infection_data <- x

# Rename the column RegionCode to NUTS_ID
colnames(TBEV_infection_data)[colnames(TBEV_infection_data) == "RegionCode"] <- "NUTS_ID"

# Join the infection data with NUTS3 geographic information
nuts_3_TBEV <- nuts_3 %>%
  left_join(TBEV_infection_data, by = "NUTS_ID")

# Remove rows of municipalities that do not have any observed infection data 
nuts_3_TBEV <- nuts_3_TBEV[!is.na(nuts_3_TBEV$Population), ]

# Extract the centroid information of each NUTS3 municipality and append the 
# info to the data frame
nuts_3_TBEV$centroid <- st_centroid(nuts_3_TBEV$geometry)
nuts_3_TBEV$x <- st_coordinates(nuts_3_TBEV$centroid)[, 1]
nuts_3_TBEV$y <- st_coordinates(nuts_3_TBEV$centroid)[, 2]

# Remove all rows with entries of observed infection numbers below 1
nuts_3_TBEV <- nuts_3_TBEV[nuts_3_TBEV$NumValue > 1, ]

# Remove infection entries from from the years 2020 and higher as these are 
# not covered by environmental data
nuts_3_TBEV <- subset(nuts_3_TBEV, nuts_3_TBEV$Time < 2020)

# Remove infection entries within the same cell if occurring in the same 
# month of a certain year
nuts_3_TBEV <- nuts_3_TBEV %>%
  distinct(NUTS_ID, month, year, .keep_all = TRUE)

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



# Select only relevant columns of the data frame
TBEV_occurrences <- nuts_3_TBEV %>%
  rowid_to_column(var = "occ_id") %>% # create unique identifier for each occurrence
  dplyr::select(occ_id, disease, x, y, year, month, NUTS_ID)

# Plot the extracted infection occurrences
ggplot(nuts_3_cell_count_df, aes(x = x, y = y, fill = num_cells)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Cell count\nper NUTS3 municipality", option = "viridis") +
  geom_point(data = TBEV_occurrences, aes(x = x, y = y, color = occ), size = 2, inherit.aes = FALSE) +
  scale_color_manual(
    name = "TBEV infection",
    values = c("presence" = "red")
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
save(TBEV_occurrences, file = "output_data/data/TBEV_occurrences.RData")



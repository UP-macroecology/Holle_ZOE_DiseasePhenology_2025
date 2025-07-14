# ZOE disease phenology analysis


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#             Past and future duration trends (extra script)             #
# ---------------------------------------------------------------------- #








#-------------------------------------------------------------------------------

# 5. Visualise duration trends  ------------------------------------------------
# for main vectors as well as viruses
# showing future duration trends trends (hotspots) and marking the cells
# borders of the case study regions

# Load needed data
# The masks containing the countries with a mandatory reporting system of 
# disease surveillance data to the ECDC; we load different masks depending on
# the virus as Austria did not report of NUTS3 level for TBE
eu_eea_mask_WNV <- terra::rast("input_data/spatial_data/eu_eea_mask_WNV.tif") 
eu_eea_mask_TBE <- terra::rast("input_data/spatial_data/eu_eea_mask_TBE.tif")

# Load the needed data
load("output_data/results/duration_trends_fut_vector_virus_ssp370_new.RData")


# Read in the shapefiles from the case study regions Slovenia and Slovakia
shp_slovenia <- terra::vect("input_data/spatial_data/si_1km.shp")
shp_slovakia <- terra::vect("input_data/spatial_data/sk_1km.shp")

# Reproject the shapefile to the used crs (WGS84)
shp_slovenia <- project(shp_slovenia, "EPSG:4326")
shp_slovakia <- project(shp_slovakia, "EPSG:4326")

# Convert these shapefiles into rasters with a 50km resolution
r_slovenia <- rasterize(shp_slovenia, eu_eea_mask_WNV)
r_slovakia <- rasterize(shp_slovakia, eu_eea_mask_WNV)

# Crop the study region rasters to match the extent of the prediction rasters
# based on one example raster
r_slovenia <- crop(r_slovenia, I_ricinus_occ_prob_clim_landuse_ens_bin)
r_slovakia <- crop(r_slovakia, I_ricinus_occ_prob_clim_landuse_ens_bin)

# Make sure that raster cells of the case study regions are masked when
# an intersection with other non-studied cell (NA) is occurring (e.g. Austria)
# and set the respective cell to NA (take one prediction raster as template)
r_slovenia[which(is.na(values(I_ricinus_occ_prob_clim_landuse_ens_bin[[1]])))] <- NA
r_slovakia[which(is.na(values(I_ricinus_occ_prob_clim_landuse_ens_bin[[1]])))] <- NA

# Merge adjacent cells into a single polygon and convert into an sf object 
# for plotting
r_slovenia_poly <- as.polygons(r_slovenia, dissolve = TRUE)
r_slovenia_poly <- st_as_sf(r_slovenia_poly)
r_slovakia_poly <- as.polygons(r_slovakia, dissolve = TRUE)
r_slovakia_poly <- st_as_sf(r_slovakia_poly)



# Make sure the vector and diseases appear in the correct position
slope_df_futssp370$species <- factor(slope_df_futssp370$species, 
                                     levels = c("Ixodes ricinus", "TBE", "Culex pipiens", "WNV"))

# Create two different datasets 
duration_trends_past_futssp370_Ixodes_TBE <- slope_df_futssp370[slope_df_futssp370$species %in% c("Ixodes ricinus", "TBE"), ]
duration_trends_past_futssp370_Culex_WNV <- slope_df_futssp370[slope_df_futssp370$species %in% c("Culex pipiens", "WNV"), ]

# Convert Europe mask spatraster into a data frame
europe_mask_50km_df <- as.data.frame(europe_mask_50km, xy = TRUE)





# Visualize the future duration trends of Ixodes ricinus and TBE
# and highlight the cell borders of the case study regions
ggplot() +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster(data = duration_trends_past_futssp370_Ixodes_TBE,
              aes(x = lon, y = lat, fill = trend)) +
  geom_sf(data = r_slovenia_poly, fill = NA, color = "firebrick4", linewidth = 0.8) +
  geom_sf(data = r_slovakia_poly, fill = NA, color = "firebrick4", linewidth = 0.8) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        breaks = c(-0.2, 0, 0.2)) +
  facet_grid2(rows = vars(species), cols = vars(time), scales = "fixed",
              strip = strip_themed(background_y = list(
                "Ixodes ricinus" = element_rect(fill = "steelblue3"),
                "TBE" = element_rect(fill = "steelblue3")
              ))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey75", color = NA)
  )

ggsave(paste0("output_data/plots/duration_trends/duration_trends_2010s_Ixodes_TBE_ssp370_studyregions.png"), width = 5.0, height = 6.5)



# Visualize the future duration trends of Culex pipiens and WNV
# and highlight the cell borders of the case study regions
ggplot() +
  geom_raster(data = europe_mask_50km_df, aes(x = x, y = y), fill = "gray22") +
  geom_raster(data = duration_trends_past_futssp370_Culex_WNV,
              aes(x = lon, y = lat, fill = trend)) +
  geom_sf(data = r_slovenia_poly, fill = NA, color = "firebrick4", linewidth = 0.8) +
  geom_sf(data = r_slovakia_poly, fill = NA, color = "firebrick4", linewidth = 0.8) +
  scale_fill_whitebox_c(palette = "muted", name = "Duration trend\n(months per decade)",
                        breaks = c(-0.2, 0, 0.2)) +
  facet_grid2(rows = vars(species), cols = vars(time),
              strip = strip_themed(background_y = list(
                "Culex pipiens" = element_rect(fill = "lightsteelblue1"),
                "WNV" = element_rect(fill = "lightsteelblue1")
              ))) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    strip.background = element_rect(fill = "grey75", color = NA)
  )

ggsave(paste0("output_data/plots/duration_trends/duration_trends_2010s_Culex_WNV_ssp370_studyregions.png"), width = 5.0, height = 6.5)

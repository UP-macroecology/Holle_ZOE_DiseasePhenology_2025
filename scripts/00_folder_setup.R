# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                      00. Set up folder structure                       #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We create all necessary directories for the project, including input data, 
# output data, and scripts folders.


#-------------------------------------------------------------------------------

# 1. Create project folder structure -------------------------------------------



dirs <- c(
  # Base directory: scripts
  "scripts",
  
  # Base directory: input data - environmental data  
  # Historical climate and land use data
  "input_data/environmental_data/ISIMIP3a/Climate/processed_data",
  "input_data/environmental_data/ISIMIP3a/Climate/raw_data",
  "input_data/environmental_data/ISIMIP3a/CounterClim/processed_data",
  "input_data/environmental_data/ISIMIP3a/CounterClim/raw_data",
  "input_data/environmental_data/ISIMIP3a/CounterLandUse/processed_data",
  "input_data/environmental_data/ISIMIP3a/CounterLandUse/raw_data",
  "input_data/environmental_data/ISIMIP3a/LandUse/processed_data",
  "input_data/environmental_data/ISIMIP3a/LandUse/raw_data",
  
  # Base directory: input data - environmental data 
  # Future climate and land use data
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/processed_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/processed_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/processed_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/processed_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/processed_data/ukesm1-0-ll",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/raw_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/raw_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/raw_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/raw_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp126/raw_data/ukesm1-0-ll",
  
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/processed_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/processed_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/processed_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/processed_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/processed_data/ukesm1-0-ll",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/raw_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/raw_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/raw_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/raw_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp370/raw_data/ukesm1-0-ll",
  
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/processed_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/processed_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/processed_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/processed_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/processed_data/ukesm1-0-ll",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/raw_data/gfdl-esm4",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/raw_data/ipsl-cm6a-lr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/raw_data/mpi-esm1-2-hr",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/raw_data/mri-esm2-0",
  "input_data/environmental_data/ISIMIP3b/Climate/ssp585/raw_data/ukesm1-0-ll",
  
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp126/processed_data",
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp126/raw_data",
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp370/processed_data",
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp370/raw_data",
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp585/processed_data",
  "input_data/environmental_data/ISIMIP3b/LandUse/ssp585/raw_data",
  
  # Base directory: input data - raw infection data
  "input_data/raw_infection_data",
  
  # Base directory: input data - raw species data
  "input_data/raw_species_data",
  
  # Base directory: input data - spatial data
  "input_data/spatial_data/climate_regions",
  
  # Base directory: output data - data
  "output_data/data",
  
  # Base directory: output data - models
  "output_data/models",
  
  # Base directory: output data - validation
  "output_data/validation",
  
  # Base directory: output data - results
  "output_data/results/preprocessed_predictions/Culex_pipiens",
  "output_data/results/preprocessed_predictions/Ixodes_ricinus",
  "output_data/results/preprocessed_predictions/WNV",
  "output_data/results/preprocessed_predictions/TBE",
  "output_data/results/postprocessed_predictions/Culex_pipiens",
  "output_data/results/postprocessed_predictions/Ixodes_ricinus",
  "output_data/results/postprocessed_predictions/WNV",
  "output_data/results/postprocessed_predictions/TBE",
  "output_data/results/decadal_trends",
  "output_data/results/duration_trends",
  "output_data/results/climateregions_trends",
  "output_data/results/distribution_trends",
  
  # Base directory: output data - plots
  "output_data/plots/maps",
  "output_data/plots/presence_background",
  "output_data/plots/response_curves",
  "output_data/plots/overview",
  "output_data/plots/decadal_trends",
  "output_data/plots/duration_trends",
  "output_data/plots/climateregions_trends",
  "output_data/plots/distribution_trends"
)


# Create all directories if missing (no warnings if folder already exists)
for (d in dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("Folder structure created successfully!\n")
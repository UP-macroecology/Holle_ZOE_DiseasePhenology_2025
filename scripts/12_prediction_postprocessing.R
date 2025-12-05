# ZOE disease phenology analysis

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                    12. Post-processing predictions                     #
# ---------------------------------------------------------------------- #


# What is done within this script:

# We post-process the ensemble prediction rasters for the viruses (TBE, WNV) and 
# their primary vector species (Ixodes ricinus, Culex pipiens) by averaging the 
# future layers of the same scenario across the five different climate models. 
# Additionally. we mask the vector prediction rasters to include only 
# cells within EU/EEA countries to ensure comparability between vector and 
# virus results. For the virus layers, cells are set to 0 whenever the 
# corresponding vector species is not predicted to be present, ensuring for 
# ecological realism. Finally, all post-processed rasters for historical and 
# future predictions are saved for further analysis.



# Load needed packages
library(terra) # terra_1.7-55

# Load needed data
eu_eea_mask <- terra::rast("input_data/spatial_data/eu_eea_mask.tif") # Mask of EU/EEA countries at a 0.5° resolution



#-------------------------------------------------------------------------------

# 1. Generate averaged future predictions --------------------------------------
# Generate averaged future predictions per cell over the five different used climate 
# models, for each investigated environmental scenario, and each of the investigated 
# species/pathogens

# Load needed data - monthly predicted ensemble occurrence probability of main vectors 
# and the respective viruses, under climate and land use change for the
# future years 2020 to 2059 based on the five different climate models and three
# environmental scenarios
# Ixodes ricinus
# ssp126
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp126.tif")) # Ensemble predictions based on climate model 1 under env. scenario ssp126
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp126.tif")) # Climate model 2
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp126.tif")) # Climate model 3
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp126.tif")) # Climate model 4
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp126.tif")) # Climate model 5
# ssp370
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp370.tif"))
# ssp585
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp585.tif"))

# Culex pipiens
# spp126
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp126.tif"))
# ssp370
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp370.tif"))
# ssp585
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp585.tif"))


# TBE
# ssp126
TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp126.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp126.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp126.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp126.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp126.tif"))
# ssp370
TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp370.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp370.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp370.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp370.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp370.tif"))
# ssp585
TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp585.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp585.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp585.tif"))
TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp585.tif"))
TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp585.tif"))

# WNV
# ssp126
WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp126.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp126.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp126.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp126.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp126.tif"))
# ssp370
WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp370.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp370.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp370.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp370.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp370.tif"))
# ssp585
WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_gfdl-esm4_ssp585.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ipsl-cm6a-lr_ssp585.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mpi-esm1-2-hr_ssp585.tif"))
WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_mri-esm2-0_ssp585.tif"))
WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ukesm1-0-ll_ssp585.tif"))



# Average the predictions per cell over the five different climate models,
# do this for each environmental scenario and for each studied species/pathogen
# Ixodes ricinus
# ssp126
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mean(I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mean(I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mean(I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585,
                                                              I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585)

# Culex pipiens
# ssp126
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mean(C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mean(C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mean(C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585,
                                                              C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585)

# TBE
# ssp126
TBE_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mean(TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
TBE_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mean(TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
TBE_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mean(TBE_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585,
                                                        TBE_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585,
                                                        TBE_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585)

# WNV
# ssp126
WNV_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mean(WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
WNV_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mean(WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
WNV_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mean(WNV_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585,
                                                        WNV_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585,
                                                        WNV_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585)





#-------------------------------------------------------------------------------

# 2. Mask species prediction rasters (to EU/EEA countries) ---------------------
# Adjust the extent of the mask and mask the cells with predicted values of countries
# that do not belong to the EU/EEA for the two main vector species
# Mask cells that are not within the EU/EEA countries - to make results comparable
# between vectors and viruses


# (a) Predictions based on historical data -------------------------------------

# Load the needed data - monthly predicted ensemble occurrence probability of 
# main vectors and the respective viruses, under climate and land use change, 
# as well as under the counterfactual scenarios for the years 1970 to 2019
 
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_1970_2019.tif") # Based on climate and land use change
I_ricinus_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_landuse_ens_1970_2019.tif") # Based on counterfactual climate and factual land use change
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_nolanduse_ens_1970_2019.tif") # Based on counterfactual land use and factual climate change
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_nolanduse_ens_1970_2019.tif") # Based on counterfactual climate and counterfactual land use

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif")

# TBE
TBE_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_1970_2019.tif")
TBE_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE/TBE_preds_noclim_landuse_ens_1970_2019.tif")
TBE_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE/TBE_preds_clim_nolanduse_ens_1970_2019.tif")
TBE_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE/TBE_preds_noclim_nolanduse_ens_1970_2019.tif")

# WNV
WNV_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_1970_2019.tif")
WNV_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV/WNV_preds_noclim_landuse_ens_1970_2019.tif")
WNV_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV/WNV_preds_clim_nolanduse_ens_1970_2019.tif")
WNV_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV/WNV_preds_noclim_nolanduse_ens_1970_2019.tif")



# Adjust the extent of the mask first, to fit the raster extent of the
# species/pathogen predictions (use one prediction raster as template/example)
eu_eea_mask <- terra::crop(eu_eea_mask, I_ricinus_occ_prob_clim_landuse_ens)


# Mask the values in the cells of the species predicitons that do not belong to
# countries within the EU/EEA
I_ricinus_occ_prob_clim_landuse_ens <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens, eu_eea_mask) # Ixodes ricinus
I_ricinus_occ_prob_noclim_landuse_ens <- terra::mask(I_ricinus_occ_prob_noclim_landuse_ens, eu_eea_mask) # Ixodes ricinus
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::mask(I_ricinus_occ_prob_clim_nolanduse_ens, eu_eea_mask) # Ixodes ricinus
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::mask(I_ricinus_occ_prob_noclim_nolanduse_ens, eu_eea_mask) # Ixodes ricinus

C_pipiens_occ_prob_clim_landuse_ens <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens, eu_eea_mask) # Culex pipiens
C_pipiens_occ_prob_noclim_landuse_ens <- terra::mask(C_pipiens_occ_prob_noclim_landuse_ens, eu_eea_mask) # Culex pipiens
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::mask(C_pipiens_occ_prob_clim_nolanduse_ens, eu_eea_mask) # Culex pipiens
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::mask(C_pipiens_occ_prob_noclim_nolanduse_ens, eu_eea_mask) # Culex pipiens


# Adjust the extent of the rasters again to match the raster extent of the
# virus predictions
I_ricinus_occ_prob_clim_landuse_ens <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens, TBE_occ_prob_clim_landuse_ens)
I_ricinus_occ_prob_noclim_landuse_ens <- terra::crop(I_ricinus_occ_prob_noclim_landuse_ens, TBE_occ_prob_noclim_landuse_ens)
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::crop(I_ricinus_occ_prob_clim_nolanduse_ens, TBE_occ_prob_clim_nolanduse_ens)
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::crop(I_ricinus_occ_prob_noclim_nolanduse_ens, TBE_occ_prob_noclim_nolanduse_ens)

C_pipiens_occ_prob_clim_landuse_ens <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens, WNV_occ_prob_clim_landuse_ens)
C_pipiens_occ_prob_noclim_landuse_ens <- terra::crop(C_pipiens_occ_prob_noclim_landuse_ens, WNV_occ_prob_noclim_landuse_ens)
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::crop(C_pipiens_occ_prob_clim_nolanduse_ens, WNV_occ_prob_clim_nolanduse_ens)
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::crop(C_pipiens_occ_prob_noclim_nolanduse_ens, WNV_occ_prob_noclim_nolanduse_ens)




# (b) Predictions based on future data -----------------------------------------

# Adjust the extent of the mask first, to fit the raster extent of the
# species/pathogen predictions (use one prediction raster as template/example)
eu_eea_mask <- terra::crop(eu_eea_mask, I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126)


# Mask the values in the cells of the species predicitons that do not belong to
# countries within the EU/EEA (for the different species based on the
# three different studied environmental scenarios)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, eu_eea_mask) # Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, eu_eea_mask) # Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, eu_eea_mask) # Ixodes ricinus

C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, eu_eea_mask) # Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, eu_eea_mask) # Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, eu_eea_mask) # Culex pipiens



# Adjust the extent of the rasters again to match the raster extent of the
# virus predictions
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, TBE_occ_prob_clim_landuse_ens_fut_ssp126)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, TBE_occ_prob_clim_landuse_ens_fut_ssp370)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, TBE_occ_prob_clim_landuse_ens_fut_ssp585)

C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, WNV_occ_prob_clim_landuse_ens_fut_ssp126)
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, WNV_occ_prob_clim_landuse_ens_fut_ssp370)
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, WNV_occ_prob_clim_landuse_ens_fut_ssp585)




#-------------------------------------------------------------------------------

# 3. Binarise post-processed continuous prediction data  -----------------------

# Load the validation data of the models to extract the ensemble thresholds
# to obtain binary predictions from averaged ensemble continuous predictions
# for the different rasters
load("output_data/validation/I_ricinus_validation.RData") # Load validation results
I_ricinus_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"] # Extract threshold

load("output_data/validation/C_pipiens_validation.RData")
C_pipiens_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]

load("output_data/validation/TBE_validation.RData") 
TBE_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]

load("output_data/validation/WNV_validation.RData")
WNV_thresh <- comp_perf[comp_perf$alg == "mean_prob", "thresh"]



# (a) Predictions based on historical data -------------------------------------

# Apply the threshold to all predictions cells of prediction rasters for
# each species and the different scenarios
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_bin <- I_ricinus_occ_prob_clim_landuse_ens >= I_ricinus_thresh # Apply corresponding threshold
I_ricinus_occ_prob_clim_landuse_ens_bin <- I_ricinus_occ_prob_clim_landuse_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

I_ricinus_occ_prob_noclim_landuse_ens_bin <- I_ricinus_occ_prob_noclim_landuse_ens >= I_ricinus_thresh
I_ricinus_occ_prob_noclim_landuse_ens_bin <- I_ricinus_occ_prob_noclim_landuse_ens_bin * 1

I_ricinus_occ_prob_clim_nolanduse_ens_bin <- I_ricinus_occ_prob_clim_nolanduse_ens >= I_ricinus_thresh
I_ricinus_occ_prob_clim_nolanduse_ens_bin <- I_ricinus_occ_prob_clim_nolanduse_ens_bin * 1

I_ricinus_occ_prob_noclim_nolanduse_ens_bin <- I_ricinus_occ_prob_noclim_nolanduse_ens >= I_ricinus_thresh
I_ricinus_occ_prob_noclim_nolanduse_ens_bin <- I_ricinus_occ_prob_noclim_nolanduse_ens_bin * 1

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_bin <- C_pipiens_occ_prob_clim_landuse_ens >= C_pipiens_thresh # Apply corresponding threshold
C_pipiens_occ_prob_clim_landuse_ens_bin <- C_pipiens_occ_prob_clim_landuse_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

C_pipiens_occ_prob_noclim_landuse_ens_bin <- C_pipiens_occ_prob_noclim_landuse_ens >= C_pipiens_thresh
C_pipiens_occ_prob_noclim_landuse_ens_bin <- C_pipiens_occ_prob_noclim_landuse_ens_bin * 1

C_pipiens_occ_prob_clim_nolanduse_ens_bin <- C_pipiens_occ_prob_clim_nolanduse_ens >= C_pipiens_thresh
C_pipiens_occ_prob_clim_nolanduse_ens_bin <- C_pipiens_occ_prob_clim_nolanduse_ens_bin * 1

C_pipiens_occ_prob_noclim_nolanduse_ens_bin <- C_pipiens_occ_prob_noclim_nolanduse_ens >= C_pipiens_thresh
C_pipiens_occ_prob_noclim_nolanduse_ens_bin <- C_pipiens_occ_prob_noclim_nolanduse_ens_bin * 1




# Ensure for ecological realism and set cells of virus predictions to 0 when 
# vector is not occurring there
# TBE based on predicted presence/absence of Ixodes ricinus
TBE_occ_prob_clim_landuse_ens <- TBE_occ_prob_clim_landuse_ens * I_ricinus_occ_prob_clim_landuse_ens_bin
TBE_occ_prob_noclim_landuse_ens <- TBE_occ_prob_noclim_landuse_ens * I_ricinus_occ_prob_noclim_landuse_ens_bin
TBE_occ_prob_clim_nolanduse_ens <- TBE_occ_prob_clim_nolanduse_ens * I_ricinus_occ_prob_clim_nolanduse_ens_bin
TBE_occ_prob_noclim_nolanduse_ens <- TBE_occ_prob_noclim_nolanduse_ens * I_ricinus_occ_prob_noclim_nolanduse_ens_bin

# WNV based on predicted presence/absence of Culex pipiens
WNV_occ_prob_clim_landuse_ens <- WNV_occ_prob_clim_landuse_ens * C_pipiens_occ_prob_clim_landuse_ens_bin
WNV_occ_prob_noclim_landuse_ens <- WNV_occ_prob_noclim_landuse_ens * C_pipiens_occ_prob_noclim_landuse_ens_bin
WNV_occ_prob_clim_nolanduse_ens <- WNV_occ_prob_clim_nolanduse_ens * C_pipiens_occ_prob_clim_nolanduse_ens_bin
WNV_occ_prob_noclim_nolanduse_ens <- WNV_occ_prob_noclim_nolanduse_ens * C_pipiens_occ_prob_noclim_nolanduse_ens_bin



# Apply the threshold to all predictions cells of prediction rasters for
# each virus and the different scenarios
# TBE
TBE_occ_prob_clim_landuse_ens_bin <- TBE_occ_prob_clim_landuse_ens >= TBE_thresh # Apply corresponding threshold
TBE_occ_prob_clim_landuse_ens_bin <- TBE_occ_prob_clim_landuse_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

TBE_occ_prob_noclim_landuse_ens_bin <- TBE_occ_prob_noclim_landuse_ens >= TBE_thresh
TBE_occ_prob_noclim_landuse_ens_bin <- TBE_occ_prob_noclim_landuse_ens_bin * 1

TBE_occ_prob_clim_nolanduse_ens_bin <- TBE_occ_prob_clim_nolanduse_ens >= TBE_thresh
TBE_occ_prob_clim_nolanduse_ens_bin <- TBE_occ_prob_clim_nolanduse_ens_bin * 1

TBE_occ_prob_noclim_nolanduse_ens_bin <- TBE_occ_prob_noclim_nolanduse_ens >= TBE_thresh
TBE_occ_prob_noclim_nolanduse_ens_bin <- TBE_occ_prob_noclim_nolanduse_ens_bin * 1

# WNV
WNV_occ_prob_clim_landuse_ens_bin <- WNV_occ_prob_clim_landuse_ens >= WNV_thresh # Apply corresponding threshold
WNV_occ_prob_clim_landuse_ens_bin <- WNV_occ_prob_clim_landuse_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

WNV_occ_prob_noclim_landuse_ens_bin <- WNV_occ_prob_noclim_landuse_ens >= WNV_thresh
WNV_occ_prob_noclim_landuse_ens_bin <- WNV_occ_prob_noclim_landuse_ens_bin * 1

WNV_occ_prob_clim_nolanduse_ens_bin <- WNV_occ_prob_clim_nolanduse_ens >= WNV_thresh
WNV_occ_prob_clim_nolanduse_ens_bin <- WNV_occ_prob_clim_nolanduse_ens_bin * 1

WNV_occ_prob_noclim_nolanduse_ens_bin <- WNV_occ_prob_noclim_nolanduse_ens >= WNV_thresh
WNV_occ_prob_noclim_nolanduse_ens_bin <- WNV_occ_prob_noclim_nolanduse_ens_bin * 1





# (b) Predictions based on future data -----------------------------------------

# Apply the threshold to all prediction cells of prediction rasters for
# each species and the different studied future env. scenarios
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 >= I_ricinus_thresh # Apply corresponding threshold
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126 * 1 # Convert TRUE/FALSE to 1/0 values

I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 >= I_ricinus_thresh
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370 * 1

I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 >= I_ricinus_thresh
I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585 * 1

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 >= C_pipiens_thresh
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126 * 1

C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 >= C_pipiens_thresh
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370 * 1

C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 >= C_pipiens_thresh
C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585 * 1




# Ensure for ecological realism and set cells of virus predictions to 0 when 
# vector is not occurring there
# TBE based on predicted presence/absence of Ixodes ricinus
TBE_occ_prob_clim_landuse_ens_fut_ssp126 <- TBE_occ_prob_clim_landuse_ens_fut_ssp126 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126
TBE_occ_prob_clim_landuse_ens_fut_ssp370 <- TBE_occ_prob_clim_landuse_ens_fut_ssp370 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370
TBE_occ_prob_clim_landuse_ens_fut_ssp585 <- TBE_occ_prob_clim_landuse_ens_fut_ssp585 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585

# WNV based on predicted presence/absence of Culex pipiens
WNV_occ_prob_clim_landuse_ens_fut_ssp126 <- WNV_occ_prob_clim_landuse_ens_fut_ssp126 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126
WNV_occ_prob_clim_landuse_ens_fut_ssp370 <- WNV_occ_prob_clim_landuse_ens_fut_ssp370 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370
WNV_occ_prob_clim_landuse_ens_fut_ssp585 <- WNV_occ_prob_clim_landuse_ens_fut_ssp585 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585





# Apply the threshold to all prediction cells of prediction rasters for
# each virus and the different scenarios
# TBE
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- TBE_occ_prob_clim_landuse_ens_fut_ssp126 >= TBE_thresh
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126 * 1

TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- TBE_occ_prob_clim_landuse_ens_fut_ssp370 >= TBE_thresh
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370 * 1

TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- TBE_occ_prob_clim_landuse_ens_fut_ssp585 >= TBE_thresh
TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585 * 1

# WNV
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_landuse_ens_fut_ssp126 >= WNV_thresh
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126 * 1

WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_landuse_ens_fut_ssp370 >= WNV_thresh
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370 * 1

WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_landuse_ens_fut_ssp585 >= WNV_thresh
WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585 * 1






#-------------------------------------------------------------------------------

# 5. Save post-processed prediction rasters ------------------------------------
# Saving continuous and binary predictions of different past and future scenarios

# (a) Predictions based on historical data -------------------------------------

# Continuous and binary predictions of different scenarios
# Ixodes ricinus
writeRaster(I_ricinus_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# Culex pipiens
writeRaster(C_pipiens_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# TBE
writeRaster(TBE_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/TBE/TBE_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/TBE/TBE_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(TBE_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/TBE/TBE_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/TBE/TBE_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# WNV
writeRaster(WNV_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/WNV/WNV_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/WNV/WNV_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(WNV_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/WNV/WNV_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/WNV/WNV_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)




# (b) Predictions based on future data -----------------------------------------

# Ixodes ricinus
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_2020_2059_ssp585.tif", overwrite = TRUE)

writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/Ixodes_ricinus/I_ricinus_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif", overwrite = TRUE)

# Culex pipiens
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_2020_2059_ssp585.tif", overwrite = TRUE)

writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/Culex_pipiens/C_pipiens_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif", overwrite = TRUE)

# TBE
writeRaster(TBE_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_2020_2059_ssp585.tif", overwrite = TRUE)

writeRaster(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/TBE/TBE_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif", overwrite = TRUE)

# WNV
writeRaster(WNV_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_2020_2059_ssp585.tif", overwrite = TRUE)

writeRaster(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp126.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp370.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/WNV/WNV_preds_clim_landuse_ens_bin_2020_2059_ssp585.tif", overwrite = TRUE)



# ZOE disease phenology analysis

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                    12. Post-processing predictions                     #
# ---------------------------------------------------------------------- #


# Load needed packages
library(terra)

# Load needed data
# The masks containing the countries with a mandatory reporting system of 
# disease surveillance data to the ECDC; we load different masks depending on
# the virus as Austria did not report of NUTS3 level for TBE
eu_eea_mask_WNV <- terra::rast("input_data/spatial_data/eu_eea_mask_WNV.tif") 
eu_eea_mask_TBE <- terra::rast("input_data/spatial_data/eu_eea_mask_TBE.tif")


#-------------------------------------------------------------------------------

# 1. Generate averaged future predictions --------------------------------------
# Generate averaged future predictions per cell over the five different used climate 
# models, for each investigated environmental scenario, and each of the investigated 
# species/pathogens

# Load needed data - monthly predicted ensemble occurrence probability of main vectors 
# and the respective viruses, under climate and land use change for the
# future years 2030 to 2070 based on the five different climate models and three
# environmental scenarios
# Ixodes ricinus
# ssp126
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp126.tif")) # Ensemble predictions based on climate model 1 under env. scenario ssp126
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp126.tif")) # Climate model 2
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp126.tif")) # Climate model 3
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp126.tif")) # Climate model 4
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp126.tif")) # Climate model 5
# ssp370
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp370.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp370.tif"))
# ssp585
I_ricinus_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp585.tif"))
I_ricinus_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp585.tif"))

# Culex pipiens
# spp126
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp126.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp126.tif"))
# ssp370
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp370.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp370.tif"))
# ssp585
C_pipiens_occ_prob_clim_landuse_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_gfdl-esm4_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ipsl-cm6a-lr_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mpi-esm1-2-hr_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_mri-esm2-0_ssp585.tif"))
C_pipiens_occ_prob_clim_landuse_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ukesm1-0-ll_ssp585.tif"))


# TBE
# ssp126
TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_gfdl-esm4_ssp126.tif"))
TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp126.tif"))
TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp126.tif"))
TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mri-esm2-0_ssp126.tif"))
TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp126.tif"))
# ssp370
TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_gfdl-esm4_ssp370.tif"))
TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp370.tif"))
TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp370.tif"))
TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mri-esm2-0_ssp370.tif"))
TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp370.tif"))
# ssp585
TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_gfdl-esm4_ssp585.tif"))
TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp585.tif"))
TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp585.tif"))
TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_mri-esm2-0_ssp585.tif"))
TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp585.tif"))

# WNV
# ssp126
WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_gfdl-esm4_ssp126.tif"))
WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp126.tif"))
WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp126.tif"))
WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mri-esm2-0_ssp126.tif"))
WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp126 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp126.tif"))
# ssp370
WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_gfdl-esm4_ssp370.tif"))
WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp370.tif"))
WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp370.tif"))
WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mri-esm2-0_ssp370.tif"))
WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp370 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp370.tif"))
# ssp585
WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_gfdl-esm4_ssp585.tif"))
WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ipsl-cm6a-lr_ssp585.tif"))
WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mpi-esm1-2-hr_ssp585.tif"))
WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_mri-esm2-0_ssp585.tif"))
WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp585 <- terra::rast(paste0("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_2030_2070_ukesm1-0-ll_ssp585.tif"))



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
TBE_occ_prob_clim_ens_fut_ssp126 <- terra::mean(TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp126,
                                                TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp126,
                                                TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp126,
                                                TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp126,
                                                TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
TBE_occ_prob_clim_ens_fut_ssp370 <- terra::mean(TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp370,
                                                TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp370,
                                                TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp370,
                                                TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp370,
                                                TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
TBE_occ_prob_clim_ens_fut_ssp585 <- terra::mean(TBE_occ_prob_clim_ens_fut_gfdl_esm4_ssp585,
                                                TBE_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp585,
                                                TBE_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp585,
                                                TBE_occ_prob_clim_ens_fut_mri_esm2_0_ssp585,
                                                TBE_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp585)

# WNV
# ssp126
WNV_occ_prob_clim_ens_fut_ssp126 <- terra::mean(WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp126,
                                                WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp126,
                                                WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp126,
                                                WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp126,
                                                WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp126)
# ssp370
WNV_occ_prob_clim_ens_fut_ssp370 <- terra::mean(WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp370,
                                                WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp370,
                                                WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp370,
                                                WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp370,
                                                WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp370)
# ssp585
WNV_occ_prob_clim_ens_fut_ssp585 <- terra::mean(WNV_occ_prob_clim_ens_fut_gfdl_esm4_ssp585,
                                                WNV_occ_prob_clim_ens_fut_ipsl_cm6a_lr_ssp585,
                                                WNV_occ_prob_clim_ens_fut_mpi_esm1_2_hr_ssp585,
                                                WNV_occ_prob_clim_ens_fut_mri_esm2_0_ssp585,
                                                WNV_occ_prob_clim_ens_fut_ukesm1_0_ll_ssp585)





#-------------------------------------------------------------------------------

# 2. Adjust extent of species prediction rasters  ------------------------------
# Adjust the extent of the mask and mask the cells with predicted values of countries
# that do not report disease surveillance data to the ECDC for the two main vector species
# Mask cells that are not within the EU/EEA countries - to make results comparable
# vector and virus of a disease 

# (a) Predictions based on historical data -------------------------------------

# Load the needed data - monthly predicted ensemble occurrence probability of 
# main vectors and the respective viruses, under climate and land use change, 
# as well as under the counterfactual scenarios for the years 1970 to 2019
 
# Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_landuse_ens_1970_2019.tif") # Based on climate and land use change
I_ricinus_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/I_ricinus_preds_noclim_landuse_ens_1970_2019.tif") # Based on counterfactual climate and land use change
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/I_ricinus_preds_clim_nolanduse_ens_1970_2019.tif") # Based on counterfactual land use and climate change
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/I_ricinus_preds_noclim_nolanduse_ens_1970_2019.tif") # Based on counterfactual climate and counterfactual land use

# Culex pipiens
C_pipiens_occ_prob_clim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_landuse_ens <- terra::rast("output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif")
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif")
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif")

# TBE
TBE_occ_prob_clim_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE_preds_clim_ens_1970_2019.tif")
TBE_occ_prob_noclim_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE_preds_noclim_ens_1970_2019.tif")
TBE_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE_preds_clim_nolanduse_ens_1970_2019.tif")
TBE_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/TBE_preds_noclim_nolanduse_ens_1970_2019.tif")

# WNV
WNV_occ_prob_clim_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV_preds_clim_ens_1970_2019.tif")
WNV_occ_prob_noclim_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV_preds_noclim_ens_1970_2019.tif")
WNV_occ_prob_clim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV_preds_clim_nolanduse_ens_1970_2019.tif")
WNV_occ_prob_noclim_nolanduse_ens <- terra::rast("output_data/results/preprocessed_predictions/WNV_preds_noclim_nolanduse_ens_1970_2019.tif")



# Adjust the extent of the mask first, to fit the raster extent of the
# species/pathogen predictions (use one prediction raster as template/example)
eu_eea_mask_TBE <- terra::crop(eu_eea_mask_TBE, I_ricinus_occ_prob_clim_landuse_ens)
eu_eea_mask_WNV <- terra::crop(eu_eea_mask_WNV, C_pipiens_occ_prob_clim_landuse_ens)


# Mask the values in the cells of the species predicitons that do not belong to
# countries with a mandatory reporting system to ECDC based on the corresponding
# mask
I_ricinus_occ_prob_clim_landuse_ens <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens, eu_eea_mask_TBE) # Ixodes ricinus
I_ricinus_occ_prob_noclim_landuse_ens <- terra::mask(I_ricinus_occ_prob_noclim_landuse_ens, eu_eea_mask_TBE) # Ixodes ricinus
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::mask(I_ricinus_occ_prob_clim_nolanduse_ens, eu_eea_mask_TBE) # Ixodes ricinus
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::mask(I_ricinus_occ_prob_noclim_nolanduse_ens, eu_eea_mask_TBE) # Ixodes ricinus

C_pipiens_occ_prob_clim_landuse_ens <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens, eu_eea_mask_WNV) # Culex pipiens
C_pipiens_occ_prob_noclim_landuse_ens <- terra::mask(C_pipiens_occ_prob_noclim_landuse_ens, eu_eea_mask_WNV) # Culex pipiens
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::mask(C_pipiens_occ_prob_clim_nolanduse_ens, eu_eea_mask_WNV) # Culex pipiens
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::mask(C_pipiens_occ_prob_noclim_nolanduse_ens, eu_eea_mask_WNV) # Culex pipiens


# Adjust the extent of the rasters again to match the raster extent of the
# virus predictions
I_ricinus_occ_prob_clim_landuse_ens <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens, TBE_occ_prob_clim_ens)
I_ricinus_occ_prob_noclim_landuse_ens <- terra::crop(I_ricinus_occ_prob_noclim_landuse_ens, TBE_occ_prob_noclim_ens)
I_ricinus_occ_prob_clim_nolanduse_ens <- terra::crop(I_ricinus_occ_prob_clim_nolanduse_ens, TBE_occ_prob_clim_nolanduse_ens)
I_ricinus_occ_prob_noclim_nolanduse_ens <- terra::crop(I_ricinus_occ_prob_noclim_nolanduse_ens, TBE_occ_prob_noclim_nolanduse_ens)

C_pipiens_occ_prob_clim_landuse_ens <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens, WNV_occ_prob_clim_ens)
C_pipiens_occ_prob_noclim_landuse_ens <- terra::crop(C_pipiens_occ_prob_noclim_landuse_ens, WNV_occ_prob_noclim_ens)
C_pipiens_occ_prob_clim_nolanduse_ens <- terra::crop(C_pipiens_occ_prob_clim_nolanduse_ens, WNV_occ_prob_clim_nolanduse_ens)
C_pipiens_occ_prob_noclim_nolanduse_ens <- terra::crop(C_pipiens_occ_prob_noclim_nolanduse_ens, WNV_occ_prob_noclim_nolanduse_ens)




# (b) Predictions based on future data -----------------------------------------

# Adjust the extent of the mask first, to fit the raster extent of the
# species/pathogen predictions (use one prediction raster as template/example)
eu_eea_mask_TBE <- terra::crop(eu_eea_mask_TBE, I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126)
eu_eea_mask_WNV <- terra::crop(eu_eea_mask_WNV, C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126)

# Mask the values in the cells of the species predicitons that do not belong to
# countries with a mandatory reporting system to ECDC based on the corresponding
# mask (for the different species based on the
# three different studied environmental scenarios)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, eu_eea_mask_TBE) # Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, eu_eea_mask_TBE) # Ixodes ricinus
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mask(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, eu_eea_mask_TBE) # Ixodes ricinus

C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, eu_eea_mask_WNV) # Ixodes ricinus
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, eu_eea_mask_WNV) # Ixodes ricinus
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::mask(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, eu_eea_mask_WNV) # Ixodes ricinus



# Adjust the extent of the rasters again to match the raster extent of the
# virus predictions
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, TBE_occ_prob_clim_ens_fut_ssp126)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, TBE_occ_prob_clim_ens_fut_ssp370)
I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::crop(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, TBE_occ_prob_clim_ens_fut_ssp585)

C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, WNV_occ_prob_clim_ens_fut_ssp126)
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, WNV_occ_prob_clim_ens_fut_ssp370)
C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 <- terra::crop(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, WNV_occ_prob_clim_ens_fut_ssp585)






#-------------------------------------------------------------------------------

# 3. Multiplication of vector and virus predictions  ---------------------------
# Multiply the predicted suitability values of vectors with their
# corresponding viruses (based on the same scenario)
# this controls for ecological plausibility as virus suitability, therefore,
# cannot be higher than vector probability (virus is only there present, where
# the respective vector species occurs)


# (a) Predictions based on historical data -------------------------------------

# # TBE based on predicted suitability values of Ixodes ricinus
# TBE_occ_prob_clim_ens <- I_ricinus_occ_prob_clim_landuse_ens * TBE_occ_prob_clim_ens # Based on climate and land use change
# TBE_occ_prob_noclim_ens <- I_ricinus_occ_prob_noclim_landuse_ens * TBE_occ_prob_noclim_ens # Based on counterfactual climate and land use change
# TBE_occ_prob_clim_nolanduse_ens <- I_ricinus_occ_prob_clim_nolanduse_ens * TBE_occ_prob_clim_nolanduse_ens # Based on counterfactual land use and climate change
# TBE_occ_prob_noclim_nolanduse_ens <- I_ricinus_occ_prob_noclim_nolanduse_ens * TBE_occ_prob_noclim_nolanduse_ens # Based on counterfactual climate and counterfactual land use
# 
# # WNV based on predicted suitability values of Culex pipiens
# WNV_occ_prob_clim_ens <- C_pipiens_occ_prob_clim_landuse_ens * WNV_occ_prob_clim_ens 
# WNV_occ_prob_noclim_ens <- C_pipiens_occ_prob_noclim_landuse_ens * WNV_occ_prob_noclim_ens 
# WNV_occ_prob_clim_nolanduse_ens <- C_pipiens_occ_prob_clim_nolanduse_ens * WNV_occ_prob_clim_nolanduse_ens 
# WNV_occ_prob_noclim_nolanduse_ens <- C_pipiens_occ_prob_noclim_nolanduse_ens * WNV_occ_prob_noclim_nolanduse_ens 
# 
# 
# 
# 
# # (b) Predictions based on future data -----------------------------------------
# 
# # TBE based on predicted suitability values of Ixodes ricinus
# TBE_occ_prob_clim_ens_fut_ssp126 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126 * TBE_occ_prob_clim_ens_fut_ssp126 # Based on climate and land use change of env. scenario ssp126
# TBE_occ_prob_clim_ens_fut_ssp370 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370 * TBE_occ_prob_clim_ens_fut_ssp370 # Based on climate and land use change of env. scenario ssp370
# TBE_occ_prob_clim_ens_fut_ssp585 <- I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585 * TBE_occ_prob_clim_ens_fut_ssp585 # Based on climate and land use change of env. scenario ssp585
# 
# 
# # WNV based on predicted suitability values of Culex pipiens
# WNV_occ_prob_clim_ens_fut_ssp126 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126 * WNV_occ_prob_clim_ens_fut_ssp126
# WNV_occ_prob_clim_ens_fut_ssp370 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370 * WNV_occ_prob_clim_ens_fut_ssp370
# WNV_occ_prob_clim_ens_fut_ssp585 <- C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585 * WNV_occ_prob_clim_ens_fut_ssp585




#-------------------------------------------------------------------------------

# 4. Binarise post-processed continuous prediction data  -----------------------

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






# # For TBE and WNV, we calculate a new threshold based on multiplied cross-validated
# # predictions of virus and vector (based on known presences and absences)
# # Load needed objects
# source("scripts/00_functions.R") # Get the function for SDM evaluation
# 
# # Ixodes ricinus
# load("output_data/validation/TBE_validation.RData") # Load validation results
# load("output_data/models/TBE_SDMs.RData") # Data frame containing all presences and absences matched with occ. prob. of vector
# TBE_multiplied_ens_preds_cv <- m_ens_preds_cv * TBE_occ_env$I_ricinus # Multiply occurrence probabilities
# TBE_multiplied_ens_perf_cv <- evalSDM(TBE_occ_env$occ, TBE_multiplied_ens_preds_cv, weights = weights) # Calculate threshold
# TBE_thresh <- TBE_multiplied_ens_perf_cv$thresh # Extract threshold
# 
# # Culex pipiens
# load("output_data/validation/WNV_validation.RData")
# load("output_data/models/WNV_SDMs.RData")
# WNV_multiplied_ens_preds_cv <- m_ens_preds_cv * WNV_occ_env$C_pipiens
# WNV_multiplied_ens_perf_cv <- evalSDM(WNV_occ_env$occ, WNV_multiplied_ens_preds_cv, weights = weights)
# WNV_thresh <- WNV_multiplied_ens_perf_cv$thresh




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
TBE_occ_prob_clim_ens <- TBE_occ_prob_clim_ens * I_ricinus_occ_prob_clim_landuse_ens_bin
TBE_occ_prob_noclim_ens <- TBE_occ_prob_noclim_ens * I_ricinus_occ_prob_noclim_landuse_ens_bin
TBE_occ_prob_clim_nolanduse_ens <- TBE_occ_prob_clim_nolanduse_ens * I_ricinus_occ_prob_clim_nolanduse_ens_bin
TBE_occ_prob_noclim_nolanduse_ens <- TBE_occ_prob_noclim_nolanduse_ens * I_ricinus_occ_prob_noclim_nolanduse_ens_bin

# WNV based on predicted presence/absence of Culex pipiens
WNV_occ_prob_clim_ens <- WNV_occ_prob_clim_ens * C_pipiens_occ_prob_clim_landuse_ens_bin
WNV_occ_prob_noclim_ens <- WNV_occ_prob_noclim_ens * C_pipiens_occ_prob_noclim_landuse_ens_bin
WNV_occ_prob_clim_nolanduse_ens <- WNV_occ_prob_clim_nolanduse_ens * C_pipiens_occ_prob_clim_nolanduse_ens_bin
WNV_occ_prob_noclim_nolanduse_ens <- WNV_occ_prob_noclim_nolanduse_ens * C_pipiens_occ_prob_noclim_nolanduse_ens_bin



# Apply the threshold to all predictions cells of prediction rasters for
# each virus and the different scenarios
# TBE
TBE_occ_prob_clim_ens_bin <- TBE_occ_prob_clim_ens >= TBE_thresh # Apply corresponding threshold
TBE_occ_prob_clim_ens_bin <- TBE_occ_prob_clim_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

TBE_occ_prob_noclim_ens_bin <- TBE_occ_prob_noclim_ens >= TBE_thresh
TBE_occ_prob_noclim_ens_bin <- TBE_occ_prob_noclim_ens_bin * 1

TBE_occ_prob_clim_nolanduse_ens_bin <- TBE_occ_prob_clim_nolanduse_ens >= TBE_thresh
TBE_occ_prob_clim_nolanduse_ens_bin <- TBE_occ_prob_clim_nolanduse_ens_bin * 1

TBE_occ_prob_noclim_nolanduse_ens_bin <- TBE_occ_prob_noclim_nolanduse_ens >= TBE_thresh
TBE_occ_prob_noclim_nolanduse_ens_bin <- TBE_occ_prob_noclim_nolanduse_ens_bin * 1

# WNV
WNV_occ_prob_clim_ens_bin <- WNV_occ_prob_clim_ens >= WNV_thresh # Apply corresponding threshold
WNV_occ_prob_clim_ens_bin <- WNV_occ_prob_clim_ens_bin * 1 # Convert TRUE/FALSE to 1/0 values

WNV_occ_prob_noclim_ens_bin <- WNV_occ_prob_noclim_ens >= WNV_thresh
WNV_occ_prob_noclim_ens_bin <- WNV_occ_prob_noclim_ens_bin * 1

WNV_occ_prob_clim_nolanduse_ens_bin <- WNV_occ_prob_clim_nolanduse_ens >= WNV_thresh
WNV_occ_prob_clim_nolanduse_ens_bin <- WNV_occ_prob_clim_nolanduse_ens_bin * 1

WNV_occ_prob_noclim_nolanduse_ens_bin <- WNV_occ_prob_noclim_nolanduse_ens >= WNV_thresh
WNV_occ_prob_noclim_nolanduse_ens_bin <- WNV_occ_prob_noclim_nolanduse_ens_bin * 1






# (b) Predictions based on future data -----------------------------------------

# Apply the threshold to all predictions cells of prediction rasters for
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
TBE_occ_prob_clim_ens_fut_ssp126 <- TBE_occ_prob_clim_ens_fut_ssp126 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126
TBE_occ_prob_clim_ens_fut_ssp370 <- TBE_occ_prob_clim_ens_fut_ssp370 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370
TBE_occ_prob_clim_ens_fut_ssp585 <- TBE_occ_prob_clim_ens_fut_ssp585 * I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585

# WNV based on predicted presence/absence of Culex pipiens
WNV_occ_prob_clim_ens_fut_ssp126 <- WNV_occ_prob_clim_ens_fut_ssp126 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126
WNV_occ_prob_clim_ens_fut_ssp370 <- WNV_occ_prob_clim_ens_fut_ssp370 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370
WNV_occ_prob_clim_ens_fut_ssp585 <- WNV_occ_prob_clim_ens_fut_ssp585 * C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585





# Apply the threshold to all predictions cells of prediction rasters for
# each virus and the different scenarios
# TBE
TBE_occ_prob_clim_ens_bin_fut_ssp126 <- TBE_occ_prob_clim_ens_fut_ssp126 >= TBE_thresh
TBE_occ_prob_clim_ens_bin_fut_ssp126 <- TBE_occ_prob_clim_ens_bin_fut_ssp126 * 1

TBE_occ_prob_clim_ens_bin_fut_ssp370 <- TBE_occ_prob_clim_ens_fut_ssp370 >= TBE_thresh
TBE_occ_prob_clim_ens_bin_fut_ssp370 <- TBE_occ_prob_clim_ens_bin_fut_ssp370 * 1

TBE_occ_prob_clim_ens_bin_fut_ssp585 <- TBE_occ_prob_clim_ens_fut_ssp585 >= TBE_thresh
TBE_occ_prob_clim_ens_bin_fut_ssp585 <- TBE_occ_prob_clim_ens_bin_fut_ssp585 * 1

# WNV
WNV_occ_prob_clim_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_ens_fut_ssp126 >= WNV_thresh
WNV_occ_prob_clim_ens_bin_fut_ssp126 <- WNV_occ_prob_clim_ens_bin_fut_ssp126 * 1

WNV_occ_prob_clim_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_ens_fut_ssp370 >= WNV_thresh
WNV_occ_prob_clim_ens_bin_fut_ssp370 <- WNV_occ_prob_clim_ens_bin_fut_ssp370 * 1

WNV_occ_prob_clim_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_ens_fut_ssp585 >= WNV_thresh
WNV_occ_prob_clim_ens_bin_fut_ssp585 <- WNV_occ_prob_clim_ens_bin_fut_ssp585 * 1






#-------------------------------------------------------------------------------

# 5. Save post-processed prediction rasters ------------------------------------
# Saving continuous and binary predictions of different past and future scenarios

# (a) Predictions based on historical data -------------------------------------

# Continuous and binary predictions of different scenarios
# Ixodes ricinus
writeRaster(I_ricinus_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/I_ricinus_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# Culex pipiens
writeRaster(C_pipiens_occ_prob_clim_landuse_ens, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_landuse_ens, "output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_landuse_ens_bin, "output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_landuse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/C_pipiens_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# TBE
writeRaster(TBE_occ_prob_clim_ens, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_ens, "output_data/results/postprocessed_predictions/TBE_preds_noclim_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/TBE_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/TBE_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(TBE_occ_prob_clim_ens_bin, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_ens_bin, "output_data/results/postprocessed_predictions/TBE_preds_noclim_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/TBE_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/TBE_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)

# WNV
writeRaster(WNV_occ_prob_clim_ens, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_ens, "output_data/results/postprocessed_predictions/WNV_preds_noclim_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_nolanduse_ens, "output_data/results/postprocessed_predictions/WNV_preds_clim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_nolanduse_ens, "output_data/results/postprocessed_predictions/WNV_preds_noclim_nolanduse_ens_1970_2019.tif", overwrite = TRUE)

writeRaster(WNV_occ_prob_clim_ens_bin, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_ens_bin, "output_data/results/postprocessed_predictions/WNV_preds_noclim_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/WNV_preds_clim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_noclim_nolanduse_ens_bin, "output_data/results/postprocessed_predictions/WNV_preds_noclim_nolanduse_ens_bin_1970_2019.tif", overwrite = TRUE)




# (b) Predictions based on future data -----------------------------------------

# Ixodes ricinus
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_2030_2070_ssp585.tif", overwrite = TRUE)

writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(I_ricinus_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/I_ricinus_preds_clim_landuse_ens_bin_2030_2070_ssp585.tif", overwrite = TRUE)

# Culex pipiens
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp126, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp370, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_fut_ssp585, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_2030_2070_ssp585.tif", overwrite = TRUE)

writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(C_pipiens_occ_prob_clim_landuse_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/C_pipiens_preds_clim_landuse_ens_bin_2030_2070_ssp585.tif", overwrite = TRUE)

# TBE
writeRaster(TBE_occ_prob_clim_ens_fut_ssp126, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_ens_fut_ssp370, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_ens_fut_ssp585, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_2030_2070_ssp585.tif", overwrite = TRUE)

writeRaster(TBE_occ_prob_clim_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(TBE_occ_prob_clim_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/TBE_preds_clim_ens_bin_2030_2070_ssp585.tif", overwrite = TRUE)

# WNV
writeRaster(WNV_occ_prob_clim_ens_fut_ssp126, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_ens_fut_ssp370, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_ens_fut_ssp585, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_2030_2070_ssp585.tif", overwrite = TRUE)

writeRaster(WNV_occ_prob_clim_ens_bin_fut_ssp126, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp126.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_ens_bin_fut_ssp370, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp370.tif", overwrite = TRUE)
writeRaster(WNV_occ_prob_clim_ens_bin_fut_ssp585, "output_data/results/postprocessed_predictions/WNV_preds_clim_ens_bin_2030_2070_ssp585.tif", overwrite = TRUE)



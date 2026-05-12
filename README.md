# Past and future phenology changes of zoonotic vector-borne diseases under climate and land-use change

Valén Holle<sup>1</sup>, Raphaëlle Klitting<sup>2</sup>, Nadja Kabisch<sup>3</sup>, Damaris Zurell<sup>1</sup>

1. University of Potsdam, Institute of Biochemistry and Biology, Potsdam, Germany
2. Aix-Marseille University, Emerging Viruses Unit, Marseille, France
3. Leibniz University Hannover, Institute of Earth System Sciences, Hannover, Germany



### ABSTRACT:
Environmental changes are reshaping the distribution and seasonal dynamics of vector-borne diseases, with important implications for public health. Tick-borne encephalitis virus (TBEV) and West Nile virus (WNV) cause growing concern in Europe, with rising case numbers and ever-expanding circulation areas. The transmission risk of TBEV and WNV follows characteristic seasonal patterns, driven largely by weather-dependent activity of their arthropod vectors. The relative roles of climate and land-use change on the seasonal dynamics and spread of these diseases and their vectors remain, however, poorly quantified. Here, we assess the spread and phenology of TBEV and WNV in response to historical and future climate and land-use changes across Europe. We developed spatiotemporal species distribution models (SDMs) for the viruses and their primary vector species, generating monthly environmental suitability predictions from the 1970s to 2050s. Virus models incorporated vector suitability as a nested predictor to capture the dependence of virus occurrence on vector presence. To disentangle drivers of observed changes, we applied counterfactual historical simulations, attributing shifts in seasonal transmission risk to climate or land-use changes.  Historical attribution results show that land-use changes mainly affected vector suitability, whereas climatic changes drove shifts in seasonal transmission risk. Transmission risk is projected to rise continent-wide for both TBEV and WNV over the coming decades. Further, TBEV is projected to undergo pronounced phenological shifts, with a dominant spring peak and a delayed autumn peak extending into October. Prolonged seasonal transmission windows are projected to create hotspots that both intensify and expand across large regions. Taken together, our findings underscore the need for coordinated transnational efforts to manage the projected health burden of TBEV and WNV across Europe, and support upstream prevention by providing climate-informed guidance on intervention timing and spatial prioritisation.

Keywords: Climate change, Detection and attribution, Europe, Land-use change, Disease phenology, Species distribution models (SDMs), Tick-borne encephalitis virus (TBEV), Seasonality, Vector-borne diseases, West Nile virus (WNV)

Funding: This work was supported by the European Union through the project Zoonosis Emergence across Degraded and Restored Forest Ecosystems (project no. 101135094).

This repository contains the R scripts needed to reproduce all results and plots.

---------------------------------------------------------------
**Workflow**
---------------------------------------------------------------

All data preparation and modelling steps are detailed in the [ODMAP protocol].

### 00 - Setup
scripts [folder setup](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/00_folder_setup.R) and [functions](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/00_functions.R)

We create all necessary directories for the project and list the needed functions.


### 01 - Species data preparation
scripts [01a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/01a_Ixodes_ricinus_species_data_prep.R), [01b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/01b_Culex_pipiens_species_data_prep.R)

We process the raw species occurrence data for *Ixodes ricinus* and *Culex pipiens* ([GBIF](https://www.gbif.org), [VectorMap](https://experience.arcgis.com/experience/5f95c3edfbea4634b8347fec0bd1dcd6)) by removing records with erroneous timestamps or coordinates and excluding duplicate records within the same calendar month of a given year. This ensures a clean dataset with a target temporal resolution of one month.


### 02 - Environmental data preparation
scripts [02a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/02a_env_data_prep_ISIMIP3a.R), [2b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/02b_env_data_prep_ISIMIP3b.R)

We process [ISIMIP3 (Inter-Sectoral Impact Model Intercomparison Project Phase 3)](https://www.isimip.org/protocol/3/) climate and land-use data, including historical data (ISIMIP3a) and future projections (ISIMIP3b). The climate data are provided as daily outputs at a 0.5° resolution. To match our target temporal resolution, daily temperature and humidity outputs are aggregated into monthly mean values, and daily precipitation outputs into monthly totals. Yearly historical land-use data at a 0.5° resolution are summarised into eight distinct land-use categories per year. Future land-use data are provided by [LUH2 (Land Use Harmonization 2)](https://luh.umd.edu/data.shtml) at 0.25° resolution; these rasters are aggregated by a factor of 2 to achieve a 0.5° spatial resolution and similarly summarised into the eight distinct land-use categories for each year.


### 03 - Species background data generation
scripts [03a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/03a_Ixodes_ricinus_background_data_prep.R), [03b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/03b_Culex_pipiens_background_data_prep.R)

We generate background data for our vector species *Ixodes ricinus* and *Culex pipiens* by randomly selecting locations within a specified buffer distance around the presence points, aiming for a presence-to-background ratio of 1:10, while excluding cells containing the actual presence records. This procedure is conducted separately based on occurrences within the same month of a given year, resulting in temporally matched background data. To reduce spatial autocorrelation, we thin both the monthly presence and background data of the species using a 50 km threshold. Finally, we match the species data - comprising both presence and background data - with the month- and year-specific climate predictors, as well as year-specific land-use predictors.


### 04 - Species model fitting
scripts [04a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/04a_Ixodes_ricinus_models.R), [04b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/04b_Culex_pipiens_models.R)

For the species *Ixodes ricinus* and *Culex pipiens*, we identify the most important and weakly correlated predictor variables to include in the SDM construction. Spatiotemporal models are built using four different algorithms: Generalised Linear Model (GLM), Generalised Additive Model (GAM), Random Forest (RF), and Boosted Regression Tree (BRT). 


### 05 - Species model validation
scripts [05a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/05a_Ixodes_ricinus_validation.R), [05b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/05b_Culex_pipiens_validation.R)

We evaluate the model performance of all algorithms for *Ixodes ricinus* and *Culex pipiens* using a 5-fold cross-validation approach, focusing on performance metrics such as AUC and the Boyce index. To further assess the ensemble model, we calculate the average of the continuous cross-validated predictions. Additionally, we extract monthly performance measures and generate plots of response curves for each predictor variable.  


### 06 - Species model predictions
scripts [06a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/06a_Ixodes_ricinus_prediction.R), [06b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/06b_Culex_pipiens_prediction.R)

We generate continuous ensemble predictions for historical time periods (1970-2019) using both observed and counterfactual environmental data, allowing us to attribute changes in predictions to climate and land-use factors. For future time periods (2030-2059), we incorporate different socio-economic forcing scenarios (ssp126, ssp370, ssp585) and, for climate, consider projections from five distinct climate models.


### 07 - Virus data preparation
scripts [07a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/07a_TBE_infection_data_prep.R), [07b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/07b_WNV_infection_data_prep.R)

We process locally acquired, confirmed human TBEV and WNV cases in Europe, reported at the NUTS3 level (provided by [TESSy/ECDC](https://atlas.ecdc.europa.eu/public/index.aspx)) to generate spatially explicit infection data. This is done by rasterising the European NUTS3 municipalities to the target spatial resolution of 0.5° and extracting the central coordinates of municipalities with observed infections. For NUTS3 municipalities that consist of only one cell after rasterisation, these coordinates are used as the final location for the infection data point. In cases where NUTS3 municipalities are too small to be represented by a 0.5° cell in the rasterisation process but had reported infections, their central coordinates are likewise considered as the infection data point. Infection records stemming from large NUTS3 municipalities spanning more than one cell are excluded to minimise the spatial uncertainty of infection locations.


### 08 - Virus absence data generation
scripts [08a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/08a_TBE_absence_data_prep.R), [08b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/08b_WNV_absence_data_prep.R)

We generate absence data for the TBEV and WNV occurrences by drawing points from all cells within NUTS3 municipalities in EU/EEA countries with a mandatory surveillance system, selecting only cells where no infections had been reported to the ECDC. Because the set of EU/EEA countries with mandatory reporting changed over time, we adjust the list of eligible countries for each year based on information from the [corresponding Annual Epidemiological Reports](https://www.ecdc.europa.eu/en/publications-data/monitoring/all-annual-epidemiological-reports). To ensure consistency with our temporal resolution, absence points are generated separately from occurrences within the same month of a given year. To minimise spatial autocorrelation, a 50 km thinning threshold was applied to both the monthly infection presence data and the corresponding absence data. Because pathogen distribution depends on the presence of their vector species, we include, in addition to climate variables, the predicted suitability of the main vector species as a predictor variable, following a nested modelling approach. 


### 09 - Virus model fitting
scripts [09a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/09a_TBE_models.R), [09b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/09b_WNV_models.R)

For the target viruses, we create a balanced dataset of thinned presence and absence points. We then identify the most important and weakly correlated predictor variables to include in the SDM construction. Spatiotemporal models are built using four different algorithms: Generalised Linear Model (GLM), Generalised Additive Model (GAM), Random Forest (RF), and Boosted Regression Tree (BRT). 


### 10 - Virus model validation
scripts [10a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/10a_TBE_validation.R), [10b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/10b_WNV_validation.R), [10c](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/10c_part_resp_curve_plot_vector_virus.R)

We evaluate the model performance of all algorithms for TBEV and WNV using a 5-fold cross-validation approach, focusing on performance metrics such as AUC and the Boyce index. To further assess the ensemble model, we calculate the average of the continuous cross-validated predictions. Monthly performance measures are extracted and response curves are generated for each predictor variable. Finally, we create a combined plot showing the partial response curves for both the target vectors and the associated viruses.


### 11 - Virus model predictions
scripts [11a](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/11a_TBE_prediction.R), [11b](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/11b_WNV_prediction.R)

We generate continuous ensemble predictions for historical time periods (1970-2019) using both observed and counterfactual environmental data, allowing us to attribute changes in predictions. For future time periods (2030-2059), we incorporate different socio-economic forcing scenarios (ssp126, ssp370, ssp585) based five distinct climate models. Although land-use variables were not directly included in the virus models, the effects of land-use change are implicitly accounted for through the incorporation of the corresponding suitability predictions of the main vector species.


### 12 - Vector and virus prediction postprocessing
script [12](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/12_prediction_postprocessing.R)

We post-process the ensemble prediction rasters for the viruses and their primary vector species in three main steps. First, future layers from the same socio-economic scenario are averaged across the five different climate models. Second, vector prediction rasters are masked to include only cells within EU/EEA countries to ensure comparability between vector and virus outputs. Third, for the virus layers, cells are set to 0 whenever the corresponding vector species is not predicted to be present, ensuring for ecological realism.


### 13 - Decadal trends in distribution
script [13](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/13_distribution_trends_vector_virus.R)

We visualise the distribution trends of the vectors and viruses for specific target months - May and October for *Ixodes ricinus* and TBEV, and July for *Culex pipiens* and WNV - across three different target decades: 1970s, 2010, 2050s. For the future target decade (2050s), we look at predictions the were derived from the three different socio-economic scenarios. Only grid cells with vector or virus suitability values indicating at least one predicted presence within the selected month of a given decade are displayed. 


### 14 - Decadal trends in phenology intensity
script [14](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/14_decadal_trends_vector_virus.R)

We calculate and visualise decadal trends in the timing of mean and peak vector and virus suitability throughout the year from historical to future time periods. In doing so, we compare historical phenology trends based on factual and counterfactual scenarios, allowing us to disentangle the relative impacts of climate and land-use changes. We also compare future phenology trends across three different socio-economic scenarios (ssp126, ssp37, ssp585), representing different potential environmental trajectories.


### 15 - Decadal trends in phenology intensity per Köppen-Geiger climate region
script [15](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/15_climateregions_trends_vector_virus.R)

We calculate and visualise the decadal trends in the timing of mean and peak vector and virus suitability throughout the year, separately for the major Köppen-Geiger climate classes. Historical predictions are based on factual climate and land-use changes, while future predictions consider the three different socio-economic scenarios. As climate zones in Europe have shifted over the past decades, we use the Köppen-Geiger classification map corresponding to each target decade (1970s, 2010s, 2050s).

### 16 - Decadal trends in phenology duration
script [16](https://github.com/UP-macroecology/Holle_ZOE_DiseasePhenology_2025/blob/main/scripts/16_duration_trends_vector_virus.R)

We calculate and visualise temporal trends in the duration of vector activity and potential virus transmission periods throughout the year. Using the ensemble predictions based on the factual historical climate and land-use changes, as well as on the three future socio-economic scenarios, we first determine the number of months with predicted presence in each cell. Next, we compute the mean duration per decade for each cell. To assess temporal historical and future trends per cell, we fit two separate linear models for each cell: one based on the decadal mean durations of the historical decades (1970s - 2010s), and one based on those of the future decades (2010s - 2050s). The slopes of these models serve as an indicator of change. 



---------------------------------------------------------------
**Required folder structure**
---------------------------------------------------------------

```

scripts

input_data
├── environmental_data
    ├── ISIMIP3a
        ├── Climate
            ├── processed_data
            ├── raw_data
        ├── CounterClim
            ├── processed_data
            ├── raw_data
        ├── CounterLandUse
            ├── processed_data
            ├── raw_data
        ├── LandUse
            ├── processed_data
            ├── raw_data
    ├── ISIMIP3b
        ├── Climate
            ├── ssp126
                ├── processed_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
                ├── raw_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
            ├── ssp370
                ├── processed_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
                ├── raw_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
            ├── ssp585
                ├── processed_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
                ├── raw_data
                    ├── gfdl-esm4
                    ├── ipsl-cm6a-lr
                    ├── mpi-esm1-2-hr
                    ├── mri-esm2-0
                    ├── ukesm1-0-ll
        ├── LandUse
            ├── ssp126
                ├── processed_data
                ├── raw_data
            ├── ssp370
                ├── processed_data
                ├── raw_data
            ├── ssp585
                ├── processed_data
                ├── raw_data
├── raw_infection_data
├── raw_species_data
├── spatial_data
    ├──climate_regions
    

output_data
├── data
├── models
├── validation
├── results
    ├── preprocessed_predictions
        ├── Culex_pipiens
        ├── Ixodes_ricinus
        ├── WNV
        ├── TBE
    ├── postprocessed_predictions
        ├── Culex_pipiens
        ├── Ixodes_ricinus
        ├── WNV
        ├── TBE
    ├── decadal_trends
    ├── duration_trends
    ├── climateregions_trends
    ├── distribution_trends
├── plots
    ├── maps
    ├── presence_background
    ├── response_curves
    ├── distribution_trends
    ├── decadal_trends
    ├── duration_trends
    ├── climateregions_trends
    
```


---------------------------------------------------------------
**Required data**
---------------------------------------------------------------

* Vector species occurrence data are available from [GBIF](https://www.gbif.org) and [VectorMap](https://experience.arcgis.com/experience/5f95c3edfbea4634b8347fec0bd1dcd6)
* Limited human TBE and WNV case infection data are available from [TESSy/ECDC](https://atlas.ecdc.europa.eu/public/index.aspx); comprehensive data are available upon request (subject to non-redistribution conditions)
* European administrative boundary data (shapefile) is available [here](https://hub.arcgis.com/datasets/bdcb40c0b6124f6d99f10b9b23647712/explore)
* Historical and future climate data, as well as historical land-use data, are available from [ISIMIP3](https://data.isimip.org/search/)
* Future land-use data are available from [LUH2](https://luh.umd.edu/data.shtml)
* Köppen-Geiger climate classification maps are available from [GLOH2O](https://www.gloh2o.org/koppen/)



---------------------------------------------------------------
**Operating system info**
---------------------------------------------------------------

* R version 4.3.1 (2023-06-16 ucrt)
* Platform: x86_64-w64-mingw32/x64 (64-bit)
* Running under: Windows 11 x64 (build 26100)

* Attached packages:
[1] CoordinateCleaner_3.0.1 [2] corrplot_0.92 [3] countrycode_1.6.0 [4] dplyr_1.1.3 [5] dismo_1.3-14 [6] gbm_2.1.8.1 [7] ggh4x_0.3.0 [8] giscoR_0.6.0 [9] ggnewscale_0.5.1 [10] ggplot2_4.0.0 [11] lubridate_1.9.3 [12] maps_3.4.1 [13] mgcv_1.8-42 [14] PresenceAbsence_1.1.11 [15] randomForest_4.7-1.1 [16] readr_2.1.4 [17] sf_1.0-16 [18] sfheaders_0.4.3 [19] stringr_1.5.0 [20] terra_1.7-55 [21] tidyr_1.3.0 [22] tibble_3.2.1 [23] tidyterra_0.6.1 [24] tidyverse_2.0.0 [25] viridis_0.6.4


# Working title: Historical and future disease phenologies of TBE and WNV in Europe

### RESEARCH AIM: 
The analysis is part of the EU-funded project [ZOE (Zoonoses Emergence across Degraded and Restored Forest Ecosystems)](https://www.zoe-project.eu), which investigates the relationship between ecosystem degradation, biodiversity loss, and the associated risk of zoonotic disease emergence. As contributors to Work Package 5 of the project, our goal is to analyse how distinct seasonal emergence patterns of arthropod species can significantly impact disease transmission risks in Europe. Using a spatiotemporal modelling framework, we aim to understand how climate warming and land use changes, as well as interannual variation, can alter disease occurrence in space and time, the timing of peak infection risk and the duration of transmission season by affecting vector and virus distribution. Understanding these trends is crucial for effectively managing disease risks and implementing appropriate public health measures. 

This repository contains the R scripts needed to reproduce all results and plots.

## Workflow
### 00 - Data setup
We list the needed functions.

### 01 - Species and infection data preparation
We process raw species (*Ixodes ricinus*, *Culex pipiens*) and infection (TBEV, WNV) data to obtain a target spatial resolution of 50 km and a target temporal resolution of one month. 

### 02 - Environmental data preparation
We process [ISIMIP3 (Inter-Sectoral Impact Model Intercomparison Project Phase 3)](https://www.isimip.org/protocol/3/) climate and land-use data, including historical simulations (ISIMIP3a) and future projections (ISIMIP3b). The climate data is provided as daily outputs at a 50 km resolution, matching the spatial resolution of the species and infection data. To additionally align with our target temporal resolution, we aggregate the daily temperature and humidity outputs into monthly mean values and daily precipitation outputs into monthly totals. Additionally, yearly land-use data at a spatial resolution of 50 km is aggregated into eight distinct land-use categories for each year.

### 03 - Species background data generation
We generate background data for our vector species *Ixodes ricinus* and *Culex pipiens* by randomly selecting locations within a specified buffer distance from the presence points using a presence-background ratio of 1:10, excluding cells containing the actual presence locations. This process is conducted separately based on occurrences within the same month of a given year, resulting in temporally matched background data. To avoid spatial autocorrelation, we thin both the presence and background data of the species using a 100 km threshold. Finally, we match the species data - comprising both presence and background data - with the month- and year-specific climate predictors, as well as year-specific land-use predictors.

### 04 - Species model fitting
For the species *Ixodes ricinus* and *Culex pipiens*, we identify the most important and weakly correlated predictor variables to include in model construction. Models are built using four different algorithms: Generalised Linear Model (GLM), Generalised Additive Model (GAM), Random Forest (RF), and Boosted Regression Tree (BRT). 

### 05 - Species model validation
We evaluate the model performance of all algorithms for *Ixodes ricinus* and *Culex pipiens* using a 5-fold cross-validation approach, focusing on performance metrics such as AUC and the Boyce index. To further assess the ensemble model, we calculate the average of the continuous cross-validated predictions. Additionally, we extract monthly performance measures and generate plots of response curves for each predictor variable.  

### 06 - Species model predictions
We generate continuous ensemble predictions for historical time periods under various environmental scenarios (observed climate and land use change, observed land use change and detrended climate, observed climate change and no land use). For future time periods, we consider different environmental scenarios (projected climate and land use change, projected climate change and steady land use, projected climate change and no land use) and forcings (ssp126, ssp370, ssp585).

### 07 - Infection background data generation
We generate background data for the human-case infection data of TBEV and WNV by sampling cell locations where the main vector species were present, but the disease itself was absent. These background points are created separately based on disease and vector occurrences within the same month of a given year to ensure alignment with our temporal resolution. After thinning the infection data using a 100 km threshold, we match the time-specific infection data with the corresponding environmental data. This time, we include the predicted monthly occurrence probability of the main vector species as a predictor variable.

### XX - Decadal occurrence probability trends

### XX - Latitudinal occurrence probability trends



# ZOE project 
# Disease phenology analysis of Culex pipiens in Europe (primary transmitter of WNV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                          05b. Model validation                         #
# ---------------------------------------------------------------------- #

# Load needed packages
library(mgcv)
library(maxnet)
library(randomForest)
library(gbm)
library(dismo)
library(ggplot2)
library(PresenceAbsence)
library(tidyverse)

# Load needed objects
source("scripts/00_functions.R") # Get the function for SDM evaluation,
# Boyce index with smoothing methods (Liu et al. (2024)), and predict function

# Read in SDM models
load(paste0("output_data/models/C_pipiens_SDMs.RData"))



#-------------------------------------------------------------------------------

# 1. 5-fold cross-validation ---------------------------------------------------

# Part the presence-background dataset into 5 folds
kfolds <- 5
ks <- dismo::kfold(seq_len(nrow(C_pipiens_occ_env)), k = kfolds)


# GLM
print("GLM")
m_glm_preds_cv <- rep(NA, nrow(C_pipiens_occ_env))
for(i in seq_len(kfolds)) { # Start of the loop over k folds
  cv_train <- C_pipiens_occ_env[ks != i,]
  cv_test <- C_pipiens_occ_env[ks == i,]
  cv_weights <- weights[ks != i]
  
  cv_glm <- update(m_glm, data = cv_train, weights = cv_weights)
  m_glm_preds_cv[ks == i] <- predict(cv_glm, cv_test, type = "response")
}

# Calculate performance measures
m_glm_perf_cv <- evalSDM(C_pipiens_occ_env$occ, m_glm_preds_cv)

# Calculate Boyce index with smoothing methods
m_glm_preds_cv_presences <- m_glm_preds_cv[presences] # Just retain the predictions of the presences based on indices
absences <- which(C_pipiens_occ_env$occ == 0) # Extract the position index of presences
m_glm_preds_cv_absences <- m_glm_preds_cv[absences] # Just retain the predictions of the absences based on indices
m_glm_boyce_cv <- sfbi(m_glm_preds_cv_presences, m_glm_preds_cv_absences, ktry = 10) # Apply sfbi function

boyce_index_m_glm <- m_glm_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers

# Add the Boyce index to the performance metrics data frame
m_glm_perf_cv$Boyce <- boyce_index_m_glm




# GAM
print("GAM")
m_gam_preds_cv <- rep(NA, nrow(C_pipiens_occ_env))
for(i in seq_len(kfolds)) { # Start of the loop over k folds
  cv_train <- C_pipiens_occ_env[ks != i,]
  cv_test <- C_pipiens_occ_env[ks == i,]
  cv_weights <- weights[ks != i]
  
  cv_gam <- update(m_gam, data = cv_train, weights = cv_weights)
  m_gam_preds_cv[ks == i] <- predict(cv_gam, cv_test, type = "response")
}

# Calculate performance measures
m_gam_perf_cv <- evalSDM(C_pipiens_occ_env$occ, m_gam_preds_cv)

# Calculate Boyce index with smoothing methods
m_gam_preds_cv_presences <- m_gam_preds_cv[presences] # Just retain the predictions of the presences based on indices
absences <- which(C_pipiens_occ_env$occ == 0) # Extract the position index of presences
m_gam_preds_cv_absences <- m_gam_preds_cv[absences] # Just retain the predictions of the absences based on indices
m_gam_boyce_cv <- sfbi(m_gam_preds_cv_presences, m_gam_preds_cv_absences, ktry = 10) # Apply sfbi function

boyce_index_m_gam <- m_gam_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers

# Add the Boyce index to the performance metrics data frame
m_gam_perf_cv$Boyce <- boyce_index_m_gam




# RF
print("RF")
m_rf_preds_matrix_cv <- matrix(nrow=nrow(C_pipiens_occ_env), ncol=ratio_presence_background)
for(i in seq_len(kfolds)){ # Start of the loop over k folds
  print(i)
  hold_in <- which(ks != i)
  hold_out <- which(ks == i)
  
  for (y in 1:ratio_presence_background) { # Start the loop over 10 models
    print(y)
    data_index <- c(presences, which(C_pipiens_occ_env$abs_index == i))
    
    cv_train <- C_pipiens_occ_env[hold_in[hold_in %in% data_index],]
    cv_test <- C_pipiens_occ_env[hold_out,]
    
    cv_rf <- update(m_rf[[y]], data = cv_train)
    m_rf_preds_matrix_cv[hold_out, y] <- predict(cv_rf, cv_test, type = "response")
  }
}

m_rf_preds_cv <- rowMeans(m_rf_preds_matrix_cv)

# Calculate performance measures
m_rf_perf_cv <- evalSDM(C_pipiens_occ_env$occ, m_rf_preds_cv)

# Calculate Boyce index with smoothing methods
m_rf_preds_cv_presences <- m_rf_preds_cv[presences] # Just retain the predictions of the presences based on indices
absences <- which(C_pipiens_occ_env$occ == 0) # Extract the position index of presences
m_rf_preds_cv_absences <- m_rf_preds_cv[absences] # Just retain the predictions of the absences based on indices
m_rf_boyce_cv <- sfbi(m_rf_preds_cv_presences, m_rf_preds_cv_absences, ktry = 10) # Apply sfbi function

boyce_index_m_rf <- m_rf_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers

# Add the Boyce index to the performance metrics data frame
m_rf_perf_cv$Boyce <- boyce_index_m_rf





# BRT
print("BRT")
m_brt_preds_matrix_cv <- matrix(nrow=nrow(C_pipiens_occ_env), ncol=ratio_presence_background)
for(i in seq_len(kfolds)){ # Start of the loop over k folds
  print(i)
  hold_in <- which(ks != i)
  hold_out <- which(ks == i)
  
  for (y in 1:ratio_presence_background) { # Start the loop over 10 models
    print(y)
    data_index <- c(presences, which(C_pipiens_occ_env$abs_index == i))
    
    cv_train <- C_pipiens_occ_env[hold_in[hold_in %in% data_index],]
    names(cv_train)[names(cv_train)=='occ'] <- m_brt[[y]]$response.name
    cv_test <- C_pipiens_occ_env[hold_out,]
    
    cv_brt <- gbm::gbm(m_brt[[y]]$call, 'bernoulli', data = cv_train[,c( m_brt[[y]]$response.name, my_preds)],
                       n.trees=m_brt[[y]]$gbm.call$best.trees,
                       shrinkage=m_brt[[y]]$gbm.call$learning.rate,
                       bag.fraction=m_brt[[y]]$gbm.call$bag.fraction,
                       interaction.depth=m_brt[[y]]$gbm.call$tree.complexity)
    m_brt_preds_matrix_cv[hold_out,y] <- predict(cv_brt, cv_test, type='response', n.trees=m_brt[[y]]$gbm.call$best.trees)
  }
}

m_brt_preds_cv <- rowMeans(m_brt_preds_matrix_cv)

# Calculate performance measures
m_brt_perf_cv <- evalSDM(C_pipiens_occ_env$occ, m_brt_preds_cv)

# Calculate Boyce index with smoothing methods
m_brt_preds_cv_presences <- m_brt_preds_cv[presences] # Just retain the predictions of the presences based on indices
absences <- which(C_pipiens_occ_env$occ == 0) # Extract the position index of presences
m_brt_preds_cv_absences <- m_brt_preds_cv[absences] # Just retain the predictions of the absences based on indices
m_brt_boyce_cv <- sfbi(m_brt_preds_cv_presences, m_brt_preds_cv_absences, ktry = 10) # Apply sfbi function

boyce_index_m_brt <- m_brt_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers

# Add the Boyce index to the performance metrics data frame
m_brt_perf_cv$Boyce <- boyce_index_m_brt



#-------------------------------------------------------------------------------

# 2. Make ensemble prediction --------------------------------------------------

# Make ensemble
print("Ensemble")
m_ens_preds_cv <- rowMeans(data.frame(m_glm_preds_cv, m_gam_preds_cv, m_brt_preds_cv, m_rf_preds_cv))
m_ens_perf_cv <- evalSDM(C_pipiens_occ_env$occ, m_ens_preds_cv)
m_ens_preds_cv_presences <- m_ens_preds_cv[presences]
absences <- which(C_pipiens_occ_env$occ == 0)
m_ens_preds_cv_absences <- m_ens_preds_cv[absences]
m_ens_boyce_cv <- sfbi(m_ens_preds_cv_presences, m_ens_preds_cv_absences, ktry = 10)
boyce_index_m_ens <- m_ens_boyce_cv[6]
m_ens_perf_cv$Boyce <- boyce_index_m_ens

# Compare cross-validated model performance across algorithms
comp_perf <- rbind(glm = m_glm_perf_cv, gam = m_gam_perf_cv, rf = m_rf_perf_cv, brt = m_brt_perf_cv, mean_prob = m_ens_perf_cv)


# Save cross-validation outcomes
save(m_glm_preds_cv, m_glm_perf_cv, m_gam_preds_cv, m_gam_perf_cv, m_brt_preds_cv, m_brt_perf_cv, 
     m_rf_preds_cv, m_rf_perf_cv, m_ens_preds_cv, m_ens_perf_cv, ks, comp_perf,
     file = "output_data/validation/C_pipiens_validation.RData")



#-------------------------------------------------------------------------------

# 3. Monthly model performances ------------------------------------------------

# Read in the needed data
load("output_data/models/C_pipiens_SDMs.RData") # The occurrence data frame
load("output_data/validation/C_pipiens_validation.RData") # The cross-validated ensemble predictions

# Create a vector containing the months of a year
month <- str_pad(1:12, width = 2, pad = "0")


for (m in month) { # Start of the loop over all months
  
  print(m)
  
  # Subset the occurrence data frame to only contain data of the respective month
  C_pipiens_occ_env_month <- subset(C_pipiens_occ_env, C_pipiens_occ_env$month == m)
  
  # Get the indices of the presences and absences of the respective month
  monthly_data_indices <- which(C_pipiens_occ_env$month == m)
  
  # Subset the cross-validated ensemble predictions by the extracted indices for that month
  monthly_data <- m_ens_preds_cv[monthly_data_indices]
  
  # Calculate the performance measures (AUC, TSS, Sensitivity, Specificity) for the respective month
  m_ens_perf_cv_month <- evalSDM(C_pipiens_occ_env_month$occ, monthly_data)
  
  # Calculate the Boyce index with smoothing methods for the respective month
  presences_month <- which(C_pipiens_occ_env_month$occ == 1)
  m_ens_preds_cv_presences_month <- monthly_data[presences_month]
  absences_month <- which(C_pipiens_occ_env_month$occ == 0)
  m_ens_preds_cv_absences_month <- monthly_data[absences_month] 
  m_ens_boyce_cv_month  <- sfbi(m_ens_preds_cv_presences_month, m_ens_preds_cv_absences_month, ktry = 10) 
  boyce_index_m_ens_month <- m_ens_boyce_cv_month[6]
  
  # Add the Boyce index to the performance metrics data frame
  m_ens_perf_cv_month$Boyce <- boyce_index_m_ens_month
  
  
  # Save the data frame with performance measures
  save(m_ens_perf_cv_month,  file = paste0("output_data/validation/C_pipiens_monthly_validation_",m,".RData"))
  
  
} # Close the loop over all months




#-------------------------------------------------------------------------------

# 4. Create response curves of ensemble predictions ----------------------------

# Read in the needed data
load("output_data/models/C_pipiens_SDMs.RData") # The different models and selected predictor variables


for (m in my_preds) { # Loop through the predictor variables
  
  print(m)
  
  # Name the other predictors in a vector
  my_preds_minus_m <- setdiff(my_preds, m)
  
  # Create an environmental dummy dataset (keeping the other predictors at their mean)
  dummy_data <- data.frame(seq(min(C_pipiens_occ_env[,m], na.rm = TRUE), max(C_pipiens_occ_env[,m], na.rm = TRUE), length = 100))
  for (pred in my_preds_minus_m) {
    dummy_data[[pred]] <- mean(C_pipiens_occ_env[[pred]], na.rm = TRUE)
  }
  
  names(dummy_data)[1] <- m
  
  # Generate predictions to dummy data for each model
  print("GLM")
  response_preds_glm <- predictSDM(m_glm, dummy_data)
  
  print("GAM") 
  response_preds_gam <- predictSDM(m_gam, dummy_data)
  
  print("RF")
  response_preds_rf_all <- matrix(nrow=nrow(dummy_data), ncol=10)
  
  for (y in 1:ratio_presence_background) { # Start of the loop over the 10 models
    print(y)
    response_preds_rf_all[,y] <- predictSDM(m_rf[[y]], dummy_data)
  }
  
  response_preds_rf <- rowMeans(response_preds_rf_all)
  
  print("BRT")
  response_preds_brt_all <- matrix(nrow=nrow(dummy_data), ncol=10)
  
  for (y in 1:ratio_presence_background) { # Start of the loop over the 10 models
    print(y)
    response_preds_brt_all[,y] <- predictSDM(m_brt[[y]], dummy_data)
  }
  response_preds_brt <- rowMeans(response_preds_brt_all)
  
  # Combine the cross-validated predictions into an ensemble and calculate the sd
  all_preds <- cbind(response_preds_glm, response_preds_gam, response_preds_rf, response_preds_brt)
  ensemble_sd <- apply(all_preds, 1, sd)
  ensemble_preds <- rowMeans(all_preds)
  
  # Prepare a data frame to plot response curves
  plot_response <- data.frame(environmental_values = dummy_data[,m], predicted_values = ensemble_preds,
                              sd = ensemble_sd)
  
  # Add bounds to the data frame
  plot_response$upper_bound <- plot_response$predicted_values + plot_response$sd
  plot_response$lower_bound <- plot_response$predicted_values - plot_response$sd
  
  # Plot the response curve
  ggplot(plot_response, aes(x = environmental_values, y = predicted_values)) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", size = 0.65, se = FALSE) +
    # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.65, se = FALSE) +
    geom_ribbon(aes(ymin = predicted_values - sd, ymax = predicted_values + sd), 
                fill = "grey24", alpha = 0.2) +
    ylim(pmin(0, min(plot_response$lower_bound, na.rm = TRUE)), pmax(1, max(plot_response$upper_bound, na.rm = TRUE))) +
    xlim(min(plot_response$environmental_values), max(plot_response$environmental_values)) +
    labs(title = paste0(m), x = "Environmental values", y = "Predicted values") +
    theme_minimal() +
    theme(axis.title = element_text(size = 8), axis.text = element_text(size = 5))
  
  ggsave(paste0("output_data/plots/response_curves/C_pipiens_response_curve_smoothed_",m,".png"), width = 5, height = 5, units = "cm")
  
} # Close the loop over the predictors




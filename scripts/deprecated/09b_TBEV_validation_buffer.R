# ZOE project 
# Disease phenology analysis of TBEV in Europe 

# ---------------------------------------------------------------------- #
#                          09a. Model validation                         #
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
load("output_data/models/TBEV_SDMs_buffer.RData")

load("output_data/models/TBEV_SDMs_buffer_presenceabsence.RData")



#-------------------------------------------------------------------------------

# 1. 5-fold cross-validation ---------------------------------------------------

# Part the presence-background dataset into 5 folds
kfolds <- 5
ks <- dismo::kfold(seq_len(nrow(TBEV_occ_env)), k = kfolds)


# GLM
print("GLM")

# Initialise a list to store the performance measures for each model
glm_performances <- list()

# Create a matrix to store the cross-validated predictions of the three models
m_glm_preds_cv_all <- matrix(nrow=nrow(TBEV_occ_env), ncol=length(models_glm))

for (glm_index in 1:length(models_glm)) { # Start of the loop over all three models
  
  print(paste("Evaluating GLM", glm_index))
  m_glm_preds_cv <- rep(NA, nrow(TBEV_occ_env))
  
  # Extract predictors from model
  model_name <- names(models_glm)[glm_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)) { # Start of the loop over k folds
    cv_train <- TBEV_occ_env[ks != i,]
    cv_test <- TBEV_occ_env[ks == i,]
    cv_weights <- weights[ks != i]
    
    cv_glm <- update(models_glm[[glm_index]], data = cv_train, weights = cv_weights)
    m_glm_preds_cv[ks == i] <- predict(cv_glm, cv_test, type = "response")
  }
  
  
  # Calculate performance measures
  m_glm_perf_cv <- evalSDM(TBEV_occ_env$occ, m_glm_preds_cv, weights = weights)
  
  # Calculate Boyce index with smoothing methods
  m_glm_preds_cv_presences <- m_glm_preds_cv[presences] # Just retain the predictions of the presences based on indices
  absences <- which(TBEV_occ_env$occ == 0) # Extract the position index of background data
  m_glm_preds_cv_absences <- m_glm_preds_cv[absences] # Just retain the predictions of the absences based on indices
  m_glm_boyce_cv <- sfbi(m_glm_preds_cv_presences, m_glm_preds_cv_absences, ktry = 10) # Apply sfbi function
  
  boyce_index_m_glm <- m_glm_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers
  
  # Add the Boyce index to the performance metrics data frame
  m_glm_perf_cv$Boyce <- boyce_index_m_glm
  
  # Store the performance measures for the respective models
  glm_performances[[glm_index]] <- m_glm_perf_cv
  
  # Store the cross-validated predictions of the three models into the prepared data frame
  m_glm_preds_cv_all[, glm_index] <- m_glm_preds_cv
  
} # Close the loop over all three models

# Extract average performance measures
avg_glm_performances <- sapply(c("AUC", "TSS", "Kappa", "Sens", "Spec", "PCC", "D2", "thresh", "Boyce"), function(col) mean(sapply(glm_performances, function(df) df[[col]])))
avg_glm_performances <- as.data.frame(t(avg_glm_performances))




# GAM
print("GAM")

# Initialise a list to store the performance measures for each model
gam_performances <- list()

# Create a matrix to store the cross-validated predictions of the three models
m_gam_preds_cv_all <- matrix(nrow=nrow(TBEV_occ_env), ncol=length(models_gam))

for (gam_index in 1:length(models_gam)) { # Start of the loop over all three models
  
  print(paste("Evaluating GAM", gam_index))
  m_gam_preds_cv <- rep(NA, nrow(TBEV_occ_env))
  
  # Extract predictors from model
  model_name <- names(models_gam)[gam_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)) { # Start of the loop over k folds
    cv_train <- TBEV_occ_env[ks != i,]
    cv_test <- TBEV_occ_env[ks == i,]
    cv_weights <- weights[ks != i]
    
    cv_gam <- update(models_gam[[gam_index]], data = cv_train, weights = cv_weights)
    m_gam_preds_cv[ks == i] <- predict(cv_gam, cv_test, type = "response")
  }
  
  
  # Calculate performance measures
  m_gam_perf_cv <- evalSDM(TBEV_occ_env$occ, m_gam_preds_cv, weights = weights)
  
  # Calculate Boyce index with smoothing methods
  m_gam_preds_cv_presences <- m_gam_preds_cv[presences] # Just retain the predictions of the presences based on indices
  absences <- which(TBEV_occ_env$occ == 0) # Extract the position index of background data
  m_gam_preds_cv_absences <- m_gam_preds_cv[absences] # Just retain the predictions of the absences based on indices
  m_gam_boyce_cv <- sfbi(m_gam_preds_cv_presences, m_gam_preds_cv_absences, ktry = 10) # Apply sfbi function
  
  boyce_index_m_gam <- m_gam_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers
  
  # Add the Boyce index to the performance metrics data frame
  m_gam_perf_cv$Boyce <- boyce_index_m_gam
  
  # Store the performance measures for the respective models
  gam_performances[[gam_index]] <- m_gam_perf_cv
  
  # Store the cross-validated predictions of the three models into the prepared data frame
  m_gam_preds_cv_all[, gam_index] <- m_gam_preds_cv
  
} # Close the loop over all three models

# Extract average performance measures
avg_gam_performances <- sapply(c("AUC", "TSS", "Kappa", "Sens", "Spec", "PCC", "D2", "thresh", "Boyce"), function(col) mean(sapply(gam_performances, function(df) df[[col]])))
avg_gam_performances <- as.data.frame(t(avg_gam_performances))





# RF
print("RF")

# Initialise a list to store the performance measures for each model
rf_performances <- list()

# Create a matrix to store the cross-validated predictions of the three models
m_rf_preds_cv_all <- matrix(nrow=nrow(TBEV_occ_env), ncol=length(models_rf))


for (rf_index in 1:length(models_rf)) { # Start of the loop over all models
  
  print(paste("Evaluating RF", rf_index))
  
  m_rf_preds_matrix_cv <- matrix(nrow=nrow(TBEV_occ_env), ncol=background_presence_ratio)
  
  # Extract predictors from model
  model_name <- names(models_rf)[rf_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)){ # Start of the loop over k folds
    print(i)
    hold_in <- which(ks != i)
    hold_out <- which(ks == i)
    
    for (y in 1:background_presence_ratio) { # Start the loop over 10 models
      print(y)
      m_rf <- models_rf[[rf_index]]
      
      data_index <- c(presences, which(TBEV_occ_env$abs_index == y))
      
      cv_train <- TBEV_occ_env[hold_in[hold_in %in% data_index],]
      cv_test <- TBEV_occ_env[hold_out,]
      
      cv_rf <- update(m_rf[[y]], data = cv_train)
      m_rf_preds_matrix_cv[hold_out, y] <- predict(cv_rf, cv_test, type = "response")
    }
  }
  
  m_rf_preds_cv <- rowMeans(m_rf_preds_matrix_cv)
  
  # Calculate performance measures
  m_rf_perf_cv <- evalSDM(TBEV_occ_env$occ, m_rf_preds_cv, weights = weights)
  
  # Calculate Boyce index with smoothing methods
  m_rf_preds_cv_presences <- m_rf_preds_cv[presences] # Just retain the predictions of the presences based on indices
  absences <- which(TBEV_occ_env$occ == 0) # Extract the position index of background data
  m_rf_preds_cv_absences <- m_rf_preds_cv[absences] # Just retain the predictions of the absences based on indices
  m_rf_boyce_cv <- sfbi(m_rf_preds_cv_presences, m_rf_preds_cv_absences, ktry = 10) # Apply sfbi function
  
  boyce_index_m_rf <- m_rf_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers
  
  # Add the Boyce index to the performance metrics data frame
  m_rf_perf_cv$Boyce <- boyce_index_m_rf
  
  # Store the performance measures for the respective models
  rf_performances[[rf_index]] <- m_rf_perf_cv
  
  # Store the cross-validated predictions of the three models into the prepared data frame
  m_rf_preds_cv_all[, rf_index] <- m_rf_preds_cv
  
} # Close the loop over all three models

# Extract average performance measures
avg_rf_performances <- sapply(c("AUC", "TSS", "Kappa", "Sens", "Spec", "PCC", "D2", "thresh", "Boyce"), function(col) mean(sapply(rf_performances, function(df) df[[col]])))
avg_rf_performances <- as.data.frame(t(avg_rf_performances))





# BRT
print("BRT")

# Initialise a list to store the performance measures for each model
brt_performances <- list()

# Create a matrix to store the cross-validated predictions of the three models
m_brt_preds_cv_all <- matrix(nrow=nrow(TBEV_occ_env), ncol=length(models_brt))

for (brt_index in 1:length(models_brt)) { # Start of the loop over all models
  
  print(paste("Evaluating BRT", brt_index))
  
  m_brt_preds_matrix_cv <- matrix(nrow=nrow(TBEV_occ_env), ncol=background_presence_ratio)
  
  # Extract predictors from model
  model_name <- names(models_brt)[brt_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)){ # Start of the loop over k folds
    print(i)
    hold_in <- which(ks != i)
    hold_out <- which(ks == i)
    
    for (y in 1:background_presence_ratio) { # Start the loop over 10 models
      print(y)
      m_brt <- models_brt[[brt_index]]
      
      data_index <- c(presences, which(TBEV_occ_env$abs_index == y))
      
      cv_train <- TBEV_occ_env[hold_in[hold_in %in% data_index],]
      names(cv_train)[names(cv_train)=='occ'] <- m_brt[[y]]$response.name
      cv_test <- TBEV_occ_env[hold_out,]
      
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
  m_brt_perf_cv <- evalSDM(TBEV_occ_env$occ, m_brt_preds_cv, weights = weights)
  
  # Calculate Boyce index with smoothing methods
  m_brt_preds_cv_presences <- m_brt_preds_cv[presences] # Just retain the predictions of the presences based on indices
  absences <- which(TBEV_occ_env$occ == 0) # Extract the position index of background data
  m_brt_preds_cv_absences <- m_brt_preds_cv[absences] # Just retain the predictions of the absences based on indices
  m_brt_boyce_cv <- sfbi(m_brt_preds_cv_presences, m_brt_preds_cv_absences, ktry = 10) # Apply sfbi function
  
  boyce_index_m_brt <- m_brt_boyce_cv[6] # Extract values of smoothed Boyce index using the mean of all smoothers
  
  # Add the Boyce index to the performance metrics data frame
  m_brt_perf_cv$Boyce <- boyce_index_m_brt
  
  # Store the performance measures for the respective models
  brt_performances[[brt_index]] <- m_brt_perf_cv
  
  # Store the cross-validated predictions of the three models into the prepared data frame
  m_brt_preds_cv_all[, brt_index] <- m_brt_preds_cv
  
} # Close the loop over all three models

# Extract average performance measures
avg_brt_performances <- sapply(c("AUC", "TSS", "Kappa", "Sens", "Spec", "PCC", "D2", "thresh", "Boyce"), function(col) mean(sapply(brt_performances, function(df) df[[col]])))
avg_brt_performances <- as.data.frame(t(avg_brt_performances))




#-------------------------------------------------------------------------------

# 2. Make ensemble prediction --------------------------------------------------

# Make ensemble
print("Ensemble")

# Calculate the average cross-validated predictions of the three models per algorithm
m_glm_preds_cv_all_avg <- rowMeans(m_glm_preds_cv_all)
m_gam_preds_cv_all_avg <- rowMeans(m_gam_preds_cv_all)
m_rf_preds_cv_all_avg <- rowMeans(m_rf_preds_cv_all)
m_brt_preds_cv_all_avg <- rowMeans(m_brt_preds_cv_all)

# Calculate mean cross-validated prediction over all algorithms (create ensemble)
# and assess performance measures
m_ens_preds_cv <- rowMeans(data.frame(m_glm_preds_cv_all_avg, m_gam_preds_cv_all_avg, m_brt_preds_cv_all_avg, m_rf_preds_cv_all_avg))
m_ens_perf_cv <- evalSDM(TBEV_occ_env$occ, m_ens_preds_cv, weights = weights)
m_ens_preds_cv_presences <- m_ens_preds_cv[presences]
absences <- which(TBEV_occ_env$occ == 0)
m_ens_preds_cv_absences <- m_ens_preds_cv[absences]
m_ens_boyce_cv <- sfbi(m_ens_preds_cv_presences, m_ens_preds_cv_absences, ktry = 10)
boyce_index_m_ens <- m_ens_boyce_cv[6]
m_ens_perf_cv$Boyce <- boyce_index_m_ens

# Compare cross-validated model performance across algorithms
comp_perf <- rbind(glm = avg_glm_performances, gam = avg_gam_performances, rf = avg_rf_performances, brt = avg_brt_performances, mean_prob = m_ens_perf_cv)

# Add a column containing the names of the algorithm
comp_perf <- data.frame(alg=rownames(comp_perf),comp_perf)


# Save cross-validation outcomes
save(m_glm_preds_cv_all, glm_performances, m_gam_preds_cv_all, gam_performances,
     m_rf_preds_cv_all, rf_performances, m_brt_preds_cv_all, brt_performances, 
     m_ens_preds_cv, m_ens_perf_cv, ks, comp_perf,
     file = "output_data/validation/TBEV_validation_buffer_presenceabsence.RData")



#-------------------------------------------------------------------------------

# 3. Monthly model performances ------------------------------------------------

# Read in the needed data
load("output_data/models/TBEV_SDMs_buffer.RData") # The occurrence data frame
load("output_data/validation/TBEV_validation_buffer_thresh.RData") # The cross-validated ensemble predictions

# Identify the earliest and latest month of observation
start_month <- min(TBEV_occ_env$month)
end_month <- max(TBEV_occ_env$month)


# Create a vector containing the months of a year
month <- str_pad(paste0(start_month:end_month), width = 2, pad = "0")


for (m in month) { # Start of the loop over all months
  
  print(m)
  
  # Subset the occurrence data frame to only contain data of the respective month
  TBEV_occ_env_month <- subset(TBEV_occ_env, TBEV_occ_env$month == m)
  
  # Get the indices of the presences and absences of the respective month
  monthly_data_indices <- which(TBEV_occ_env$month == m)
  
  # Subset the cross-validated ensemble predictions by the extracted indices for that month
  monthly_data <- m_ens_preds_cv[monthly_data_indices]
  
  # Extract the weights of the respective monthly cross-validations
  monthly_weights <- weights[monthly_data_indices]
  
  # Calculate the performance measures (AUC, TSS, Sensitivity, Specificity) for the respective month
  m_ens_perf_cv_month <- evalSDM(TBEV_occ_env_month$occ, monthly_data, weights = monthly_weights)
  
  # Calculate the Boyce index with smoothing methods for the respective month
  presences_month <- which(TBEV_occ_env_month$occ == 1)
  m_ens_preds_cv_presences_month <- monthly_data[presences_month]
  absences_month <- which(TBEV_occ_env_month$occ == 0)
  m_ens_preds_cv_absences_month <- monthly_data[absences_month] 
  m_ens_boyce_cv_month  <- sfbi(m_ens_preds_cv_presences_month, m_ens_preds_cv_absences_month, ktry = 10) 
  boyce_index_m_ens_month <- m_ens_boyce_cv_month[6]
  
  # Add the Boyce index to the performance metrics data frame
  m_ens_perf_cv_month$Boyce <- boyce_index_m_ens_month
  
  
  # Save the data frame with performance measures
  save(m_ens_perf_cv_month,  file = paste0("output_data/validation/TBEV_monthly_validation_",m,"_buffer.RData"))
  
  
} # Close the loop over all months




#-------------------------------------------------------------------------------

# 4. Create response curves of ensemble predictions ----------------------------

# Read in the needed data
load("output_data/models/TBEV_SDMs_buffer.RData") # The different models and selected predictor variables

# Create a vector with all predictors within the models, soley the temperature
# variables, and the remaining variables that are used in all models
my_preds_all <- c("tas", "tasmin", "tasmax", "primary_openland", "pasture", "pr", "hurs", "rangeland", "primary_forest", "urban", "secondary_forest", "cropland", "I_ricinus")
temp_var <- c("tas", "tasmin", "tasmax")
remain_var <- c("primary_openland", "pasture", "pr", "hurs", "rangeland", "primary_forest", "urban", "secondary_forest", "cropland", "I_ricinus")


for (p in my_preds_all) { # Loop through all predictor variables
  
  # Name the other predictors in a vector
  my_preds_minus_p <- setdiff(my_preds_all, p)
  
  # Create an environmental dummy dataset (keeping the other predictors at their mean)
  dummy_data <- data.frame(seq(min(TBEV_occ_env[,p], na.rm = TRUE), max(TBEV_occ_env[,p], na.rm = TRUE), length = 100))
  for (pred in my_preds_minus_p) {
    dummy_data[[pred]] <- mean(TBEV_occ_env[[pred]], na.rm = TRUE)
  }
  
  names(dummy_data)[1] <- p
  
  # Create a matrix to store the predictions for the response plots per variable,
  # for the temperature variables solely one column is created as they are just 
  # used in one model each, for the other variables we additionally create a matrix
  # to store the ensemble sd values resulting from the application of the different algorithms
  if (p %in% temp_var) { print(p)
  } else if (p %in% remain_var) { ens_preds_p <- matrix(nrow=nrow(dummy_data), ncol=length(models_glm))
  ens_sd_p <- matrix(nrow=nrow(dummy_data), ncol=length(models_glm))
  print(p)
  }
  
  # Use an if condition to assign the different temperature variables to the
  # respective model numbers
  if (p == "tas") { s <- 1
  print(s)
  } else if (p == "tasmin") { s <- 2
  print(s)
  } else if (p == "tasmax") { s <- 3
  print(s)
  }
  
  # Start the predictions to dummy data based on one models for each temperature variable
  if (p %in% temp_var) { 
    
    # Extract the predictors within that model (use GLM as example model)
    model_name <- names(models_glm)[s]
    my_preds <- unlist(strsplit(model_name, "\\+"))
    
    print("GLM")
    response_preds_glm <- predictSDM(models_glm[[s]], dummy_data[, my_preds])
    
    print("GAM") 
    response_preds_gam <- predictSDM(models_gam[[s]], dummy_data[, my_preds])
    
    print("RF")
    response_preds_rf_all <- matrix(nrow=nrow(dummy_data[, my_preds]), ncol=background_presence_ratio)
    
    for (y in 1:background_presence_ratio) {
      print(y)
      response_preds_rf_all[,y] <- predictSDM(models_rf[[s]][[y]], dummy_data[, my_preds])
    }
    
    response_preds_rf <- rowMeans(response_preds_rf_all)
    
    print("BRT")
    response_preds_brt_all <- matrix(nrow=nrow(dummy_data), ncol=background_presence_ratio)
    
    for (y in 1:background_presence_ratio) { 
      print(y)
      response_preds_brt_all[,y] <- predictSDM(models_brt[[s]][[y]], dummy_data[, my_preds])
    }
    response_preds_brt <- rowMeans(response_preds_brt_all)
    
    # Combine the predictions into an ensemble and calculate the sd
    all_preds <- cbind(response_preds_glm, response_preds_gam, response_preds_rf, response_preds_brt)
    ens_sd_p <- apply(all_preds, 1, sd)
    ens_preds_p <- rowMeans(all_preds)
    
    
  } else if (p %in% remain_var) { # If the predictor variable is used in all three models, we start looping
    # over the number of different models
    
    for (s in 1:length(models_glm)) { # Start to loop through the number of models with different sets of variables
      
      print(s)
      
      # Extract the predictors within that model (use GLM as example model)
      model_name <- names(models_glm)[s]
      my_preds <- unlist(strsplit(model_name, "\\+"))
      
      # Generate predictions to dummy data for each model
      print("GLM")
      response_preds_glm <- predictSDM(models_glm[[s]], dummy_data[, my_preds])
      
      print("GAM") 
      response_preds_gam <- predictSDM(models_gam[[s]], dummy_data[, my_preds])
      
      print("RF")
      response_preds_rf_all <- matrix(nrow=nrow(dummy_data), ncol=background_presence_ratio)
      
      for (y in 1:background_presence_ratio) {
        print(y)
        response_preds_rf_all[,y] <- predictSDM(models_rf[[s]][[y]], dummy_data[, my_preds])
      }
      
      response_preds_rf <- rowMeans(response_preds_rf_all)
      
      print("BRT")
      response_preds_brt_all <- matrix(nrow=nrow(dummy_data), ncol=background_presence_ratio)
      
      for (y in 1:background_presence_ratio) { 
        print(y)
        response_preds_brt_all[,y] <- predictSDM(models_brt[[s]][[y]], dummy_data[, my_preds])
      }
      response_preds_brt <- rowMeans(response_preds_brt_all)
      
      # Combine the predictions into an ensemble and calculate the sd
      all_preds <- cbind(response_preds_glm, response_preds_gam, response_preds_rf, response_preds_brt)
      ensemble_sd <- apply(all_preds, 1, sd)
      ensemble_preds <- rowMeans(all_preds)
      
      # Insert the predictions into the created results matrix
      ens_preds_p[,s] <- ensemble_preds
      ens_sd_p[,s] <- ensemble_sd
      
      
    } # Close the loop over the three different models
    
    # Calculate the mean predictions over the three different models
    ens_preds_p <- rowMeans(ens_preds_p)
    ens_sd_p <- rowMeans(ens_sd_p)
    
  } # Close if-condition
  
  # Prepare a data frame to plot response curves
  plot_response <- data.frame(environmental_values = dummy_data[,p], predicted_values = ens_preds_p,
                              sd = ens_sd_p)
  
  # Add bounds to the data frame
  plot_response$upper_bound <- plot_response$predicted_values + plot_response$sd
  plot_response$lower_bound <- plot_response$predicted_values - plot_response$sd
  
  
  # Plot the response curve
  ggplot(plot_response, aes(x = environmental_values, y = predicted_values)) +
    geom_smooth(method = "lm", formula = y ~ x + I(x^2), color = "black", linewidth = 0.65, se = FALSE) +
    # geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black", size = 0.65, se = FALSE) +
    geom_ribbon(aes(ymin = predicted_values - sd, ymax = predicted_values + sd), 
                fill = "grey24", alpha = 0.2) +
    ylim(pmin(0, min(plot_response$lower_bound, na.rm = TRUE)), pmax(1, max(plot_response$upper_bound, na.rm = TRUE))) +
    xlim(min(plot_response$environmental_values), max(plot_response$environmental_values)) +
    labs(title = paste0(p), x = "Environmental values", y = "Predicted values") +
    theme_minimal() +
    theme(axis.title = element_text(size = 8), axis.text = element_text(size = 5))
  
  ggsave(paste0("output_data/plots/response_curves/TBEV_response_curve_smoothed_",p,"_buffer_presenceabsence.png"), width = 5, height = 5, units = "cm")
  
  
} # Close the loop over the different predictor variables




#-------------------------------------------------------------------------------

# 5. Create scatterplot of cross-validated ensemble predictions (I_ricinus) ----


my_preds_all <- c("tas", "tasmin", "tasmax", "primary_openland", "pasture", "pr", "hurs", "rangeland", "primary_forest", "urban", "secondary_forest", "cropland", "I_ricinus")

TBEV_occ_env_scat <- TBEV_occ_env

# Add the cross-validated ensemble predictions as column to the data frame containing 
# the predictor values per presence/absence location
TBEV_occ_env_scat$ens_preds_cv <- m_ens_preds_cv


for (p in my_preds_all) { # Loop through all predictor variables
  
  print(p)

  # Plot as scatterplot
  ggplot(TBEV_occ_env_scat, aes(x = TBEV_occ_env_scat[,p], y = ens_preds_cv, color = factor(occ))) +
    geom_point(alpha = 0.5, size = 0.005) + 
    geom_smooth(method = "loess", color = "black", se = FALSE, linewidth = 0.65) + 
    scale_color_manual(values = c("firebrick4", "darkblue"), labels = c("Pseudo-Absence", "Presence")) +
    labs(title = paste0(p), x = "Environmental values", y = "CV Ensemble Predictions",
         color = "Occurrence") +
    theme_minimal()  +
    theme(axis.title = element_text(size = 8), axis.text = element_text(size = 5))
  
  ggsave(paste0("output_data/plots/response_curves/TBEV_response_scatterplot_",p,"_buffer_presenceabsence.png"), width = 10, height = 5, units = "cm")


} # Close the loop over predictors

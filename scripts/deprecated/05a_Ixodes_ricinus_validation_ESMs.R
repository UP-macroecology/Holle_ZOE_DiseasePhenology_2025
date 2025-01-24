# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                         05a. Model validation                          #
# ---------------------------------------------------------------------- #

# Load needed packages
library(mgcv)
library(maxnet)
library(randomForest)
library(gbm)
library(dismo)
library(ggplot2)
library(PresenceAbsence)

# Load needed objects
source("scripts/00_functions.R") # Get the function for SDM evaluation

# Read in SDM models
load("output_data/models/I_ricinus_SDMs.RData")




#-------------------------------------------------------------------------------

# 1. 5-fold cross-validation ---------------------------------------------------

# Part the presence-background dataset into 5 folds
kfolds <- 5
ks <- dismo::kfold(seq_len(nrow(I_ricinus_occ_env)), k = kfolds)





# GLM
print("GLM")

# Initialise a list to store the performance measures for each model
glm_performances <- list()

for (glm_index in 1:length(models_glm)) { # Start of the loop over all models
  
  print(paste("Evaluating GLM", glm_index))
  m_glm_preds_cv <- rep(NA, nrow(I_ricinus_occ_env))
  
  # Extract predictors from model
  model_name <- names(models_glm)[glm_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)) { # Start of the loop over k folds
    cv_train <- I_ricinus_occ_env[ks != i,]
    cv_test <- I_ricinus_occ_env[ks == i,]
    cv_weights <- weights[ks != i]
    
    cv_glm <- update(models_glm[[glm_index]], data = cv_train, weights = cv_weights)
    m_glm_preds_cv[ks == i] <- predict(cv_glm, cv_test, type = "response")
  }
  
  # Calculate performance measures
  m_glm_perf_cv <- evalSDM(I_ricinus_occ_env$occ, m_glm_preds_cv)
  
  # Store the performance measures for the respective models
  glm_performances[[glm_index]] <- m_glm_perf_cv
  
}

# Extract the AUC performance values of all models
AUC_all_glm <- sapply(glm_performances, function(perf) perf$AUC)

# Calculate the mean of the performance measure
AUC_mean_performance_glm <- mean(AUC_all_glm, na.rm = TRUE)




# GAM
print("GAM")

# Initialise a list to store the performance measures for each model
gam_performances <- list()

for (gam_index in 1:length(models_gam)) { # Start of the loop over all models
  
  print(paste("Evaluating GAM", gam_index))
  m_gam_preds_cv <- rep(NA, nrow(I_ricinus_occ_env))
  
  # Extract predictors from model
  model_name <- names(models_gam)[gam_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)) { # Start of the loop over k folds
    cv_train <- I_ricinus_occ_env[ks != i,]
    cv_test <- I_ricinus_occ_env[ks == i,]
    cv_weights <- weights[ks != i]
    
    cv_gam <- update(models_gam[[gam_index]], data = cv_train, weights = cv_weights)
    m_gam_preds_cv[ks == i] <- predict(cv_gam, cv_test, type = "response")
  }
  
  # Calculate performance measures
  m_gam_perf_cv <- evalSDM(I_ricinus_occ_env$occ, m_gam_preds_cv)
  
  # Store the performance measures for the respective models
  gam_performances[[gam_index]] <- m_gam_perf_cv
  
}

# Extract the AUC performance values of all models
AUC_all_gam <- sapply(gam_performances, function(perf) perf$AUC)

# Calculate the mean of the performance measure
AUC_mean_performance_gam <- mean(AUC_all_gam, na.rm = TRUE)




# RF
print("RF")

# Initialise a list to store the performance measures for each model
rf_performances <- list()

for (rf_index in 1:length(models_rf)) { # Start of the loop over all models
  
  print(paste("Evaluating RF", rf_index))
  
  m_rf_preds_matrix_cv <- matrix(nrow=nrow(I_ricinus_occ_env), ncol=4)
  
  # Extract predictors from model
  model_name <- names(models_rf)[rf_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  for(i in seq_len(kfolds)){ # Start of the loop over k folds
    print(i)
    hold_in <- which(ks != i)
    hold_out <- which(ks == i)
    
    for (y in 1:4) { # Start the loop over 10 models
      print(y)
      m_rf <- models_rf[[rf_index]]
      
      data_index <- c(presences, which(I_ricinus_occ_env$abs_index == y))
      
      cv_train <- I_ricinus_occ_env[hold_in[hold_in %in% data_index],]
      cv_test <- I_ricinus_occ_env[hold_out,]
      
      cv_rf <- update(m_rf[[y]], data = cv_train)
      m_rf_preds_matrix_cv[hold_out, y] <- predict(cv_rf, cv_test, type = "response")
    }
  }
  
  m_rf_preds_cv <- rowMeans(m_rf_preds_matrix_cv)
  
  # Calculate performance measures
  m_rf_perf_cv <- evalSDM(I_ricinus_occ_env$occ, m_rf_preds_cv)
  
  # Store the performance measures for the respective models
  rf_performances[[rf_index]] <- m_rf_perf_cv
  
}

# Extract the AUC performance values of all models
AUC_all_rf <- sapply(rf_performances, function(perf) perf$AUC)

# Calculate the mean of the performance measure
AUC_mean_performance_rf <- mean(AUC_all_rf, na.rm = TRUE)




# BRT
print("BRT")

# Initialise a list to store the performance measures for each model
brt_performances <- list()

for (brt_index in 1:length(models_brt)) { # Start of the loop over all models
  
  print(paste("Evaluating BRT", brt_index))
  
  m_brt_preds_matrix_cv <- matrix(nrow=nrow(I_ricinus_occ_env), ncol=4)
  
  # Extract predictors from model
  model_name <- names(models_brt)[brt_index]
  my_preds <- unlist(strsplit(model_name, "\\+"))
  print(my_preds)
  
  
  for(i in seq_len(kfolds)){ # Start of the loop over k folds
    print(i)
    hold_in <- which(ks != i)
    hold_out <- which(ks == i)
    
    for (y in 1:4) { # Start the loop over 10 models
      print(y)
      m_brt <- models_brt[[brt_index]]
      
      data_index <- c(presences, which(I_ricinus_occ_env$abs_index == y))
      
      cv_train <- I_ricinus_occ_env[hold_in[hold_in %in% data_index],]
      names(cv_train)[names(cv_train)=='occ'] <- m_brt[[y]]$response.name
      cv_test <- I_ricinus_occ_env[hold_out,]
      
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
  m_brt_perf_cv <- evalSDM(I_ricinus_occ_env$occ, m_brt_preds_cv)
  
  # Store the performance measures for the respective models
  brt_performances[[brt_index]] <- m_brt_perf_cv
  
}

# Extract the AUC performance values of all models
AUC_all_brt <- sapply(brt_performances, function(perf) perf$AUC)

# Calculate the mean of the performance measure
AUC_mean_performance_brt <- mean(AUC_all_brt, na.rm = TRUE)



# Save performances
save(glm_performances, gam_performances, rf_performances, brt_performances, ks, 
     file = "output_data/validation/I_ricinus_validation_ESMs.RData")

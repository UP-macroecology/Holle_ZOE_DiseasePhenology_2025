# ZOE project 
# Disease phenology analysis of West Nile Fever in Europe

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                          08b. Model fitting                            #
# ---------------------------------------------------------------------- #


# Load needed packages
library(mgcv)
library(maxnet)
library(randomForest)
library(gbm)
library(dismo)
library(tidyverse)
library(corrplot)

# Load needed objects
source("scripts/00_functions.R") # Get the select07_cv function (explained deviance function)

# Read in presence and background data
load("output_data/data/WNF_occ_env_vector.RData")



#-------------------------------------------------------------------------------

# 1. Variable selection --------------------------------------------------------

# Retrieve predictors (excluding predictors related to temperature, as we will include
# all of these in different models, because tmin, tmean, and tmax might be more relevant
# for different months of a year. We will include them in different models because they are
# highly correlated.)
predictors <- names(WNF_occ_env[, c(7, 11:20)])

# Check for collinearity in correlation matrix
cor_mat <- cor(WNF_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Generate weights
wgt <- rep(1, times = nrow(WNF_occ_env))

# Run select07_cv function
var_sel <- select07_cv(X = WNF_occ_env[,predictors], 
                       y = WNF_occ_env$occ, 
                       threshold = 0.7,
                       weights = wgt)

# Extract most important and weakly correlated predictors
my_preds <- var_sel$pred_sel

# Create a list with three different vectors combining the three different 
# temperature variables (tas, tasmin, tasmax) separately with the selected 
# variables
my_preds_list <- list(tas_mypreds = c(my_preds, "tas"), tasmin_mypreds = c(my_preds, "tasmin"),
                      tasmax_mypreds = c(my_preds, "tasmax"))






#-------------------------------------------------------------------------------

# 2. Model fitting -------------------------------------------------------------

# Fit GLM (including linear and quadratic terms, AIC-based stepwise variable selection, equal weights)
print("GLM")

# Create a list to store the models
models_glm <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_glm <- step(glm(as.formula(paste('occ~',paste(c(my_preds, paste0('I(',my_preds,'^2)')), collapse ='+'))),
                    family='binomial', data = WNF_occ_env, weights = wgt))
  
  models_glm[[m]] <- m_glm
  
} # Close the loop over the three different predictor combinations

names(models_glm) <- sapply(my_preds_list, paste, collapse = "+")



# Fit GAM (cubic smoothing splines, equal weights)
print("GAM")

# Create a list to store the models
models_gam <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                     family='binomial', data = WNF_occ_env, weights = wgt)
  
  models_gam[[m]] <- m_gam
  
} # Close the loop over the three different predictor combinations

names(models_gam) <- sapply(my_preds_list, paste, collapse = "+")



# Fit RF (same number of presences and background data, ten models in total)
print("RF")

models_rf <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_rf <- randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                       data = WNF_occ_env, ntree = 1000, nodesize = 10, importance = T)
  
  
  models_rf[[m]] <- m_rf
  
} # Close the loop over the three different predictor combinations

names(models_rf) <- sapply(my_preds_list, paste, collapse = "+")



# Fit BRT (same number of presences and background data, ten models in total, adaptable learning rate to fit model with 1000 and 10000 trees)
print("BRT")

models_brt <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  opt.LR <- TRUE
  LR = 0.01
  
  while(opt.LR){
    m_brt <- try(gbm.step(data = WNF_occ_env, gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
    if (class(m_brt) == "try-error" | class(m_brt) == "NULL"){
      LR <- LR/2
    } else
      if(m_brt$gbm.call$best.trees<1000){
        LR <- LR/2
      } else 
        if(m_brt$gbm.call$best.trees>10000){
          LR <- LR*2
        } else { 
          opt.LR <- FALSE}}
  
  models_brt[[m]] <- m_brt
  
  
} # Close the loop over the three different predictor combinations

names(models_brt) <- sapply(my_preds_list, paste, collapse = "+")



# Save the models
save(models_glm, models_gam, models_rf, models_brt, predictors, my_preds_list, WNF_occ_env, 
     file = "output_data/models/WNF_SDMs_vector.RData")

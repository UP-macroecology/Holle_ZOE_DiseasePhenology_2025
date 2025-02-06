# ZOE project 
# Disease phenology analysis of Culex pipiens in Europe (primary transmitter of WNV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                          04b. Model fitting                            #
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
load("output_data/data/C_pipiens_occ_env.RData")



#-------------------------------------------------------------------------------

# 1. Variable selection --------------------------------------------------------

# Retrieve predictors (excluding predictors related to temperature, as we will include
# all of these in different models, because tmin, tmean, and tmax might be more relevant
# for different months of a year. We will include them in different models because they are
# highly correlated.)
predictors <- names(C_pipiens_occ_env[, c(7, 11:19)])

# Check for collinearity in correlation matrix
cor_mat <- cor(C_pipiens_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Define equal weights for presences and background data
weights <- ifelse(C_pipiens_occ_env$occ == 1, 1, sum(C_pipiens_occ_env$occ == 1) / sum(C_pipiens_occ_env$occ == 0))

# Run select07_cv function
var_sel <- select07_cv(X = C_pipiens_occ_env[,predictors], 
                       y = C_pipiens_occ_env$occ, 
                       threshold = 0.7,
                       weights = weights)

# Extract most important and weakly correlated predictors
my_preds <- var_sel$pred_sel

# Create a list with three different vectors combining the three different 
# temperature variables (tas, tasmin, tasmax) separately with the selected 
# variables
my_preds_list <- list(tas_mypreds = c(my_preds, "tas"), tasmin_mypreds = c(my_preds, "tasmin"),
                      tasmax_mypreds = c(my_preds, "tasmax"))

# Extract the ratio of presence to background data to know how many machine learning models we can build
background_presence_ratio <- round(sum(C_pipiens_occ_env$occ == 0) / sum(C_pipiens_occ_env$occ == 1), 0)


# Machine-learning methods should be fitted with equal number of presences and background and repeated 10 times (Barbet-Massin et al. (2012))
# As we don't have enough background points, we fit the amount of machine-learning models that our data allows for
presences <- which(C_pipiens_occ_env$occ == 1) # Extract the position index of presences
C_pipiens_occ_env$abs_index <- NA # Create a new column to insert the absence index
C_pipiens_occ_env$abs_index[C_pipiens_occ_env$occ!=1] <- sample(1:background_presence_ratio, sum(C_pipiens_occ_env$occ!=1), replace = TRUE) # Insert sampled numbers of background-presence-ratio




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
                    family='binomial', data = C_pipiens_occ_env, weights = weights))
  
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
                     family='binomial', data = C_pipiens_occ_env, weights = weights)
  
  models_gam[[m]] <- m_gam
  
} # Close the loop over the three different predictor combinations

names(models_gam) <- sapply(my_preds_list, paste, collapse = "+")




# Fit RF (same number of presences and background data, ten models in total)
print("RF")

models_rf <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_rf <- lapply(1:background_presence_ratio, FUN=function(i){sp_train <- C_pipiens_occ_env[c(presences, which(C_pipiens_occ_env$abs_index == i)),]; print(i); randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                                                                                                                                                                            data = sp_train, ntree = 1000, nodesize = 10, importance = T)})
  
  models_rf[[paste(my_preds, collapse = "+")]] <- m_rf
  
} # Close the loop over the three different predictor combinations




# Fit BRT (same number of presences and background data, ten models in total, adaptable learning rate to fit model with 1000 and 10000 trees)
print("BRT")

models_brt <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_brt = lapply(1:background_presence_ratio, FUN=function(i) {
    print(i);
    opt.LR <- TRUE;
    LR = 0.005;
    while(opt.LR){
      m.brt <- try(gbm.step(data = C_pipiens_occ_env[c(presences, which(C_pipiens_occ_env$abs_index == i)),], gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
      if (class(m.brt) == "try-error" | class(m.brt) == "NULL"){
        LR <- LR/2
      } else
        if(m.brt$gbm.call$best.trees<1000){
          LR <- LR/2
        } else 
          if(m.brt$gbm.call$best.trees>10000){
            LR <- LR*2
          } else { 
            opt.LR <- FALSE}}; 
    return(m.brt)})
  
  models_brt[[paste(my_preds, collapse = "+")]] <- m_brt
  
} # Close the loop over the three different predictor combinations



# Save the models
save(models_glm, models_gam, models_rf, models_brt, weights, predictors, my_preds_list, presences, C_pipiens_occ_env, background_presence_ratio,
     file = "output_data/models/C_pipiens_SDMs.RData")

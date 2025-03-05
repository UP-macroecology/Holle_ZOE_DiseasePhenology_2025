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
load("output_data/data/WNF_occ_env.RData")




#-------------------------------------------------------------------------------

# 1. Variable selection --------------------------------------------------------

# Create a subset to have a balanced presence-absence ratio per month and year in the
# data frame containing the occurrences
years <- c(2008:2019) # Create a vector containing the years
months <- str_pad(1:12, width = 2, pad = "0") # Create a vector containing the months

# Crate an empty data frame having the same columns as the occurrence data frame
WNF_occ_env_subset <- data.frame(matrix(ncol = 20, nrow = 0))
colnames(WNF_occ_env_subset) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "primary_forest", "primary_openland", "secondary_forest", 
                                  "secondary_openland", "pasture", "rangeland", "cropland", "urban", "C_pipiens")

for (y in years) { # Start the loop over the years
  
  print(y)
  
  for (m in months) { # Start the loop over the months
    
    print(m)
    
    # Create a subset of occurrence per month and year
    subset_month_year <- subset(WNF_occ_env, WNF_occ_env$year == y & WNF_occ_env$month == m)
    
    # Retain the rows containing presences
    presences <- subset_month_year[subset_month_year$occ == 1, ]
    
    # Get number of presences
    presence_numbers <- nrow(presences)
    
    # Randomly sample absence rows, equal number as presences
    absence_indices <- which(subset_month_year$occ == 0)
    indices_sampled <- sample(absence_indices, presence_numbers, replace = FALSE)
    absences_sampled <- subset_month_year[indices_sampled, ]
    
    # Bind presences and samples absences
    df_sampled <- rbind(presences, absences_sampled)
    
    # Add data to empty data frame
    WNF_occ_env_subset <- rbind(WNF_occ_env_subset, df_sampled)
    
  } # Close loop over months
} # Close loop over years

# Replace original occurrence data frame with subsetted data frame
WNF_occ_env <- WNF_occ_env_subset 


# Retrieve predictors (excluding predictors related to temperature, as we will include
# all of these in different models, because tmin, tmean, and tmax might be more relevant
# for different months of a year. We will include them in different models because they are
# highly correlated.)
predictors <- names(WNF_occ_env[, c(7, 11:20)])

# Check for collinearity in correlation matrix
cor_mat <- cor(WNF_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Generate weights
weights <- rep(1, times = nrow(WNF_occ_env))

# Run select07_cv function
var_sel <- select07_cv(X = WNF_occ_env[,predictors], 
                       y = WNF_occ_env$occ, 
                       threshold = 0.7,
                       weights = weights)

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
                    family='binomial', data = WNF_occ_env, weights = weights))
  
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
                     family='binomial', data = WNF_occ_env, weights = weights)
  
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
save(models_glm, models_gam, models_rf, models_brt, predictors, my_preds_list, WNF_occ_env, weights,
     file = "output_data/models/WNF_SDMs.RData")



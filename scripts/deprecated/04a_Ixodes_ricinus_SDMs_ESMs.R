# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

# ---------------------------------------------------------------------- #
#                          04a. Model fitting                            #
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
load("output_data/data/I_ricinus_occ_env.RData")



#-------------------------------------------------------------------------------

# 1. Bivariate variable selection ----------------------------------------------

# Retrieve predictors
predictors <- names(I_ricinus_occ_env[, c(7:19)])

# Define equal weights for presences and background data
weights <- ifelse(I_ricinus_occ_env$occ == 1, 1, sum(I_ricinus_occ_env$occ == 1) / sum(I_ricinus_occ_env$occ == 0))

# Machine-learning methods should be fitted with equal number of presences and background and repeated 10 times (Barbet-Massin et al. (2012))
presences <- which(I_ricinus_occ_env$occ == 1) # Extract the position index of presences
I_ricinus_occ_env$abs_index <- NA # Create a new column to insert the absence index
I_ricinus_occ_env$abs_index[I_ricinus_occ_env$occ!=1] <- sample(1:4, sum(I_ricinus_occ_env$occ!=1), replace = TRUE) # Insert sampled numbers of 1 to 10

# Find all possible bivariate combinations of predictor variables
predictors_combinations <- combn(predictors, 2, simplify = FALSE)

# Find out which variables are highly correlated (Spearman > 0.7)
cm <- cor(I_ricinus_occ_env[,predictors], method = "spearman")
pairs <- which(abs(cm)>= 0.7, arr.ind=T) # identifies correlated variable pairs
index <- which(pairs[,1] == pairs[,2]) # removes entry on diagonal (self-correlated)
pairs <- pairs[-index, ]

high_corr_pairs <- apply(pairs, 1, function(idx) { # Extract the variable names corresponding to the row and column indices
  c(rownames(cm)[idx["row"]], colnames(cm)[idx["col"]]) 
})


# Remove highly correlated variables within the list of bivariate combinations
# of predictor variables 
for (n in 1:ncol(high_corr_pairs)) {
  
  pair <- c(high_corr_pairs[1, n], high_corr_pairs[2, n])
  
  print(pair)
  
  predictors_combinations <- Filter(function(combination) {
    !all(combination == pair)
  }, predictors_combinations)
  
}



#-------------------------------------------------------------------------------

# 2. Model fitting of ESMs -----------------------------------------------------

# Fit individual GLMs (including linear and quadratic terms, AIC-based stepwise variable selection, equal weights)
# based on all bivariate combinations of predictors
print("GLM")

# Create a list to store the models
models_glm <- list()

for (c in seq_along(predictors_combinations)) { # Start of the loop over all bivariate predictor combinations
  
  # Get the respective combination of predictors
  my_preds <- predictors_combinations[[c]]
  
  print(my_preds)
  
  m_glm <- step(glm(as.formula(paste('occ~',paste(c(my_preds, paste0('I(',my_preds,'^2)')), collapse ='+'))),
                    family='binomial', data = I_ricinus_occ_env, weights = weights))
  
  models_glm[[c]] <- m_glm
  
} # Close the loop over all bivariate predictor combinations

names(models_glm) <- sapply(predictors_combinations, paste, collapse = "+")


# Fit individual GAMs (cubic smoothing splines, equal weights) based on all bivariate
# combinations of predictors
print("GAM")

models_gam <- list()

for (c in seq_along(predictors_combinations)) {
  
  my_preds <- predictors_combinations[[c]]
  
  print(my_preds)
  
  m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                     family='binomial', data = I_ricinus_occ_env, weights = weights)
  
  models_gam[[c]] <- m_gam
  
}

names(models_gam) <- sapply(predictors_combinations, paste, collapse = "+")




# Fit individual RF (same number of presences and background data, ten models in total)
# based on all bivariate combinations of predictors
print("RF")

models_rf <- list()

for (c in seq_along(predictors_combinations)) {
  
  my_preds <- predictors_combinations[[c]]
  
  print(my_preds)
  
  m_rf <- lapply(1:4, FUN=function(i){sp_train <- I_ricinus_occ_env[c(presences, which(I_ricinus_occ_env$abs_index == i)),]; print(i); randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                                                                                                                                                    data = sp_train, ntree = 1000, nodesize = 10, importance = T)})
 
   models_rf[[paste(my_preds, collapse = "+")]] <- m_rf
}



# Fit individual BRTs (same number of presences and background data, ten models in total, adaptable learning rate to fit model with 1000 and 10000 trees)
# based on all bivariate combinations of predictors
print("BRT")

models_brt <- list()

for (c in seq_along(predictors_combinations)) {
  
  my_preds <- predictors_combinations[[c]]
  
  print(my_preds)

  m_brt = lapply(1:4, FUN=function(i) {
    print(i);
    opt.LR <- TRUE;
    LR = 0.01;
    while(opt.LR){
      m.brt <- try(gbm.step(data = I_ricinus_occ_env[c(presences, which(I_ricinus_occ_env$abs_index == i)),], gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
      if (class(m.brt) == "try-error" | class(m.brt) == "NULL"){
        LR <- LR/2
      } else
        if(m.brt$gbm.call$best.trees<50){
          LR <- LR/2
        } else 
          if(m.brt$gbm.call$best.trees>10000){
            LR <- LR*2
          } else { 
            opt.LR <- FALSE}}; 
    return(m.brt)})
  
  models_brt[[paste(my_preds, collapse = "+")]] <- m_brt
  
  
}



# Save lists of models and additional data
save(models_glm, models_gam, models_rf, models_brt, weights, predictors_combinations, predictors, presences, I_ricinus_occ_env, 
     file = "output_data/models/I_ricinus_SDMs_ESMs.RData")

  

  
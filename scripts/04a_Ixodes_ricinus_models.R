# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                           04a. Model fitting                           #
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

# 1. Variable selection --------------------------------------------------------

# Retrieve predictors
predictors <- names(I_ricinus_occ_env[, c(7:19)])
#predictors <- setdiff(predictors, "hurs")

# Check for collinearity in correlation matrix
cor_mat <- cor(I_ricinus_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Define equal weights for presences and background data
weights <- ifelse(I_ricinus_occ_env$occ == 1, 1, sum(I_ricinus_occ_env$occ == 1) / sum(I_ricinus_occ_env$occ == 0))

# Run select07_cv function
var_sel <- select07_cv(X = I_ricinus_occ_env[,predictors], 
                       y = I_ricinus_occ_env$occ, 
                       threshold = 0.7,
                       weights = weights)

# Extract most important and weakly correlated predictors
my_preds <- var_sel$pred_sel

# Extract the ratio of presence to background data to know how many machine learning models we can build
ratio_presence_background <- round(sum(I_ricinus_occ_env$occ == 0) / sum(I_ricinus_occ_env$occ == 1), 0)

# Machine-learning methods should be fitted with equal number of presences and background and repeated 10 times (Barbet-Massin et al. (2012))
presences <- which(I_ricinus_occ_env$occ == 1) # Extract the position index of presences
I_ricinus_occ_env$abs_index <- NA # Create a new column to insert the absence index
I_ricinus_occ_env$abs_index[I_ricinus_occ_env$occ!=1] <- sample(1:ratio_presence_background, sum(I_ricinus_occ_env$occ!=1), replace = TRUE) # Insert sampled numbers of 1 to 10



#-------------------------------------------------------------------------------

# 2. Model fitting -------------------------------------------------------------

# Fit GLM (including linear and quadratic terms, AIC-based stepwise variable selection, equal weights)
print("GLM")
m_glm <- step(glm(as.formula(paste('occ~',paste(c(my_preds, paste0('I(',my_preds,'^2)')), collapse ='+'))),
                  family='binomial', data = I_ricinus_occ_env, weights = weights))


# Fit GAM (cubic smoothing splines, equal weights)
print("GAM")
m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                   family='binomial', data = I_ricinus_occ_env, weights = weights)



# Fit RF (same number of presences and background data, ten models in total)
print("RF")
m_rf <- lapply(1:ratio_presence_background, FUN=function(i){sp_train <- I_ricinus_occ_env[c(presences, which(I_ricinus_occ_env$abs_index == i)),]; print(i); randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                                                                                                                                                                          data = sp_train, ntree = 1000, nodesize = 10, importance = T)})

# Fit BRT (same number of presences and background data, ten models in total, adaptable learning rate to fit model with 1000 and 10000 trees)
print("BRT")
m_brt = lapply(1:ratio_presence_background, FUN=function(i) {
  print(i);
  opt.LR <- TRUE;
  LR = 0.01;
  while(opt.LR){
    m.brt <- try(gbm.step(data = I_ricinus_occ_env[c(presences, which(I_ricinus_occ_env$abs_index == i)),], gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
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


# Fit Maxent 
# print("Maxent")
# m_maxent <- lapply(1:10,FUN=function(i){sp_train <- I_ricinus_occ_env[c(presences, which(I_ricinus_occ_env$abs_index == i)),]; print(i); maxnet(p=sp_train$occ, data=sp_train[,my_preds])})



# Save the models
save(m_glm, m_gam, m_rf, m_brt, weights, predictors, my_preds, presences, I_ricinus_occ_env, ratio_presence_background, 
     file = "output_data/models/I_ricinus_SDMs.RData")
       

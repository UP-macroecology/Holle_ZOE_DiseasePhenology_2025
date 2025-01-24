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

# Retrieve predictors
predictors <- names(C_pipiens_occ_env[, c(7:19)])

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

# Machine-learning methods should be fitted with equal number of presences and background and repeated 10 times (Barbet-Massin et al. (2012))
# num_presences <- sum(C_pipiens_occ_env$occ == 1) # Count the number of presences
# presences <- which(C_pipiens_occ_env$occ == 1) # Extract the position index of presences
# 
# for (i in 1:10) {
#   # Create a column named "abs_index_i" and initialize with NA
#   C_pipiens_occ_env[[paste0("abs_index_", i)]] <- NA
#   
#   # Sample indices of absence rows (occ != 1)
#   absence_rows <- which(C_pipiens_occ_env$occ != 1)
#   sampled_rows <- sample(absence_rows, num_presences, replace = FALSE)
#   
#   # Assign the corresponding index value `i` to the sampled rows
#   C_pipiens_occ_env[[paste0("abs_index_", i)]][sampled_rows] <- i
# }

# Extract the ratio of presence to background data to know how many machine learning models we can build
ratio_presence_background <- round(sum(C_pipiens_occ_env$occ == 0) / sum(C_pipiens_occ_env$occ == 1), 0)

presences <- which(C_pipiens_occ_env$occ == 1) # Extract the position index of presences
C_pipiens_occ_env$abs_index <- NA # Create a new column to insert the absence index
C_pipiens_occ_env$abs_index[C_pipiens_occ_env$occ!=1] <- sample(1:ratio_presence_background, sum(C_pipiens_occ_env$occ!=1), replace = TRUE) # Insert sampled numbers of 1 to 10




#-------------------------------------------------------------------------------

# 2. Model fitting -------------------------------------------------------------

# Fit GLM (including linear and quadratic terms, AIC-based stepwise variable selection, equal weights)
print("GLM")
m_glm <- step(glm(as.formula(paste('occ~',paste(c(my_preds, paste0('I(',my_preds,'^2)')), collapse ='+'))),
                  family='binomial', data = C_pipiens_occ_env, weights = weights))


# Fit GAM (cubic smoothing splines, equal weights)
print("GAM")
m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                   family='binomial', data = C_pipiens_occ_env, weights = weights)

# Fit Maxent 
# print("Maxent")
# m_maxent <- lapply(1:10,FUN=function(i){sp_train <- Culex_occ_env[c(presences, which(Culex_occ_env$abs_index == i)),]; print(i); maxnet(p=sp_train$occ, data=sp_train[,my_preds])})


# Fit RF (same number of presences and background data, ten models in total)
print("RF")
m_rf <- lapply(1:ratio_presence_background, FUN=function(i){sp_train <- C_pipiens_occ_env[c(presences, which(C_pipiens_occ_env$abs_index == i)),]; print(i); randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                                                                                                                                                                    data = sp_train, ntree = 1000, nodesize = 10, importance = T)})

# Fit BRT (same number of presences and background data, ten models in total, adaptable learning rate to fit model with 1000 and 10000 trees)
print("BRT")
m_brt = lapply(1:ratio_presence_background, FUN=function(i) {
  print(i);
  opt.LR <- TRUE;
  LR = 0.008;
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



# Save the models
save(m_glm, m_gam, m_rf, m_brt, weights, predictors, my_preds, presences, C_pipiens_occ_env, ratio_presence_background,
     file = "output_data/models/C_pipiens_SDMs.RData")

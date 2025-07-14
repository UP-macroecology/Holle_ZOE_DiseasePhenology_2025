# ZOE project 
# Disease phenology analysis of TBE in Europe 

# ---------------------------------------------------------------------- #
#                          09a. Model fitting                            #
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
load("output_data/data/TBE_occ_env.RData")



#-------------------------------------------------------------------------------

# 1. Post-processing of absence generation -------------------------------------

# Create a subset to have a balanced presence-absence ratio per month and year in the
# data frame containing the occurrences
start_year <- min(TBE_occ_env$year) # Find the year of earliest observation
years <- c(start_year:2019) # Create a vector containing the years
months <- str_pad(1:12, width = 2, pad = "0") # Create a vector containing the months

# Crate an empty data frame having the same columns as the occurrence data frame
TBE_occ_env_subset <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(TBE_occ_env_subset) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "I_ricinus")

for (y in years) { # Start the loop over the years
  
  print(y)
  
  for (m in months) { # Start the loop over the months
    
    print(m)
    
    # Create a subset of occurrence per month and year
    subset_month_year <- subset(TBE_occ_env, TBE_occ_env$year == y & TBE_occ_env$month == m)
    
    # Retain the rows containing presences
    presences <- subset_month_year[subset_month_year$occ == 1, ]
    
    # Get number of presences
    presence_numbers <- nrow(presences)
    
    # Randomly sample absence rows, equal number as presences
    absence_indices <- which(subset_month_year$occ == 0)
    indices_sampled <- sample(absence_indices, presence_numbers, replace = FALSE)
    absences_sampled <- subset_month_year[indices_sampled, ]
    
    # Bind presences and sampled absences
    df_sampled <- rbind(presences, absences_sampled)
    
    # Add data to empty data frame
    TBE_occ_env_subset <- rbind(TBE_occ_env_subset, df_sampled)
    
  } # Close loop over months
} # Close loop over years

# Replace original occurrence data frame with subsetted data frame
TBE_occ_env <- TBE_occ_env_subset 

# Map the thinned presences and background data with balanced ratio
png("output_data/plots/presence_background/TBE_presence_absence.png", width = 2000, height = 2000, res = 300)


maps::map('world',xlim=c(-31,40), ylim=c(34,72))
points(TBE_occ_env$lon[TBE_occ_env$occ == 0], TBE_occ_env$lat[TBE_occ_env$occ == 0], col='steelblue4',  pch=19, cex = 0.5)
points(TBE_occ_env$lon[TBE_occ_env$occ == 1], TBE_occ_env$lat[TBE_occ_env$occ == 1], col='goldenrod',  pch=19, cex = 0.5)
legend(title = "TBE:", x = -25, y = 50, legend = c("Absence", "Presence"), col = c("steelblue4", "goldenrod"), pch = 19, pt.cex = 1, bty = "n")

dev.off()



#-------------------------------------------------------------------------------

# 2. Variable selection --------------------------------------------------------

# Retrieve predictors (excluding predictors related to temperature, as we will include
# all of these in different models, because tmin, tmean, and tmax might be more relevant
# for different months of a year. We will include them in different models because they are
# highly correlated.)
predictors <- names(TBE_occ_env[, c(7, 11, 12)])

# Check for collinearity in correlation matrix
cor_mat <- cor(TBE_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Generate weights
weights <- rep(1, times = nrow(TBE_occ_env))

# Run select07_cv function
var_sel <- select07_cv(X = TBE_occ_env[,predictors], 
                       y = TBE_occ_env$occ, 
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

# 3. Model fitting -------------------------------------------------------------

# Fit GLM (including linear and quadratic terms, AIC-based stepwise variable selection, but making sure 
# that the main vector species is included as linear term in the final model)
# for all three predictor sets
print("GLM")

# Create a list to store the models
models_glm <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  my_preds_no_vector <- setdiff(my_preds, "I_ricinus")

  # Define a model only containing Ixodes ricinus
  model_I_ricinus <- glm(as.formula('occ ~ I_ricinus'), family = 'binomial', data = TBE_occ_env, weights = weights)

  # Define a model with all predicotrs, including Ixodes ricinus
  model_full <- step(glm(as.formula(paste('occ ~ I_ricinus +', paste(c(my_preds_no_vector, paste0('I(', my_preds_no_vector, '^2)')), collapse = '+'))),
                    family='binomial', data = TBE_occ_env, weights = weights))


  # Perform step-wise selection, forcing Culex pipiens to be on of the predictors
  m_glm <- step(model_I_ricinus, scope = list(lower = model_I_ricinus, upper = model_full), direction = "both")


  models_glm[[m]] <- m_glm
  
} # Close the loop over the three different predictor combinations

names(models_glm) <- sapply(my_preds_list, paste, collapse = "+")



# Fit GAM (cubic smoothing splines) for all three predictor sets
print("GAM")

# Create a list to store the models
models_gam <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                     family='binomial', data = TBE_occ_env, weights = weights)
  
  models_gam[[m]] <- m_gam
  
} # Close the loop over the three different predictor combinations

names(models_gam) <- sapply(my_preds_list, paste, collapse = "+")



# Fit RF for all three predictor sets
print("RF")

models_rf <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_rf <- randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                       data = TBE_occ_env, ntree = 1000, nodesize = 10, importance = T)
  
  
  models_rf[[m]] <- m_rf
  
} # Close the loop over the three different predictor combinations

names(models_rf) <- sapply(my_preds_list, paste, collapse = "+")



# Fit BRT (adaptable learning rate to fit model between 1000 and 5000 trees)
# for all three predictor sets
print("BRT")

models_brt <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  opt.LR <- TRUE
  LR = 0.01
  
  while(opt.LR){
    m_brt <- try(gbm.step(data = TBE_occ_env, gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
    if (class(m_brt) == "try-error" | class(m_brt) == "NULL"){
      LR <- LR/2
    } else
      if(m_brt$gbm.call$best.trees<1000){
        LR <- LR/2
      } else 
        if(m_brt$gbm.call$best.trees>5000){
          LR <- LR*2
        } else { 
          opt.LR <- FALSE}}
  
  models_brt[[m]] <- m_brt
  
  
} # Close the loop over the three different predictor combinations

names(models_brt) <- sapply(my_preds_list, paste, collapse = "+")



# Save the models
save(models_glm, models_gam, models_rf, models_brt, predictors, my_preds_list, TBE_occ_env, weights,
     file = "output_data/models/TBE_SDMs.RData")




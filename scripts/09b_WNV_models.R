# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                         09b. Model fitting - WNV                       #
# ---------------------------------------------------------------------- #

# What is done within this script:

# In a first step, we create a balanced data set of thinned presences and absences.
# In a second step, we identify the most important and weakly correlated 
# predictor variables to include in model construction. Because the
# temperature variables (tas, tasmin, tasmax) are highly correlated, they are 
# excluded from the model selection process and later used to create three
# different predictor sets on which the final models are built. Models are built
# using four different algorithms: Generalised Linear Model (GLM), Generalised 
# Additive Model (GAM), Random Forest (RF), and Boosted Regression Tree (BRT).



# Load needed packages
library(mgcv) # mgcv_1.8-42
library(randomForest) # randomForest_4.7-1.1
library(gbm) # gbm_2.1.8.1
library(dismo) # dismo_1.3-14 
library(tidyverse) # tidyverse_2.0.0
library(corrplot) # corrplot_0.92
library(sf) # sf_1.0-16

# Load needed objects
source("scripts/00_functions.R") # Get the select07_cv function (explained deviance function)

# Read in presence and background data
load("output_data/data/WNV_occ_env.RData")




#-------------------------------------------------------------------------------

# 1. Post-processing of absence generation -------------------------------------

# Create a subset to have a balanced presence-absence ratio per month and year in the
# data frame containing the occurrences
start_year <- min(WNV_occ_env$year) # Find the year of earliest observation
years <- c(start_year:2019) # Create a vector containing the years
months <- str_pad(1:12, width = 2, pad = "0") # Create a vector containing the months

# Crate an empty data frame having the same columns as the occurrence data frame
WNV_occ_env_subset <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(WNV_occ_env_subset) <- c("lon", "lat", "occ", "year", "month", "ID", "pr", "tas", "tasmax", "tasmin", "hurs", "C_pipiens")

for (y in years) { # Start the loop over the years
  
  print(y)
  
  for (m in months) { # Start the loop over the months
    
    print(m)
    
    # Create a subset of occurrence per month and year
    subset_month_year <- WNV_occ_env[
      WNV_occ_env$year == y & 
        WNV_occ_env$month == m, 
    ]
    
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
    WNV_occ_env_subset <- rbind(WNV_occ_env_subset, df_sampled)
    
  } # Close loop over months
} # Close loop over years

# Replace original occurrence data frame with subsetted data frame
WNV_occ_env <- WNV_occ_env_subset 



#-------------------------------------------------------------------------------

# 2. Visualise thinned presences and absences ----------------------------------


# Map the thinned presences and absence data with balanced ratio
png("output_data/plots/presence_background/WNV_presence_absence.png", width = 2000, height = 2000, res = 300)


maps::map('world',xlim=c(-31,40), ylim=c(34,72), 
          col = "gray97",
          fill = TRUE,
          border = "gray30")

maps::map.axes(cex.axis = 0.75)


points(WNV_occ_env$lon[WNV_occ_env$occ == 0], WNV_occ_env$lat[WNV_occ_env$occ == 0], col='steelblue4',  pch=19, cex = 0.5)
points(WNV_occ_env$lon[WNV_occ_env$occ == 1], WNV_occ_env$lat[WNV_occ_env$occ == 1], col='goldenrod',  pch=19, cex = 0.5)

legend(title = "WNV:", x = -28, y = 50, legend = c("Absence", "Presence"), col = c("steelblue4", "goldenrod"), pch = 19, pt.cex = 1, bty = "n")

dev.off()



#-------------------------------------------------------------------------------

# 3. Variable selection --------------------------------------------------------

# Retrieve predictors (excluding predictors related to temperature, as we will include
# all of these in different models, because tmin, tmean, and tmax might be more relevant
# for different months of a year. We will include them in different models because they are
# highly correlated.)
predictors <- names(WNV_occ_env[, c(7, 11, 12)])

# Check for collinearity in correlation matrix
cor_mat <- cor(WNV_occ_env[,predictors], method='spearman')
corrplot.mixed(cor_mat, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

# Generate weights
weights <- rep(1, times = nrow(WNV_occ_env))

# Run select07_cv function
var_sel <- select07_cv(X = WNV_occ_env[,predictors], 
                       y = WNV_occ_env$occ, 
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

# 4. Model fitting -------------------------------------------------------------
# Fit models based on four different algorithms



# (a) Generalised linear models ------------------------------------------------

# Fit GLM (including linear and quadratic terms, AIC-based stepwise variable selection, but making sure 
# that the main vector species is included as linear term in the final model)
# for all three predictor sets
print("GLM")

# Create a list to store the models
models_glm <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  my_preds_no_vector <- setdiff(my_preds, "C_pipiens")

  # Define a model only containing Culex pipiens
  model_C_pipiens <- glm(as.formula('occ ~ C_pipiens'), family = 'binomial', data = WNV_occ_env, weights = weights)
  
  # Define a model with all predictors, including Culex pipiens
  model_full <- glm(as.formula(paste('occ ~ C_pipiens +', paste(c(my_preds_no_vector, paste0('I(', my_preds_no_vector, '^2)')), collapse = '+'))),
                    family='binomial', data = WNV_occ_env, weights = weights)

  # Perform step-wise selection, forcing Culex pipiens to be one of the predictors
  m_glm <- step(model_C_pipiens, scope = list(lower = model_C_pipiens, upper = model_full), direction = "both")


  models_glm[[m]] <- m_glm
  
} # Close the loop over the three different predictor combinations

names(models_glm) <- sapply(my_preds_list, paste, collapse = "+")

# par(mfrow=c(2,2))
# partial_response(models_glm[[1]], predictors = WNV_occ_env[,c("hurs", "C_pipiens", "pr", "tas")], ylab='Occurrence probability')

# (b) Generalised additive models ----------------------------------------------

# Fit GAM (cubic smoothing splines) for all three predictor sets
print("GAM")

# Create a list to store the models
models_gam <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_gam <- mgcv::gam(as.formula(paste('occ~',paste(paste0('s(',my_preds,',k=4)'), collapse='+'))),
                     family='binomial', data = WNV_occ_env, weights = weights)
  
  models_gam[[m]] <- m_gam
  
} # Close the loop over the three different predictor combinations

names(models_gam) <- sapply(my_preds_list, paste, collapse = "+")



# (c) Random forests -----------------------------------------------------------

# Fit RF for all three predictor sets
print("RF")

models_rf <- list()

for (m in seq_along(my_preds_list)) { # Start the loop over the three different predictor combinations
  
  my_preds <- my_preds_list[[m]] # Retain the predictors
  print(my_preds)
  
  m_rf <- randomForest(as.formula(paste('occ~',paste(my_preds, collapse='+'))), 
                       data = WNV_occ_env, ntree = 1000, nodesize = 10, importance = T)
  
  
  models_rf[[m]] <- m_rf
  
} # Close the loop over the three different predictor combinations

names(models_rf) <- sapply(my_preds_list, paste, collapse = "+")



# (d) Boosted regression trees -------------------------------------------------

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
    m_brt <- try(gbm.step(data = WNV_occ_env, gbm.x = my_preds, gbm.y = "occ", family = 'bernoulli', tree.complexity = 2, bag.fraction = 0.75, learning.rate = LR, verbose=F, plot.main=F))
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



# (e) # Save the models --------------------------------------------------------

save(models_glm, models_gam, models_rf, models_brt, predictors, my_preds_list, WNV_occ_env, weights,
     file = "output_data/models/WNV_SDMs.RData")





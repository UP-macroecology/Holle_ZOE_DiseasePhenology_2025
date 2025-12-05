# ZOE project 


#-------------------------------------------------------------------------------

# ---------------------------------------------------------------------- #
#                         10c. Response curve plot                       #
# ---------------------------------------------------------------------- #

# What is done within this script:

# We plot the partial response curves for Ixodes ricinus, Culex pipiens, and 
# their associated diseases (TBE and WNV) in a single combined plot



# Load needed packages
library(ggplot2) # ggplot2_4.0.0

# Load needed data
load("output_data/validation/C_pipiens_response_data.RData") # Predictor response data from Culex pipiens
load("output_data/validation/I_ricinus_response_data.RData") # Predictor response data from Ixodes ricinus
load("output_data/validation/WNV_response_data.RData") # Predictor response data from WNV
load("output_data/validation/TBE_response_data.RData") # Predictor response data from TBE



#-------------------------------------------------------------------------------

# 1. Visualise partial response plots  -----------------------------------------

# Including Ixodes ricinus, TBE, Culex pipiens, and WNV

# Bind the predictor response data from each species into one data frame
response_data <- rbind(C_pipiens_response_data, I_ricinus_response_data, WNV_response_data, TBE_response_data)

# Arrange order of predictors to group them in land use and climate ones
response_data$predictor <- factor(response_data$predictor, levels = c("cropland", "pasture", "primary_forest",
                                                                      "primary_openland", "rangeland", "secondary_forest", 
                                                                      "secondary_openland","urban", 
                                                                      "hurs", "pr", "tas", "tasmax", "tasmin",
                                                                      "C_pipiens", "I_ricinus"))

# Create a named vector to rename predictors
my_labels <- c("cropland" = "Cropland", 
               "hurs" = "Relative humidity", 
               "pasture" = "Pasture",
               "pr" = "Precipitation", 
               "primary_forest" = "Primary forest", 
               "primary_openland" = "Primary open land",
               "rangeland" = "Rangeland", 
               "secondary_forest" = "Secondary forest",
               "secondary_openland" = "Secondary open land", 
               "tas" = "Mean temperature",
               "tasmax" = "Max temperature", 
               "tasmin" = "Min temperature", 
               "urban" = "Urban",
               "I_ricinus" = "Ixodes ricinus",
               "C_pipiens" = "Culex pipiens")

# Arrange order of disease vector and pathogens
response_data$species <- factor(response_data$species, levels = c("Ixodes ricinus", "TBE",
                                                                  "Culex pipiens", "WNV"))

# Plot the response curves of the different predictors for the different
# disease vectors and pathogens
ggplot(data = response_data, aes(x = environmental_values, y = predicted_values, color = species)) +
  geom_line(linewidth = 0.7, alpha = 0.75) + 
  facet_wrap(~ predictor, scales = "free_x", labeller = labeller(predictor = my_labels)) + 
  labs(x = "Environmental Values", y = "Predicted Values", color = "Disease vectors and pathogens") + 
  ylim(0,1) +
  theme_bw() +
  scale_color_manual(values = c("Ixodes ricinus" = "springgreen4", "Culex pipiens" = "royalblue4",
                                "TBE" = "darkseagreen3", "WNV" = "lightsteelblue3")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 7.5),
        axis.title = element_text(size = 7.5), 
        axis.text = element_text(size = 4),
        strip.text = element_text(size = 6.5, face = "bold"),
        strip.background = element_rect(fill = "grey75", color = NA)) +
  guides(color = guide_legend(title.position = "top", nrow = 1, byrow = TRUE))

# Save the plot containing the repsonse curves of the different vectors and
# pathogens for each of their predictors
ggsave("output_data/plots/response_curves/partial_response_curves_all.png", width = 15, height = 14, units = "cm")




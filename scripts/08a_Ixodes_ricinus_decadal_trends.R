# ZOE project 
# Disease phenology analysis of Ixodes ricinus in Europe (primary transmitter of TBEV)

# ---------------------------------------------------------------------- #
#       07a. Past and future decadal trends of disease phenology         #
# ---------------------------------------------------------------------- #



# Load needed packages
library(ggplot2)
library(terra)
library(tidyverse)



#-------------------------------------------------------------------------------

# 1. Past disease phenology ----------------------------------------------------
# Calculate the past decadal trends of peak and mean occurrence probability per month

# Read in past monthly prediction data from 1970 to 2019 
r_curr_preds_clim_landuse <- terra::rast("output_data/results/I_ricinus_preds_clim_landuse_1970_2019.tif") # under observed climate and land use change
r_curr_preds_noclim_landuse <- terra::rast("output_data/results/I_ricinus_preds_noclim_landuse_1970_2019.tif") # under observed land use change and counterfactual climate 

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("max", "mean")


for (o in operations) { # Loop over the max and mean functions 
  
  print(o)
  
  # Generate time information
  dates <- seq(as.Date("1970-01-01"), as.Date("2019-12-01"), by = "month")
  
  # Assign dates as layer names
  names(r_curr_preds_clim_landuse) <- dates
  names(r_curr_preds_noclim_landuse) <- dates
  
  # Calculate mean or maximum (peak) occurrence probability across the study area (Europe) for each month
  # and create a data frame with the calculated information (under the factual and counterfactual scenario)
  monthly_clim_landuse <- global(r_curr_preds_clim_landuse, fun = o, na.rm = TRUE)
  monthly_clim_landuse <- data.frame(date = dates, occurrence = monthly_clim_landuse[,1], scenario = "Factual")
  
  monthly_noclim_landuse <- global(r_curr_preds_noclim_landuse, fun = o, na.rm = TRUE)
  monthly_noclim_landuse <- data.frame(date = dates, occurrence = monthly_noclim_landuse[,1], scenario = "Counterfactual climate")
  
  
  # Combine the data frames
  combined_df <- bind_rows(monthly_clim_landuse, monthly_noclim_landuse)
  
  # Create columns for month (1:12) and decade (1970s, 1980s, etc.) to group data
  combined_df <- combined_df %>%
    mutate(
      month = month(date),
      decade = case_when(
        year(date) >= 1970 & year(date) < 1980 ~ "1970s",
        year(date) >= 1980 & year(date) < 1990 ~ "1980s",
        year(date) >= 1990 & year(date) < 2000 ~ "1990s",
        year(date) >= 2000 & year(date) < 2010 ~ "2000s",
        year(date) >= 2010 & year(date) < 2020 ~ "2010s"
      )
    )
  
  # Aggregate by decade and month to calculate the mean for each time window (decade) and month
  aggregated_df <- combined_df %>%
    group_by(decade, month, scenario) %>%
    summarise(
      occurrence = mean(occurrence, na.rm = TRUE),
    ) %>%
    ungroup()
  
  # Change month numbers to month name abbreviations
  aggregated_df <- aggregated_df %>%
    mutate(month = factor(month, levels = 1:12, labels = month.abb))
  
  # Adapt ylim for the plot based on used function
  if (o == "max") { y_lim_values <- c(0.745, 0.845) 
  } else if (o == "mean") {y_lim_values <- c(0.30, 0.435)}
  
  # Visualize the data
  ggplot(aggregated_df, aes(x = month, y = occurrence, group = interaction(decade, scenario), color = decade)) +
    geom_line(aes(linetype = scenario), size = 1) + 
    labs(
      x = "month in a year",
      y = paste(o, "occurrence probability"),
      title = paste("Ixodes ricinus - Decadal trends in monthly", o, "occurrence probability in Europe"),
      color = "Decade",
      linetype = "Scenario"
    ) +
    scale_x_discrete(labels = month.abb) + 
    ylim(y_lim_values) +
    guides(fill = "none") +
    guides(
      color = guide_legend(nrow = 3, title.position = "top"),
      linetype = guide_legend(nrow = 2, title.position = "top") 
    ) +
    scale_linetype_manual(
      values = c("Factual" = "solid", "Counterfactual climate" = "twodash")) +
    scale_color_manual(
      values = c(
        "1970s" = "black",
        "1980s" = "midnightblue",
        "1990s" = "cadetblue",
        "2000s" = "seagreen",
        "2010s" = "aquamarine3"
      )) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      text = element_text(size = 14),  
      axis.title = element_text(size = 15),  
      axis.text = element_text(size = 12),  
      legend.title = element_text(size = 14),  
      legend.text = element_text(size = 12),  
      plot.title = element_text(size = 15, face = "bold"))
  
  
  ggsave(paste0("output_data/plots/decadal_trends/I_ricinus_decadal_trends_",o,"occprob_past.png"), width = 9, height = 7.5)
  
} # Close the loop over max and mean function







#-------------------------------------------------------------------------------

# 2. Future disease phenology --------------------------------------------------
# Calculate the future decadal trends of peak and mean occurrence probability per month
# for three different environmental change scenarios

# Create a vector containing the three different environmental change scenarios
env_scenarios <- c("ssp126", "ssp370", "ssp585")

# Create a vector indicating the mathematical operations to extract peak and mean 
# occurrence probabilities
operations <- c("max", "mean")


for (s in env_scenarios) { # Start of the loop over respective environmental change scenario
  
  print(s)

  # Read in future monthly prediction data from 2030 to 2070
  r_fut_preds_clim_landuse <- terra::rast(paste0("output_data/results/I_ricinus_preds_clim_landuse_2030_2070_",s,".tif")) # under scenario of future climate and land use change
  r_fut_preds_clim_nolanduse <- terra::rast(paste0("output_data/results/I_ricinus_preds_clim_nolanduse_2030_2070_",s,".tif")) # under scenario of future climate change and steady land use
  
  
  # Generate time information
  dates <- seq(as.Date("2030-01-01"), as.Date("2070-12-01"), by = "month")
  
  # Assign dates as layer names
  names(r_fut_preds_clim_landuse) <- dates
  names(r_fut_preds_clim_nolanduse) <- dates
  
  for (o in operations) { # Loop over the max and mean functions 
    
    print(o)
    
    monthly_clim_landuse <- global(r_fut_preds_clim_landuse, fun = o, na.rm = TRUE)
    monthly_clim_landuse <- data.frame(date = dates, occurrence = monthly_clim_landuse[,1], scenario = "Climate + Land use change")
    
    monthly_clim_nolanduse <- global(r_fut_preds_clim_nolanduse, fun = o, na.rm = TRUE)
    monthly_clim_nolanduse <- data.frame(date = dates, occurrence = monthly_clim_nolanduse[,1], scenario = "Climate change")
    
    # Combine the data frames
    combined_df <- bind_rows(monthly_clim_landuse, monthly_clim_nolanduse)
    
  
    # Create columns for month (1:12) and decade (2030s, 2040s, etc.) to group data
    combined_df <- combined_df %>%
      mutate(
        month = month(date),
        decade = case_when(
          year(date) >= 2030 & year(date) < 2040 ~ "2030s",
          year(date) >= 2040 & year(date) < 2050 ~ "2040s",
          year(date) >= 2050 & year(date) < 2060 ~ "2050s",
          year(date) >= 2060 & year(date) < 2070 ~ "2060s",
        )
      )
    
    # Remove rows that conatin NA values
    combined_df <- na.omit(combined_df)
    
    # Aggregate by decade and month to calculate the mean for each time window (decade) and month
    aggregated_df <- combined_df %>%
      group_by(decade, month, scenario) %>%
      summarise(
        occurrence = mean(occurrence, na.rm = TRUE),
      ) %>%
      ungroup()
    
    # Change month numbers to month name abbreviations
    aggregated_df <- aggregated_df %>%
      mutate(month = factor(month, levels = 1:12, labels = month.abb))
    
    # Adapt ylim for the plot based on used function
    if (o == "max") { y_lim_values <- c(0.745, 0.845) 
    } else if (o == "mean") {y_lim_values <- c(0.30, 0.47)}
    
    # Visualize the data for respective scenario
    ggplot(aggregated_df, aes(x = month, y = occurrence, group = interaction(decade, scenario), color = decade)) +
      geom_line(aes(linetype = scenario), size = 1.25) + 
      labs(
        x = "month in a year",
        y = paste(o, "occurrence probability"),
        title = paste("Environmental forcing scenario", s),
        color = "Future decade",
        linetype = "Scenario"
      ) +
      scale_x_discrete(labels = month.abb) +
      ylim(y_lim_values) +
      guides(fill = "none") +
      scale_linetype_manual(
        values = c("Climate change" = "longdash", "Climate + Land use change" = "dotted")) +
      scale_color_manual(
        values = c(
          "2030s" = "#F7D13D",
          "2040s" = "#FB9B06",
          "2050s" = "#CF4446",
          "2060s" = "#A52C60"
        )) +
      theme_minimal() +
      theme(text = element_text(size = 14),  
            axis.title = element_text(size = 15),  
            axis.text = element_text(size = 12),  
            legend.title = element_text(size = 14),  
            legend.text = element_text(size = 12),  
            plot.title = element_text(size = 13))
    
    ggsave(paste0("output_data/plots/decadal_trends/I_ricinus_decadal_trends_",o,"occprob_future_",s,".png"), width = 10, height = 6)
    
  } # Close the loop over the mean and max functions
  
} # Close the loop over the three environmental change scenarios
    
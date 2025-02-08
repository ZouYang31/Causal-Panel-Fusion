#############################
# Data Cleaning for COVID Data
# Preprocessing & Transformation
#############################


# Load Required Libraries
library(dplyr)
library(ggplot2)
library(limSolve)
library(tidyr)
library(kernlab)

#Load Data
covid_demo_yr <- read.csv("/Users/leeseunghee/Dropbox/Vaccine_Causal/Data/final_transformed_data_20cities_v2.csv")
covid_vac_mo <- read.csv("/Users/leeseunghee/Dropbox/Vaccine_Causal/Data/monthly_covid_vaccination_rate_20cities_v2.csv")
# Filter data for relevant timeframe
covid_demo_yr <- covid_demo_yr %>% filter(year >= 2021)
covid_vac_mo <- covid_vac_mo %>% filter(year_month < "2022-07")

#For the future usage
covid_demo_copy <- covid_demo_yr

# -----------------------------
# Data Preprocessing Functions
# -----------------------------

# Min-Max Normalization Function
 min_max_norm <- function(x) {
  if (is.numeric(x)) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  } else {
    return(x)
  }
}

# Function to Reorder Data, Prioritizing a Specific Unit
reorder_unit_first <- function(data, unit_pattern) {
  data %>%
    filter(grepl(unit_pattern, unit)) %>%
    bind_rows(filter(data, !grepl(unit_pattern, unit)))
}

# Function to Normalize Selected Columns (Ensures Data Frame Format)
normalize_columns <- function(data, columns) {
  data[columns] <- as.data.frame(lapply(data[columns], min_max_norm))
  return(data)
}

# Helper function to combine yearly data for F Z X 
combine_yearly_data <- function(data, year1, year2, columns) {
  combined <- cbind(
    data[data$year == year1, columns, drop = FALSE],
    data[data$year == year2, columns, drop = FALSE]
  )
  colnames(combined) <- NULL  # Remove column names
  return(combined)
}

# Helper function to combine yearly data for table
combine_yearly_col_data <- function(data, year1, year2, columns) {
  
  # Filter rows for the specified years
  data_year1 <- data %>% filter(year == year1) %>% select(all_of(columns))
  data_year2 <- data %>% filter(year == year2) %>% select(all_of(columns))
  
  # Ensure columns are named uniquely for each year
  colnames(data_year1) <- paste0(columns, "_", year1)
  colnames(data_year2) <- paste0(columns, "_", year2)
  
  # Combine datasets
  combined <- cbind(data_year1, data_year2)
  
  return(combined)
}

# Function: Aggregate Data by Unit (Average Over Years)
compute_average_by_unit <- function(data,covariate_columns) {
  # Select numeric columns excluding 'year'
  numeric_cols <- setdiff(names(data), c("unit", "year"))
  
  averages_by_unit <- aggregate(
    data[, covariate_columns], 
    by = list(unit = data$unit), 
    FUN = mean, na.rm = TRUE
  )
  return(averages_by_unit)
}


compute_synth_avg <- function(data, w) {
  
  # Compute synthetic weighted averages
  synth_median_income   <- sum(data[, "median_income"] * w)
  synth_proportion      <- sum(data[, "proportion"] * w)
  synth_median_age      <- sum(data[, "median_age"] * w)
  synth_65proportion    <- sum(data[, "X65p_proportion"] * w)
  
  # Store results in a named data frame for easy access
  synth_results <- data.frame(
    median_income  = synth_median_income,
    proportion     = synth_proportion,
    median_age     = synth_median_age,
    X65p_proportion = synth_65proportion
  )
  
  return(synth_results)
}


compute_controls_avg <- function(data){
  # Compute averages for each covariate
  avg_median_income   <- mean(data[, "median_income"], na.rm = TRUE)
  avg_proportion      <- mean(data[, "proportion"], na.rm = TRUE)
  avg_median_age      <- mean(data[, "median_age"], na.rm = TRUE)
  avg_X65p_proportion <- mean(data[, "X65p_proportion"], na.rm = TRUE)
  
  # Store results in a named data frame
  controls_avg_df <- data.frame(
    median_income   = avg_median_income,
    proportion      = avg_proportion,
    median_age      = avg_median_age,
    X65p_proportion = avg_X65p_proportion
  )
  
  # Convert data frame to table format for display
  controls_avg_table <- knitr::kable(controls_avg_df, format = "simple")
  
  # Return both data frame and table
  return(list(
    controls_avg_df = controls_avg_df))
  
}

# ------------------------------------------------------------------------------#
# Function: Process Data for Synthetic, Real, and Control Groups
# ------------------------------------------------------------------------------#
process_group_data <- function(df_filtered, unit_name, w) {
  
  # Compute average data for the filtered dataset
  avg_data <- compute_average_by_unit(df_filtered, covariate_columns)
  
  # Compute synthetic data using weights
  synth_data <- compute_synth_avg(avg_data, w)
  
  # Extract Chelsea-specific row and compute its average
  chelsea_row <- covid_demo_copy[covid_demo_copy$unit == unit_name, ]
  avg_chelsea <- compute_average_by_unit(chelsea_row, covariate_columns)
  
  # Compute control group averages
  avg_controls <- compute_controls_avg(avg_data)
  
  return(list(
    avg_data = avg_data,
    synth_data = synth_data,
    avg_chelsea = avg_chelsea,
    avg_controls = avg_controls
  ))
}


# Convert Data Frames to Numeric Matrices
convert_to_numeric_matrix <- function(data) {
  data %>%
    select(-unit) %>%  # Ensure only numeric columns are processed
    as.matrix() %>%
    apply(2, as.numeric)
}


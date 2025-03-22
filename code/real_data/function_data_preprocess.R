#############################
# Data Cleaning for COVID Data
# Preprocessing & Transformation
#############################

# 
# #Load Data
# covid_demo_yr <- read.csv(paste0(dir$data, "/final_transformed_data_20cities_v2.csv"))
# covid_vac_mo <- read.csv(paste0(dir$data, "/monthly_covid_vaccination_rate_20cities_v2.csv"))
# # Filter data for relevant timeframe
# covid_demo_yr <- covid_demo_yr %>% filter(year >= 2021)
# covid_vac_mo <- covid_vac_mo %>% filter(year_month < "2022-07")
# 
# #For the future usage
# covid_demo_copy <- covid_demo_yr

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

# Function to Normalize Selected Columns (Ensures Data Frame Format)
normalize_columns <- function(data, columns) {
  data[columns] <- as.data.frame(lapply(data[columns], min_max_norm))
  return(data)
}


# Function to Reorder Data, Prioritizing a Specific Unit
reorder_unit_first <- function(data, unit_pattern) {
  data %>%
    filter(grepl(unit_pattern, unit)) %>%
    bind_rows(filter(data, !grepl(unit_pattern, unit)))
}


compute_synth_avg <- function(data, w) {
  
  # Compute synthetic weighted averages
  synth_median_income   <- sum(data[data$year == 2021, "median_income"][-1] * w)
  synth_proportion      <- sum(data[data$year == 2021, "proportion"][-1] * w)
  synth_median_age      <- sum(data[data$year == 2021, "median_age"][-1] * w)
  synth_65proportion    <- sum(data[data$year == 2021, "X65p_proportion"][-1] * w)
  
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
  avg_median_income   <- mean(data[data$year == 2021, "median_income"], na.rm = TRUE)
  avg_proportion      <- mean(data[data$year == 2021, "proportion"], na.rm = TRUE)
  avg_median_age      <- mean(data[data$year == 2021, "median_age"], na.rm = TRUE)
  avg_X65p_proportion <- mean(data[data$year == 2021, "X65p_proportion"], na.rm = TRUE)
  
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
  if(unit_name == "Chelsea-black"){
    chelsea_row <- df_filtered[df_filtered$unit == unit_name, ] 
  }
  else {
    chelsea_row <- df_filtered[df_filtered$unit == unit_name, ] 
  }
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
convert_to_matrix <- function(data) {
  data %>%
    select(-unit) %>%  # Ensure only numeric columns are processed
    as.matrix()
}

# Generate the B list
generate_b_list <- function(step = 0.01, min_value = 0.01) {
  B <- list()  # Initialize empty list
  c <- 1  # Counter
  
  for (b_F in seq(min_value, 1 - min_value, by = step)) {
    for (b_Z in seq(min_value, 1 - b_F, by = step)) {
      b_X <- 1 - b_F - b_Z  # Ensure sum to 1
      
      if (b_X >= min_value) {  # Ensure b_X meets the minimum threshold
        B[[c]] <- list(b_F = b_F, b_Z = b_Z, b_X = b_X)
        c <- c + 1
      }
    }
  }
  
  return(B)
}


## ---------------------------------------  ##
## ----         placebo test -------------  ##

Placebo_test_data <- function(F, Y, J, t_max, s_max, i_max, 
                         w, X, Z, dr, dt,
                         eta_Z = 0.1, eta_X = 0.1){ #X, Z, dr, dt
  
  # Select all control units as placebo treated units
  placebo_units <- 2:i_max  
  
  synthetic_treated_F <- t(w) %*% F[2:(J+1), ]
  synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]
  
  # Initialize lists to store synthetic placebo outcomes
  synthetic_placebo_list_F <- list()
  synthetic_placebo_list_Y <- list()
  
  F_treated <- F[1, ]
  F_control <- F[2:(J+1), ]
  # Calculate the baseline X and Z
  result_X <- optimize_w_ipop(F_treated = F_treated, F_control = F_control, 
                              X, Z, t_max, dr, dt, i_max, target = "X")
  wX <- result_X$weights
  NSE_X_baseline <- NSE_x(wX, F, X, Z, t_max, dr, dt, i_max, target = "X")
  print(NSE_X_baseline)
  
  result_Z <- optimize_w_ipop(F_treated = F_treated, F_control = F_control, 
                              X, Z, t_max, dr, dt, i_max, target = "Z")
  wZ <- result_Z$weights
  NSE_Z_baseline <- NSE_x(wZ, F, X, Z, t_max, dr, dt, i_max, target = "Z")
  print(NSE_Z_baseline)
  
  for (placebo_unit in placebo_units) {
    print(placebo_unit)
    F_placebo_treated <- F[placebo_unit, ]
    F_placebo_control <- F[-placebo_unit, ]
    
    Y_placebo_treated <- Y[placebo_unit, ]
    Y_placebo_control <- Y[-placebo_unit, ]
    
    X_placebo_treated <- X[placebo_unit, ]
    X_placebo_control <- X[-placebo_unit, ]
    
    Z_placebo_treated <- Z[placebo_unit, ]
    Z_placebo_control <- Z[-placebo_unit, ]
    
    b_list <- find_best_B(B, F_treated = F_placebo_treated, F_control = F_placebo_control,
                          X_treated = X_placebo_treated, X_control = X_placebo_control,
                          Z_treated = Z_placebo_treated, Z_control = Z_placebo_control,
                          t_max, dr, dt, i_max, 
                          NSE_Z_baseline = NSE_Z_baseline, NSE_X_baseline = NSE_X_baseline, eta_Z= eta_Z, eta_X= eta_X) 
    # Solve for W
    w_placebo <- b_list$best_w_star
    print(w_placebo)
    
    if (is.null(w_placebo)) {
      # Define default behavior when w_placebo is NULL
      warning(paste("w_placebo is NULL for unit", placebo_unit, ". Skipping this unit."))
      next # Skip to the next iteration of the loop
    }
    
    synthetic_placebo_treated_F <- t(w_placebo) %*% F[-placebo_unit, ]
    effect_F <- F[placebo_unit, ] - synthetic_placebo_treated_F
    
    synthetic_placebo_list_F[[paste0("Unit_", placebo_unit)]] <- effect_F
    
    synthetic_placebo_treated_Y <- t(w_placebo) %*% Y[-placebo_unit, ]
    
    effect_Y <- Y[placebo_unit, ] - synthetic_placebo_treated_Y
    synthetic_placebo_list_Y[[paste0("Unit_", placebo_unit)]] <- effect_Y
  }
  
  # Combine the synthetic placebo outcomes into data frames
  synthetic_placebo_df_F <- do.call(rbind, lapply(names(synthetic_placebo_list_F), function(unit) {
    data.frame(
      year = 1:t_max,
      Outcome = as.numeric(synthetic_placebo_list_F[[unit]]),
      Unit = unit
    )
  }))
  
  synthetic_placebo_df_Y <- do.call(rbind, lapply(names(synthetic_placebo_list_Y), function(unit) {
    data.frame(
      year = 1:s_max,
      Outcome = as.numeric(synthetic_placebo_list_Y[[unit]]),
      Unit = unit
    )
  }))
  
  # Highlight the treated unit (Unit_1)
  synthetic_placebo_df_F <- rbind(synthetic_placebo_df_F, data.frame(
    year = 1:t_max,
    Outcome = as.numeric(F[1,] - synthetic_treated_F),
    Unit = "Synthetic_Treated_F"
  ))
  
  synthetic_placebo_df_Y <- rbind(synthetic_placebo_df_Y, data.frame(
    year = 1:s_max,
    Outcome = as.numeric(Y[1,] - synthetic_treated_Y),
    Unit = "Synthetic_Treated_Y"
  ))
  
  color_map_F <- c("Synthetic_Treated_F" = "black")
  line_type_map_F <- c("Synthetic_Treated_F" = "solid")
  
  for (unit in unique(synthetic_placebo_df_F$Unit)) {
    if (!unit %in% names(color_map_F)) {
      color_map_F[unit] <- scales::alpha("grey", 0.6) 
      line_type_map_F[unit] <- "solid"
    }
  }
  
  color_map_Y <- c("Synthetic_Treated_Y" = "black")
  line_type_map_Y <- c("Synthetic_Treated_Y" = "solid")
  
  for (unit in unique(synthetic_placebo_df_Y$Unit)) {
    if (!unit %in% names(color_map_Y)) {
      color_map_Y[unit] <- scales::alpha("grey", 0.6)
      line_type_map_Y[unit] <- "solid"
    }
  }
  
  # Custom legend breaks and labels
  legend_breaks <- c("Synthetic_Treated_F", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("Chelsea", "19 Cities")
  
  # Create the plots
  plot_F <- ggplot(synthetic_placebo_df_F, aes(x = year, y = Outcome, color = Unit, linetype = Unit)) +
    geom_line(data = subset(synthetic_placebo_df_F, !Unit %in% c("Synthetic_Treated_F")), 
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.6) +  # Grey control lines first
    geom_line(data = subset(synthetic_placebo_df_F, Unit == "Synthetic_Treated_F"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.5) + # Line at y = 0
    scale_color_manual(values = color_map_F, 
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map_F, 
                          breaks = legend_breaks,
                          labels = legend_labels) +
    scale_y_continuous(limits = c(-0.1, 0.1)) +
    labs(title = NULL, x = "Time", y = "Causal Effect", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.1), 
          legend.justification = c(1, 1), 
          legend.text = element_text(size = 14), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "solid"), 
                                                                  color = c("black", "grey")),
                                order = 1),
           linetype = guide_legend(title = NULL, override.aes = list(color = c("black", "grey")), 
                                   order = 1))
  
  ggsave(filename = paste0(dir$output, "/placebo_ref.pdf"), width = 8, height = 6, dpi = 300)
  
  legend_breaks <- c("Synthetic_Treated_Y", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("Chelsea", "19 Cities")
  
  plot_Y <- ggplot(synthetic_placebo_df_Y, aes(x = year, y = Outcome, color = Unit, linetype = Unit)) +
    geom_line(data = subset(synthetic_placebo_df_Y, !Unit %in% c("Synthetic_Treated_Y")), 
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.6) +  # Grey control lines first
    geom_line(data = subset(synthetic_placebo_df_Y, Unit == "Synthetic_Treated_Y"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 0.5) + # Line at y = 0
    scale_color_manual(values = color_map_Y, 
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map_Y, 
                          breaks = legend_breaks,
                          labels = legend_labels) +
    labs(title = NULL, x = "Time", y = "Causal Effect", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.1), 
          legend.justification = c(1, 1), 
          legend.text = element_text(size = 14), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "solid"), 
                                                                  color = c("black", "grey")),
                                order = 1),
           linetype = guide_legend(title = NULL, override.aes = list(color = c("black", "grey")), 
                                   order = 1))
  
  ggsave(filename = paste0(dir$output, "/placebo_target.pdf"), width = 8, height = 6, dpi = 300)
  print(plot_F)
  print(plot_Y)
}


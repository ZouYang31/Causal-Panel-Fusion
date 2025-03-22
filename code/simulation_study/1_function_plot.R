##########################################
# Data Fusion Simulation Study function
# Author: Zoey Yang
# Date: 
#########################################

synthetic_reference_plot <- function(F, w, t_max, i_max){
  
  F_treated <-  F[1, ]
  synthetic_treated_F <- t(w) %*% F[2:(J+1), ]
  
  # Create a data frame for plotting
  data <- data.frame(
    year = 1:t_max,
    Real_Treated = F[1, ],
    Synthetic_Treated = t(synthetic_treated_F)
  )
  
  # Create a data frame for plotting all units in F
  data_all <- data.frame(
    year = rep(1:t_max, i_max),
    Outcome = as.vector(t(F)),
    Unit = rep(paste0("Unit_", 1:i_max), each = t_max)
  )
  
  # Add the synthetic treated unit to the data frame
  data_synthetic <- data.frame(
    year = 1:t_max,
    Outcome = as.vector(synthetic_treated_F),
    Unit = rep("Synthetic_control", t_max)
  )
  
  # Combine the data frames
  data_all <- rbind(data_all, data_synthetic)
  
  # Define colors and line types
  color_map <- c("Unit_1" = "black", "Synthetic_control" = "black")
  line_type_map <- c("Unit_1" = "solid", "Synthetic_control" = "dashed")
  
  # Set the color and linetype for other units
  for (unit in unique(data_all$Unit)) {
    if (!unit %in% names(color_map)) {
      color_map[unit] <- scales::alpha("grey", 0.5)
      line_type_map[unit] <- "solid"
    }
  }
  
  # Custom legend breaks and labels
  legend_breaks <- c("Unit_1", "Synthetic_control", "Unit_2")  # Correct the value to match the data
  legend_labels <- c("observed target unit", "synthetic control target unit", "donor pool units")
  
  #data_all$Unit <- factor(data_all$Unit, levels = c("Grey Line 1", "Grey Line 2", "Black Line")) # Adjust to your actual labels
  
  
  # Plot the results - reference domain
  ggplot(data_all, aes(x = year, y = Outcome, group = Unit)) +
    geom_line(data = subset(data_all, Unit != "Unit_1" & Unit != "Synthetic_control"), 
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.5) +  # Plot grey donor pool lines first
    geom_line(data = subset(data_all, Unit == "Synthetic_control"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Dashed black line on top
    geom_line(data = subset(data_all, Unit == "Unit_1"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    scale_color_manual(values = color_map, 
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map, 
                          breaks = legend_breaks,
                          labels = legend_labels) +
    scale_y_continuous(limits = c(5, 50)) +
    labs(title = NULL,
         x = "Time",
         y = "Outcome in the reference domain",
         color = "Legend",
         linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.15), # Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
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
    guides(color = guide_legend(title = NULL, override.aes = list(
      linetype = c("solid", "dashed", "solid"),
      color = c("black", "black", "grey")),
      order = 1),
      linetype = guide_legend(title = NULL, override.aes = list(
        color = c("black", "black", "grey")),
        order = 1))
  
  ggsave(filename = paste0(dir$output, "/F_plot_multilines.pdf"), width = 8, height = 6, dpi = 300)
}

synthetic_target_plot <- function(Y, w, s_max, i_max){
    
  synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]
  
  # Create a data frame for plotting all units in F
  data_all <- data.frame(
    year = rep(1:s_max, i_max),
    Outcome = as.vector(t(Y)),
    Unit = rep(paste0("Unit_", 1:i_max), each = s_max)
  )
  
  # Add the synthetic treated unit to the data frame
  data_synthetic <- data.frame(
    year = 1:s_max,
    Outcome = as.vector(synthetic_treated_Y),
    Unit = rep("Synthetic_control", s_max)
  )
  
  data_observed <- data.frame(
    year = 1:s_max,
    Outcome = as.vector(Y[1, ] + causal),
    Unit = rep("Observed_target", s_max)
  )
  
  # Combine the data frames
  data_all <- rbind(data_all, data_synthetic, data_observed)
  
  # Define colors and line types
  color_map <- c("Unit_1" = "red", "Synthetic_control" = "black","Observed_target" = "black")
  line_type_map <- c("Unit_1" = "solid", "Synthetic_control" = "dashed", "Observed_target" = "solid")
  
  for (unit in unique(data_all$Unit)) {
    if (!unit %in% names(color_map)) {
      color_map[unit] <- scales::alpha("grey", 0.5)
      line_type_map[unit] <- "solid"
    }
  }
  
  # Custom legend breaks and labels
  legend_breaks <- c("Unit_1", "Synthetic_control", "Observed_target","Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("counterfactual target unit", "synthetic control target unit", "observed target unit", "donor pool units")
  
  # Plot the results - reference domain
  ggplot(data_all, aes(x = year, y = Outcome, group = Unit, color = Unit, linetype = Unit)) +
    geom_line(data = subset(data_all, !Unit %in% c("Unit_1", "Synthetic_control", "Observed_target")), 
              aes(color = Unit, linetype = Unit), size = 0.8, alpha = 0.5) +  # Grey donor pool lines first
    geom_line(data = subset(data_all, Unit == "Synthetic_control"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Dashed black line on top
    geom_line(data = subset(data_all, Unit == "Observed_target"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Solid black line on top
    geom_line(data = subset(data_all, Unit == "Unit_1"), 
              aes(color = Unit, linetype = Unit), size = 0.8) +  # Red solid line on top
    scale_color_manual(values = color_map, 
                       breaks = legend_breaks,
                       labels = legend_labels) +
    scale_linetype_manual(values = line_type_map, 
                          breaks = legend_breaks,
                          labels = legend_labels) +
    scale_y_continuous(limits = c(5, 50)) +
    labs(title = NULL,
         x = "Time",
         y = "Outcome in the target domain",
         color = "Legend",
         linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.2), # Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
          legend.text = element_text(size = 12), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(
      linetype = c("solid", "dashed", "solid", "solid"),
      color = c("red", "black","black", "grey")),
      order = 1),
      linetype = guide_legend(title = NULL, override.aes = list(
        color = c("red", "black","black", "grey")),
        order = 1))
  
  ggsave(filename = paste0(dir$output, "/Y_plot_multilines.pdf"), width = 8, height = 6, dpi = 300)
}


# --------------------------------------------------------------------#
# ------------------   Assumption Plot -------------------------------#

linear_equi_conf_plot <- function(F, Y){
  
  # normalized plot: 
  nor = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ])
  
  Y_0 = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ])/nor
  Y_1 = mean(Y[1, ] + causal - varepsilon[1, ])/nor # observed target unit
  Y_0_t = mean(Y[1, ])/nor
  
  F_1 = mean(F[1, ] - epsilon[1, ])/nor
  F_0 = mean(F[2:(J+1), ] - epsilon[2:(J+1), ])/nor
  
  
  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(Y_0, F_0, Y_1, F_1, Y_0_t),
    group = factor(c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'), levels = c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'))
  )
  
  # Map labels for plotting
  labels <- c('Y_0' = "E * '[' * Y[0] * ']'", 'F_0' = "E * '[' * F[0] * ']'", 'Y_1' = "E * '[' * Y[1] * ']'", 'F_1' = "E * '[' * F[1] * ']'", 'Y_0_t' = "E * '[' * Y[1]^(0) * ']'")
  
  # Plot with modified expression
  ggplot(data, aes(x = time, y = value)) +
    geom_point(size = 1.5) +
    geom_text(aes(label = labels[group]), hjust = 0.5, vjust = 0, parse = TRUE, size = 5) +  # Adding labels to points
    geom_line(aes(linetype = group), size = 0.8) +
    geom_line(data = data.frame(time = c(0, 1), value = c(F_1, Y_0_t), group = 'Y_0_t'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    geom_line(data = data.frame(time = c(1, 0), value = c(Y_0, F_0), group = 'Y_0'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    scale_x_continuous(breaks = c(0, 1), labels = c('', '')) +
    #scale_y_continuous(limits = c(1, 1.4)) +
    scale_linetype_manual(
      values = c('Y_0' = "solid", 'F_0' = "solid", 'Y_1' = "solid", 'F_1' = "solid", 'Y_0_t' = "solid"),
      labels = c('Y_0' = expression(Y[0]), 'F_0' = expression(E * "[" * F[0] * "]"), 'Y_1' = expression(Y[1]), 'F_1' = expression(F[1]), 'Y_0_t' = expression(Y[1]^(0)))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(text = element_text(family = "sans"),
          legend.position = "none", 
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    theme(legend.title = element_blank())
  
  ggsave(filename = paste0(dir$output, "/linear_eq.pdf"), width = 8, height = 6, dpi = 300)
}

logorithmic_equi_conf_plot <- function(F, Y){
  
  nor_a = log(mean(Y[2:(J+1),] -varepsilon[2:(J+1),]))
  
  log_Y_1 = log(mean(Y[1,] + causal-varepsilon[1,]))/nor_a # observed target unit
  log_Y_0_t = log(mean(Y[1,]))/nor_a
  
  log_Y_0 = log(mean(Y[2:(J+1),] - varepsilon[2:(J+1), ]))/nor_a
  
  log_F_1 = log(mean(F[1,] - epsilon[1,]))/nor_a
  log_F_0 = log(mean(F[2:(J+1),] - epsilon[2:(J+1),]))/nor_a

  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(log_Y_0, log_F_0, log_Y_1, log_F_1, log_Y_0_t),
    group = factor(c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'), levels = c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'))
  )
  
  # Map labels for plotting
  labels <- c('Y_0' = "log(E * '[' * Y[0] * ']')", 'F_0' = "log(E * '[' * F[0] * ']')", 'Y_1' = "log(E * '[' * Y[1] * ']')", 'F_1' = "log(E * '[' * F[1] * ']')", 'Y_0_t' = "log(E * '[' * Y[1]^(0) * ']')")
  
  # Plot with modified expression
  ggplot(data, aes(x = time, y = value)) +
    geom_point(size = 1.5) +
    geom_text(aes(
      label = labels[group], 
      hjust = ifelse(time == 0, -0.1, 0),  # Left side labels move right, right side labels move left
      vjust = ifelse(time == 0, 0, -0.2)  
    ), parse = TRUE, size = 5) +  # Adding labels to points
    geom_line(aes(linetype = group), size = 0.8) +
    geom_line(data = data.frame(time = c(0, 1), value = c(log_F_1, log_Y_0_t), group = 'Y_0_t'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    geom_line(data = data.frame(time = c(1, 0), value = c(log_Y_0, log_F_0), group = 'Y_0'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    scale_x_continuous(limits = c(0, 1.25)) +
    scale_y_continuous(limits = c(1, 1.4)) +
    scale_linetype_manual(
      values = c('log_Y_0' = "solid", 'log_F_0' = "solid", 'log_Y_1' = "solid", 'log_F_1' = "solid", 'log_Y_0_t' = "solid"),
      labels = c('log_Y_0' = expression(Y[0]), 'log_F_0' = expression(E * "[" * F[0] * "]"), 
                 'log_Y_1' = expression(Y[1]), 'log_F_1' = expression(F[1]), 'log_Y_0_t' = expression(Y[1]^(0)))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = "none", 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm"),
          plot.margin = margin(30, 60, 10, 10)) +  # Added margin to prevent label cutoff
    theme(legend.title = element_blank())
  
  ggsave(filename = paste0(dir$output, "/log_eq.pdf"), width = 8, height = 6, dpi = 300)
}


# ----------------------------------------------------------------- #
# --------------------------    NSE plot   ------------------------ #
# for simulation results 
NSE_plot_ZXF <- function(results_list){
  # --------------------------NSE plot ------------------------ #
  results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
    res <- results_list[[i]]
    data.frame(
      t_max = t_max_values[i],
      mean_mu_nse = res$mean_mu_nse,
      mean_Z_nse = res$mean_Z_nse,
      mean_X_nse = res$mean_X_nse,
      mean_F_nse = res$mean_F_nse
    )
  }))
  
  # Reshape the data to long format
  results_long <- melt(results_df, id.vars = "t_max", 
                       measure.vars = c("mean_mu_nse", "mean_Z_nse", "mean_X_nse", "mean_F_nse"), 
                       variable.name = "Metric", value.name = "Value")
  
  # Create the box plot without outliers
  ggplot(results_long, aes(x = as.factor(t_max), y = Value, fill = Metric)) +
    geom_boxplot(outlier.shape = NA) +  # Remove outliers
    scale_fill_manual(values = c("mean_mu_nse" = "green", "mean_Z_nse" = "purple", "mean_X_nse" = "orange", "mean_F_nse" = "pink"), 
                      labels = c("mean_mu_nse" = expression(NSE(mu[1])),
                                 "mean_Z_nse" = expression(NSE(Z[1])),
                                 "mean_X_nse" = expression(NSE(X[1])),
                                 "mean_F_nse" = expression(NSE(F[1])))) +
    #scale_y_continuous(limits = c(0,6)) +
    labs(title = NULL,
         x = "Time",
         y = "NSE") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.95), # Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
          legend.text = element_text(size = 16), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +  # Adjust the aspect ratio
    guides(fill = guide_legend(title = NULL), 
           color = guide_legend(title = NULL))
  
  ggsave(filename = paste0(dir$output, "/nse_mu_Z_X_F.pdf"), width = 8, height = 6, dpi = 300)
  
}

NSE_plot_ZX <- function(results_list){
  
  # --------------------------NSE plot without F ------------------------ #
  results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
    res <- results_list[[i]]
    data.frame(
      t_max = t_max_values[i],
      mean_mu_nse = res$mean_mu_nse,
      mean_Z_nse = res$mean_Z_nse,
      mean_X_nse = res$mean_X_nse
    )
  }))
  
  # Reshape the data to long format
  results_long <- melt(results_df, id.vars = "t_max", 
                       measure.vars = c("mean_mu_nse", "mean_Z_nse", "mean_X_nse"), 
                       variable.name = "Metric", value.name = "Value")
  
  # Create the box plot without outliers
  ggplot(results_long, aes(x = as.factor(t_max), y = Value, fill = Metric)) +
    geom_boxplot(outlier.shape = NA) +  # Remove outliers
    scale_fill_manual(values = c("mean_mu_nse" = "green", "mean_Z_nse" = "purple", "mean_X_nse" = "orange"), 
                      labels = c("mean_mu_nse" = expression(NSE(mu[1])),
                                 "mean_Z_nse" = expression(NSE(Z[1])),
                                 "mean_X_nse" = expression(NSE(X[1])))) +
    scale_y_continuous(limits = c(0, 0.08)) +
    labs(title = NULL,
         x = "Time",
         y = "NSE") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.95),# Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
          legend.text = element_text(size = 16), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm"))+
    guides(fill = guide_legend(title = NULL), 
           color = guide_legend(title = NULL))
  
  ggsave(filename = paste0(dir$output, "/nse_mu_Z_X.pdf"), width = 8, height = 6, dpi = 300)
  
}



# ------------------------------------------------------------------- #
# ----------------------- Placebo Test Plot  ------------------------ #

Placebo_test <- function(F, Y, J, t_max, s_max, w, eta_Z, eta_X){ #X, Z, dr, dt
  
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
  
  # Causal effect
  causal <- causal
  
  synthetic_placebo_df_Y <- rbind(synthetic_placebo_df_Y, data.frame(
    year = 1:s_max,
    Outcome = as.numeric(Y[1,] + causal - synthetic_treated_Y),
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
  legend_labels <- c("target unit", "control units")
  
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
    scale_y_continuous(limits = c(-10, 10)) +
    labs(title = NULL, x = "Time", y = "Causal Effect", color = "Legend", linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.1), 
          legend.justification = c(1, 1), 
          legend.text = element_text(size = 12), 
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
  
  ggsave(filename = paste0(dir$output, "/placebo_F.pdf"), width = 8, height = 6, dpi = 300)
  
  legend_breaks <- c("Synthetic_Treated_Y", "Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("target unit", "control units")
  
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
          legend.text = element_text(size = 12), 
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

  ggsave(filename = paste0(dir$output, "/placebo_Y.pdf"), width = 8, height = 6, dpi = 300)
  print(plot_F)
  print(plot_Y)
}


# --------------------------------------------------------------------#
# ------------------   Assumption Plot -------------------------------#

linear_equi_conf_plot <- function(F, Y){
  
  # normalized plot: 
  nor = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ])
  
  Y_0 = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ])/nor
  Y_1 = mean(Y[1, ] + causal - varepsilon[1, ])/nor # observed target unit
  Y_0_t = mean(Y[1, ])/nor
  Y_0 = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ] )/nor
  
  F_1 = mean(F[1, ] - epsilon[1, ])/nor
  F_0 = mean(F[2:(J+1), ] - epsilon[2:(J+1), ])/nor
  
  
  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(Y_0, F_0, Y_1, F_1, Y_0_t),
    group = factor(c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'), levels = c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'))
  )
  
  # Map labels for plotting
  labels <- c('Y_0' = "E * '[' * Y[0] * ']'", 'F_0' = "E * '[' * F[0] * ']'", 'Y_1' = "E * '[' * Y[1] * ']'", 'F_1' = "E * '[' * F[1] * ']'", 'Y_0_t' = "E * '[' * Y[1]^(0) * ']'")
  
  # Plot with modified expression
  ggplot(data, aes(x = time, y = value)) +
    geom_point(size = 1.5) +
    geom_text(aes(label = labels[group]), hjust = 0.2, vjust = -0.5, parse = TRUE, size = 7) +  # Adjusted vjust for better separation
    geom_line(aes(linetype = group), size = 0.8) +
    geom_line(data = data.frame(time = c(0, 1), value = c(F_1, Y_0_t), group = 'Y_0_t'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    geom_line(data = data.frame(time = c(1, 0), value = c(Y_0, F_0), group = 'Y_0'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    scale_x_continuous(limits = c(0, 1.05)) +
    scale_y_continuous(limits = c(1, 1.4)) +
    scale_linetype_manual(
      values = c('Y_0' = "solid", 'F_0' = "solid", 'Y_1' = "solid", 'F_1' = "solid", 'Y_0_t' = "solid"),
      labels = c('Y_0' = expression(Y[0]), 'F_0' = expression(E * "[" * F[0] * "]"), 'Y_1' = expression(Y[1]), 'F_1' = expression(F[1]), 'Y_0_t' = expression(Y[1]^(0)))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = "none", 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    theme(legend.title = element_blank())
  
  ggsave(filename = paste0(dir$output, "/linear_eq.pdf"), width = 8, height = 6, dpi = 300)
}

logarithmic_equi_conf_plot <- function(F, Y) {
  
  nor_a <- log(mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ]))
  
  log_Y_1 <- log(mean(Y[1, ] + causal - varepsilon[1, ])) / nor_a # observed target unit
  log_Y_0_t <- log(mean(Y[2:(J+1), ])) / nor_a
  
  log_Y_0 <- log(mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ])) / nor_a
  
  log_F_1 <- log(mean(F[1, ] - epsilon[1, ])) / nor_a
  log_F_0 <- log(mean(F[2:(J+1), ] - epsilon[2:(J+1), ])) / nor_a
  
  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(log_Y_0, log_F_0, log_Y_1, log_F_1, log_Y_0_t),
    group = factor(c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'), 
                   levels = c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'))
  )
  
  # Correct labels to match group names in data
  labels <- c(
    'log_Y_0' = "log(E * '[' * Y[0] * ']')", 
    'log_F_0' = "log(E * '[' * F[0] * ']')", 
    'log_Y_1' = "log(E * '[' * Y[1] * ']')", 
    'log_F_1' = "log(E * '[' * F[1] * ']')", 
    'log_Y_0_t' = "log(E * '[' * Y[1]^(0) * ']')"
  )
  
  ggplot(data, aes(x = time, y = value)) +
    geom_point(size = 1.5) +
    geom_text(aes(label = labels[group]), hjust = 0.5, vjust = -1, parse = TRUE, size = 7) +  # Adjusted vjust for better separation
    geom_line(aes(linetype = group), size = 0.8) +
    geom_line(data = data.frame(time = c(0, 1), value = c(log_F_1, log_Y_0_t), group = 'Y_0_t'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    geom_line(data = data.frame(time = c(1, 0), value = c(log_Y_0, log_F_0), group = 'Y_0'), 
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    #scale_x_continuous(breaks = c(0, 1), labels = c('', '')) +
    scale_x_continuous(limits = c(0, 1.25)) +
    scale_y_continuous(limits = c(1, 1.4)) +
    scale_linetype_manual(
      values = c('log_Y_0' = "solid", 'log_F_0' = "solid", 'log_Y_1' = "solid", 'log_F_1' = "solid", 'log_Y_0_t' = "solid"),
      labels = c('log_Y_0' = expression(Y[0]), 'log_F_0' = expression(E * "[" * F[0] * "]"), 'log_Y_1' = expression(Y[1]), 'log_F_1' = expression(F[1]), 'log_Y_0_t' = expression(Y[1]^(0)))
    ) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = "none", 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    theme(legend.title = element_blank())
}




# ----------------------------------------------------------------- #
# --------------------------    Bias plot   ------------------------ #
bias_plot <- function(results_list){
  # Combine results into a data frame
  results_df <- do.call(rbind, lapply(1:length(results_list), function(i) {
    res <- results_list[[i]]
    data.frame(
      t_max = t_max_values[i],
      mean_bias_linear = res$mean_bias_linear,
      mean_bias_log = res$mean_bias_log,
      mean_bias_sc = res$mean_bias_sc
    )
  }))
  
  # Reshape the data to long format
  results_long <- melt(results_df, id.vars = "t_max", 
                       measure.vars = c( "mean_bias_linear", "mean_bias_log", "mean_bias_sc"), 
                       variable.name = "Metric", value.name = "Value")
  
  # Create the box plot without outliers
  ggplot(results_long, aes(x = as.factor(t_max), y = Value, fill = Metric)) +
    geom_boxplot(outlier.shape = NA) +  # Remove outliers
    scale_fill_manual(values = c("mean_bias_log" = "red", "mean_bias_sc" = "grey", "mean_bias_linear" = "blue"), 
                      labels = c("mean_bias_linear" = expression(hat(psi)^eq1),
                                 "mean_bias_log" = expression(hat(psi)^eq2),
                                 "mean_bias_sc" = expression(hat(psi)^sc))) +
    scale_color_manual(values = c("mean_bias_log" = "red", "mean_bias_sc" = "grey", "mean_bias_linear" = "blue"), 
                       labels = c("mean_bias_linear" = expression(hat(psi)^eq1),
                                  "mean_bias_log" = expression(hat(psi)^eq2),
                                  "mean_bias_sc" = expression(hat(psi)^sc))) +
    labs(title = NULL,
         x = "Time",
         y = "Bias") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.95), # Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
          legend.text = element_text(size = 16), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha("grey", 0.5)),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(fill = guide_legend(title = NULL), 
           color = guide_legend(title = NULL))
  
  ggsave(filename = paste0(dir$output, "/bias_sc_eq.pdf"), width = 8, height = 6, dpi = 300)
  
}


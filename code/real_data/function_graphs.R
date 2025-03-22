# ------------------------------------------------------------------------------#
# Function: Generate Reference and Target Domain Plot for Synthetic Control
# ------------------------------------------------------------------------------#

synthetic_reference_plot_real_data <- function(F, w, t_max, i_max){
  
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
  legend_labels <- c("Chelsea", "Synthetic Chelsea", "19 Cities")
  
  # Plot the results - reference domain
  p <- ggplot(data_all, aes(x = year, y = Outcome, group = Unit, color = Unit, linetype = Unit)) +
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
    labs(title = NULL,
         x = "Time",
         y = "Covid-19 Vaccination Rate of Black Subpopulation",
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
          legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
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
  
  print(p)
  ggsave(filename = paste0(dir$output, "/covid_ref.pdf"), width = 8, height = 6, dpi = 300)
}

synthetic_target_plot_real_data <- function(Y, w, s_max, i_max){
  
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
    Outcome = as.vector(Y_treated),
    Unit = rep("Observed_target", s_max)
  )
  
  # Combine the data frames
  data_all <- rbind(data_all, data_synthetic, data_observed)
  
  # Define colors and line types
  color_map <- c("Unit_1" = "black", "Synthetic_control" = "black","Observed_target" = "black")
  line_type_map <- c("Unit_1" = "solid", "Synthetic_control" = "dashed", "Observed_target" = "solid")
  
  for (unit in unique(data_all$Unit)) {
    if (!unit %in% names(color_map)) {
      color_map[unit] <- scales::alpha("grey", 0.5)
      line_type_map[unit] <- "solid"
    }
  }
  
  # Custom legend breaks and labels
  legend_breaks <- c("Unit_1", "Synthetic_control","Unit_2")  # Add "Unit_2" to represent control units
  legend_labels <- c("Chelsea", "Synthetic Chelsea", "19 Cities")
  
  # Plot the results - reference domain
  p <- ggplot(data_all, aes(x = year, y = Outcome, group = Unit, color = Unit, linetype = Unit)) +
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
    labs(title = NULL,
         x = "Time",
         y = "Covid 19 Vaccination Rate of Hispanic",
         color = "Legend",
         linetype = "Legend") +
    theme_minimal() +
    theme(text = element_text(family = "sans"), 
          legend.position = c(0.95, 0.2), # Adjust this to position the legend inside the graph
          legend.justification = c(1, 1), # Adjust this to align the legend box
          legend.text = element_text(size = 14), 
          legend.key.size = unit(1, "lines"),
          legend.key.width = unit(2, "lines"), 
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill = scales::alpha('grey', 0.5)),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.ticks.length = unit(0.25, "cm")) +
    guides(color = guide_legend(title = NULL, override.aes = list(
      linetype = c("solid", "dashed", "solid"),
      color = c("black", "black","grey")),
      order = 1),
      linetype = guide_legend(title = NULL, override.aes = list(
        color = c( "black","black", "grey")),
        order = 1))
  
  print(p)
  ggsave(filename = paste0(dir$output, "/covid_target.pdf"), width = 8, height = 6, dpi = 300)
  
}

# ------------------------------------------------------------------------------#
# Function: Plot Linear Equi-Confounding Data
# ------------------------------------------------------------------------------#

linear_equi_conf_plot_real_data <- function(Y_treated, Y_control, F_treated, F_control) {
  # Compute means for observed and control units
  H_1 <- mean(Y_treated)  # Observed target unit
  H_0 <- mean(Y_control)
  B_1 <- mean(F_treated)
  B_0 <- mean(F_control)
  
  # Prepare data for visualization
  data <- data.frame(
    time = c(1, 0, 1, 0),
    value = c(H_0, B_0, H_1, B_1),
    group = factor(c('H_0', 'B_0', 'H_1', 'B_1'), levels = c('H_0', 'B_0', 'H_1', 'B_1'))
  )
  
  # Define labels for plotting
  labels <- c('H_0' = "H[0]", 'B_0' = "B[0]", 'H_1' = "H[1]", 'B_1' = "B[1]")
  
  # Generate the plot
  p <- ggplot(data, aes(x = time, y = value)) +
    geom_point(size = 1) +
    geom_text(aes(label = labels[group]), vjust = -0.2, parse = TRUE) +  # Adding labels to points
    geom_line(aes(linetype = group), size = 0.8) +
    
    # Add connecting lines
    geom_line(data = data.frame(time = c(1, 0), value = c(H_0, B_0), group = 'B_0'),
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    geom_line(data = data.frame(time = c(0, 1), value = c(B_1, H_1), group = 'B_0'),
              aes(x = time, y = value), linetype = "solid", color = "black", size = 0.8) +
    
    # Customize x-axis labels
    scale_x_continuous(breaks = c(0, 1), labels = c('Reference', 'Target')) +
    
    # Define line types and labels
    scale_linetype_manual(
      values = c('H_0' = "solid", 'B_0' = "solid", 'H_1' = "solid", 'B_1' = "solid"),
      labels = c('H_0' = expression(H[0]), 'B_0' = expression(B[0]),
                 'H_1' = expression(H[1]), 'B_1' = expression(B[1]))
    ) +
    
    # Labels and Theme
    labs(title = NULL, x = NULL, y = 'Outcome') +
    theme_minimal() +
    theme(
      text = element_text(family = "Times New Roman"),
      legend.position = "none",
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.length = unit(0.25, "cm"),
      legend.title = element_blank()
    )
  print(p)
}


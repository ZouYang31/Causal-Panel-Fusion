# ------------------------------------------------------------------------------#
# Function: Calculate Synthetic Estimated Treatment Effect (ETT)
# ------------------------------------------------------------------------------#
synthetic_ett <- function(Y_treated, Y, w, J) {
  # Calculate the synthetic treated unit
  synthetic_treated_Y <- t(w) %*% Y[2:(J + 1), ]
  
  # Compute the Estimated Treatment Effect (ETT)
  ett_synthetic <- Y_treated - synthetic_treated_Y
  
  # Return the ETT result
  return(mean(ett_synthetic))
}

# ------------------------------------------------------------------------------#
# Function: Compute Estimated Treatment Effect (ETT) using Equi-Confounding Data Fusion
# ------------------------------------------------------------------------------#
linear_equi_confounding_ett <- function(Y_treated, Y_control, F_treated, F_control) {
  
  # Compute means for observed and control units
  Y_1 <- mean(Y_treated)  # Observed target unit
  Y_0 <- mean(Y_control)
  F_1 <- mean(F_treated)
  F_0 <- mean(F_control)
  
  # Compute the Estimated Treatment Effect (ETT) for Equi-Confounding
  ett <- (Y_1 - Y_0) - (F_1 - F_0)
  
  return(ett)
}

log_equi_confounding_ett <- function(Y_treated, Y_control, F_treated, F_control) {
  
  # Normalize using log of the mean of control group
  nor_a <- log(mean(Y_control))
  
  # Compute log-transformed values
  log_Y_1 <- log(mean(Y_treated)) / nor_a  # Observed target unit
  log_Y_0 <- log(mean(Y_control)) / nor_a
  log_B_1 <- log(mean(F_treated)) / nor_a
  log_B_0 <- log(mean(F_control)) / nor_a
  
  # Compute the Estimated Treatment Effect (ETT) for Logarithm Equi-Confounding
  estimate <- mean(Y_treated) - (mean(F_treated) / mean(F_control)) * mean(Y_control)
  
  return(estimate)
}





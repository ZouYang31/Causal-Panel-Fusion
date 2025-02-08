# ------------------------------------------------------------------------------#
# Function: Perform Sensitivity Analysis for Equi-Confounding Data Fusion
# ------------------------------------------------------------------------------#
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/0_function_simulation.R")


# 
sensitivity_analysis <- function(
    B, F_treated, F_control, X_treated, X_control, Z_treated, Z_control,
    F, X, Z,
    t_max, dr, dt, num_city,
    NSE_Z_baseline, NSE_X_baseline,
    E_eta_Z = c(0.02, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0),
    E_eta_X = c(0.02, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0)
) {
  
  # Initialize output matrices
  ll1 <- length(E_eta_Z)
  ll0 <- length(E_eta_X)
  EST_NSEs <- matrix(0, nrow = ll1, ncol = ll0)
  b.lists <- vector("list", ll1 * ll0)
  
  # Iterate over eta_Z and eta_X values
  for (i in 1:ll1) {
    for (j in 1:ll0) {
      
      # Step 1: Find the best B for given eta_Z and eta_X
      b_list <- find_best_B(
        B, F_treated = F_treated, F_control = F_control,
        X_treated = as.matrix(X_treated), X_control = as.matrix(X_control),
        Z_treated = as.matrix(Z_treated), Z_control = as.matrix(Z_control),
        t_max = t_max, dr = dr, dt = dt, i_max = num_city,
        NSE_Z_baseline = NSE_Z_baseline, NSE_X_baseline = NSE_X_baseline,
        eta_Z = E_eta_Z[i], eta_X = E_eta_X[j]
      )
      
      # Store b_list
      index <- (i - 1) * ll0 + j
      b.lists[[index]] <- b_list$best_B_list
      
      # Step 2: Solve for w using the best B values
      solution.w <- b_list$best_w_star
      # Step 3: Compute the loss function
      loss.w <- as.numeric(
        b_list$best_B_list$b_F * (sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max) +
          b_list$best_B_list$b_Z * (sum((matrix(Z_treated) - t(Z_control) %*% solution.w)^2) / dr) +
          b_list$best_B_list$b_X * (sum((matrix(X_treated) - t(X_control) %*% solution.w)^2) / dt)
      )
      
      # Step 4: Store the computed loss function result
      EST_NSEs[i, j] <- round(loss.w, 2)
    }
  }
  
  # Convert results into data frames for output
  res_df_est <- data.frame(EST_NSEs)
  colnames(res_df_est) <- paste0("X_", E_eta_X)
  row.names(res_df_est) <- paste0("Z_", E_eta_Z)
  
  # Convert b.lists into a human-readable data frame
  b_lists_df <- matrix("", nrow = ll1, ncol = ll0)
  
  for (i in 1:ll1) {
    for (j in 1:ll0) {
      index <- (i - 1) * ll0 + j
      b_list <- b.lists[[index]]
      
      if (!is.null(b_list)) {
        b_F <- round(b_list$b_F, 2)
        b_Z <- round(b_list$b_Z, 2)
        b_X <- round(b_list$b_X, 2)
        b_lists_df[i, j] <- paste0("b_F=", b_F, ", b_Z=", b_Z, ", b_X=", b_X)
      } else {
        b_lists_df[i, j] <- "NA"
      }
    }
  }
  
  res_df_b_lists <- data.frame(b_lists_df)
  colnames(res_df_b_lists) <- paste0("X_", E_eta_X)
  row.names(res_df_b_lists) <- paste0("Z_", E_eta_Z)
  
  # Return results as a list
  return(list(
    EST_NSEs = res_df_est,
    B_lists = res_df_b_lists
  ))
}

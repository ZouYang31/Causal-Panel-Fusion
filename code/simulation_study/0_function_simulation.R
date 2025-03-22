##########################################
# Data Fusion Simulation Study function
# Author: Zoey Yang
# Date: 
#########################################

# Call the function sources
source(paste0(dir$code, "/2_function_optim_w.R"))


# ------------------------- define the function -------------------------
### ---------------------------------------------------------------------------------

simulation_B_list <- function(time, B, eta_Z, eta_X, alg = FALSE) {

  print(paste0("Function testing:", time))
  
  # Initialized parameter:
  # n, rhi, phi, theta
  
  # Initialize the simulation vector to store MSE values
  simulation_bias_sc <- numeric(n)
  simulation_bias_linear <- numeric(n)
  simulation_bias_log <- numeric(n)
  mu_nse_values <- numeric(n)
  Z_nse_values <- numeric(n)
  X_nse_values <- numeric(n)
  F_nse_values <- numeric(n)
  
  # Potential outcomes
  
  for (S in seq(1, n)) {
    # ------------ reference domain -----------------
    # sd = 2 
    epsilon <- matrix(rnorm(i_max * end_point, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = end_point) 
    
    # ------------ target domain --------------------
    # sd = 0.5
    varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)
    
    # Calculate F_it^N for each combination of i and t
    for (t in 1:end_point){
      for (i in 1:i_max) {
        F[i, t] <- rho[t] + sum(phi[t, ] * Z[i, ]) + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
      }
    }
    
    for (s in 1:s_max) {
      for (i in 1:i_max) {
        Y[i, s] <- varrho[s] + sum(varphi[s, ] * X[i, ]) + sum(vartheta[s, ] * mu[i, ]) + varepsilon[i, s]
      }
    }
    
    # Calculate the weights
    J <- i_max - 1
    w <- rep(0, J)
    
    # Without loss of generality, choose i = 1 as the treated unit and j = 1 as the control unit
    F_treated <- F[1, 1:time]
    F_control <- F[2:(J+1), 1:time]
    Z_treated <- Z[1, ]
    Z_control <- Z[2:(J+1), ]
    X_treated <- X[1, ]
    X_control <- X[2:(J+1), ]
    mu_treated <- mu[1, ]
    mu_control <- mu[2:(J+1), ]
    
    
    # # Apply algorithm if `alg` is TRUE
    if (alg) {
      
      # Calculate the NSE_baselines
      result_X <- optimize_w_ipop(F_treated, F_control, X, Z, time, dr, dt, i_max, target = "X")
      wX <- result_X$weights
      NSE_X_baseline <- NSE_x(wX, F[, 1:time], X, Z, time, dr, dt, i_max, target = "X")
      
      
      result_Z <- optimize_w_ipop(F_treated, F_control, X, Z, time, dr, dt, i_max, target = "Z")
      wZ <- result_Z$weights
      NSE_Z_baseline <- NSE_x(wZ, F[, 1:time], X, Z, time, dr, dt, i_max, target = "Z")
      
      
      # Call the function
      B_result <- find_best_B(B = B, F_treated = F_treated, F_control = F_control, 
                              X_treated = X_treated, X_control = X_control,
                              Z_treated = Z_treated, Z_control = Z_control,
                              t_max = time, dr = dr, dt = dt, i_max = i_max, 
                              NSE_Z_baseline = NSE_Z_baseline, 
                              NSE_X_baseline = NSE_X_baseline, 
                              eta_Z = eta_Z, eta_X = eta_X) 
      optimal_par <-B_result$best_B_list
      
      w <- B_result$best_w_star
      
      # Check if `optimal_par` contains NA values
      if (any(is.na(optimal_par))) {
        cat("Optimization failed. Using initial values for b_F, b_Z, b_X.\n")
        b_F <- initial_B[1]
        b_Z <- initial_B[2]
        b_X <- initial_B[3]
      } else {
        b_F <- optimal_par$b_F
        b_Z <- optimal_par$b_Z
        b_X <- optimal_par$b_X
      }
    
      print(S)
    
    }  
  
    # Extract the weights
    # w <- solve_w(F = F, X = X, Z = Z, 
    #              b_F = b_F, b_Z = b_Z, b_X = b_X,
    #              i_max = i_max, t_max = time, dr = dr, dt = dt)$X
    #print("Checking w is: ") 
    #print(w_check)
    
    
    # Calculate the NSEs
    mu_nse_values[S] <- sum((mu_treated - t(w) %*% mu_control)^2)/du
    Z_nse_values[S] <- sum((Z_treated - t(w) %*% Z_control)^2)/dr
    X_nse_values[S] <- sum((X_treated - t(w) %*% X_control)^2)/dt
    F_nse_values[S] <- sum((F_treated - t(w) %*% F_control)^2)/time
    
    
    # Create the synthetic treated unit for the target domain
    synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]
    
    Y_null <- mean(Y[1, ])
    Y_1 <- mean(Y[1, ] + causal)
    Y_0 <- mean(Y[2:(J+1), ])
    F_1 <- mean(F_treated)
    F_0 <- mean(F_control)
    
    groud_truth <- mean(causal)
    estimate_linear <- (Y_1 - Y_0)-(F_1 - F_0)
    estimate_log <- Y_1 - (F_1/F_0) * Y_0
    
    # Calculate average NSE for this simulation
    simulation_bias_sc[S] <- mean(Y[1, ] - synthetic_treated_Y)
    simulation_bias_linear[S] <- groud_truth - estimate_linear
    simulation_bias_log[S] <- groud_truth - estimate_log
    
  }
  
  # Calculate the mean MSE across all simulations
  mean_bias_sc <- mean(simulation_bias_sc)
  mean_bias_linear <- mean(simulation_bias_linear)
  mean_bias_log <- mean(simulation_bias_log)
  mean_mu_nse <- mean(mu_nse_values)
  mean_Z_nse <- mean(Z_nse_values)
  mean_X_nse <- mean(X_nse_values)
  mean_F_nse <- mean(F_nse_values)
  
  
  
  ci_bias_sc <- simulation_bias_sc
  ci_bias_linear <- simulation_bias_linear
  ci_bias_log <- simulation_bias_log
  ci_mu_nse <- mu_nse_values
  ci_Z_nse <- Z_nse_values
  ci_X_nse <- X_nse_values
  ci_F_nse <- F_nse_values
  
  # Note: change back;  
  return(list(mean_bias_log = simulation_bias_log, ci_bias_log = ci_bias_log, 
              mean_mu_nse = mu_nse_values, ci_mu_nse = ci_mu_nse,
              mean_Z_nse = Z_nse_values, ci_Z_nse = ci_Z_nse,
              mean_bias_sc = simulation_bias_sc, ci_bias_sc = ci_bias_sc,
              mean_bias_linear = simulation_bias_linear, ci_bias_linear = ci_bias_linear,
              mean_X_nse = X_nse_values, ci_X_nse = ci_X_nse,
              mean_F_nse = F_nse_values, ci_F_nse = ci_F_nse
#              F = F, Y = Y
              ))
}

### ---------------------------------------------------------------------------------


simulation_w_matrix <- function(time) {
  
  # Initialize the simulation vector to store MSE values
  simulation_bias_sc <- numeric(n)
  simulation_bias_linear <- numeric(n)
  simulation_bias_log <- numeric(n)
  mu_nse_values <- numeric(n)
  Z_nse_values <- numeric(n)
  X_nse_values <- numeric(n)
  F_nse_values <- numeric(n)
  
  # Potential outcomes
  
  for (S in seq(1, n)) {

    # ------------ reference domain -----------------
    # sd = 2 
    epsilon <- matrix(rnorm(i_max * end_point, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = end_point)
    
    # ------------ target domain --------------------
    # sd = 0.5
    varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)
    
    # Calculate F_it^N for each combination of i and t
    for (t in 1:end_point){
      for (i in 1:i_max) {
        F[i, t] <- rho[t] + sum(phi[t, ] * Z[i, ]) + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
      }
    }
    
    for (t in 1:s_max) {
      for (i in 1:i_max) {
        Y[i, t] <- varrho[t] + sum(varphi[t, ] * X[i, ]) + sum(vartheta[t, ] * mu[i, ]) + varepsilon[i, t]
      }
    }
    
    # Calculate the weights
    J <- i_max - 1
    w <- rep(0, J)
    
    # Without loss of generality, choose i = 1 as the treated unit and j = 1 as the control unit
    F_treated <- F[1, 1:time]
    F_control <- F[2:(J+1), 1:time]
    Z_treated <- Z[1, ]
    Z_control <- Z[2:(J+1), ]
    X_treated <- X[1, ]
    X_control <- X[2:(J+1), ]
    mu_treated <- mu[1, ]
    mu_control <- mu[2:(J+1), ]
    
    # Extract the weights
    w_matrix <- solve_w_matrix(F = F, X = X, Z = Z)
    w <- ifelse(w_matrix < 10^(-5), 0, w_matrix)
    print(w)
    
    # Calculate the NSEs
    mu_nse_values[S] <- sum((mu_treated - t(w) %*% mu_control)^2)/du
    Z_nse_values[S] <- sum((Z_treated - t(w) %*% Z_control)^2)/dr
    X_nse_values[S] <- sum((X_treated - t(w) %*% X_control)^2)/dt
    F_nse_values[S] <- sum((F_treated - t(w) %*% F_control)^2)/time
    
    # Create the synthetic treated unit for the target domain
    synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]
    
    Y_null <- mean(Y[1, ])
    Y_1 <- mean(Y[1, ] + causal)
    Y_0 <- mean(Y[2:(J+1), ])
    F_1 <- mean(F_treated)
    F_0 <- mean(F_control)
    
    groud_truth <- mean(causal)
    estimate_linear <- (Y_1 - Y_0)-(F_1 - F_0)
    estimate_log <- Y_1 - (F_1/F_0) * Y_0
    
    # Calculate average NSE for this simulation
    simulation_bias_sc[S] <- mean(Y[1, ] - synthetic_treated_Y)
    simulation_bias_linear[S] <- groud_truth - estimate_linear
    simulation_bias_log[S] <- groud_truth - estimate_log
    
  }
  
  # Calculate the mean MSE across all simulations
  mean_bias_sc <- mean(simulation_bias_sc)
  mean_bias_linear <- mean(simulation_bias_linear)
  mean_bias_log <- mean(simulation_bias_log)
  mean_mu_nse <- mean(mu_nse_values)
  mean_Z_nse <- mean(Z_nse_values)
  mean_X_nse <- mean(X_nse_values)
  mean_F_nse <- mean(F_nse_values)
  
  
  
  ci_bias_sc <- simulation_bias_sc
  ci_bias_linear <- simulation_bias_linear
  ci_bias_log <- simulation_bias_log
  ci_mu_nse <- mu_nse_values
  ci_Z_nse <- Z_nse_values
  ci_X_nse <- X_nse_values
  ci_F_nse <- F_nse_values
  
  # Note: change back;  
  return(list(mean_bias_log = simulation_bias_log, ci_bias_log = ci_bias_log, 
              mean_mu_nse = mu_nse_values, ci_mu_nse = ci_mu_nse,
              mean_Z_nse = Z_nse_values, ci_Z_nse = ci_Z_nse,
              mean_bias_sc = simulation_bias_sc, ci_bias_sc = ci_bias_sc,
              mean_bias_linear = simulation_bias_linear, ci_bias_linear = ci_bias_linear,
              mean_X_nse = X_nse_values, ci_X_nse = ci_X_nse,
              mean_F_nse = F_nse_values, ci_F_nse = ci_F_nse))
}

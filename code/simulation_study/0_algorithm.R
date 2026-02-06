##########################################
# Data Fusion Simulation Study function
# Author: Zoey Yang
# Date: 
#########################################


# Call the function sources
source(paste0(dir$code, "/2_function_optim_w.R"))
library(optimx)
#library(Synth)
library(gurobi)

# Setup the gurobi function
#Sys.setenv(GRB_LICENSE_FILE = "your_path_to_gurobi_license/gurobi.lic")
#install.packages("/Library/gurobi1200/macos_universal2/R/gurobi_12.0-0_R_4.4.1.tgz", repos = NULL, type = "source")


# ------------ define the data driven algorithm -----------------
### --------------- Algorithm B-list ------------------- ###
find_best_B <- function(B, F_treated, F_control, X_treated, X_control, Z_treated, Z_control,
                        t_max, dr, dt, i_max, 
                        NSE_Z_baseline, NSE_X_baseline, eta_Z, eta_X) {
  
  # Initialize variables to store the least loss and corresponding B
  least_loss_B <- Inf
  best_B_list <- NULL
  best_w_star <- NULL
  
  # Iterate over all combinations in result_list
  for (B_list in B) {
    # Calculate the loss for the current B_list
    result <- optimize_w_b_gurobi_B_list(
      par = B_list,  # Replace with appropriate parameters if needed
      F_treated = F_treated,
      F_control = F_control,
      X_treated = X_treated,
      X_control = X_control,
      Z_treated = Z_treated,
      Z_control = Z_control,
      t_max = t_max,
      dr = dr,
      dt = dt,
      i_max = i_max,
      NSE_Z_baseline = NSE_Z_baseline,
      NSE_X_baseline = NSE_X_baseline,
      eta_Z = eta_Z,
      eta_X = eta_X
    )
    loss_B <- result$loss.B
    
    w_star <- result$solution.w
    
    # if (S == 40){
    # print(loss_B)}
    
    # Update least loss and best B if the current loss is smaller
    if (!is.na(loss_B) && loss_B < least_loss_B) {
      least_loss_B <- loss_B
      best_B_list <- B_list
      best_w_star <- w_star
    }
  }
  
  # if (S == 40){
  # print(best_B_list)}
  
  # Return the best B_list and the corresponding least loss as a list
  return(list(best_B_list = best_B_list, least_loss_B = least_loss_B, best_w_star = best_w_star))
}

simulation_B_list_3DGP <- function(time, B, eta_Z, eta_X, alg = FALSE, timeout_seconds = 300) {
  
  print(paste0("Function testing:", time))
  
  library(parallel)
  
  # Initialize the simulation vector to store MSE values
  simulation_bias_sc <- numeric(n)
  simulation_bias_linear <- numeric(n)
  simulation_bias_log <- numeric(n)
  mu_nse_values <- numeric(n)
  Z_nse_values <- numeric(n)
  X_nse_values <- numeric(n)
  F_nse_values <- numeric(n)
  
  C <- mean(rho[1:time]) / mean(varrho)
  varrho_local <- C * varrho
  
  theta_local <- theta[1:time, ]
  delta_vec = colMeans(theta_local) - colMeans(vartheta)
  vartheta_local <- sweep(vartheta, 2, delta_vec, "+")
  
  for (S in seq(1, n)) {
    # ------------ reference domain -----------------
    epsilon <- matrix(rnorm(i_max * end_point, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = end_point) 
    
    # ------------ target domain --------------------
    varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)
    
    for (t in 1:end_point){
      for (i in 1:i_max) {
        F[i, t] <- rho[t] + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
      }
    }
    
    for (s in 1:s_max) {
      for (i in 1:i_max) {
        Y[i, s] <- varrho_local[s] + sum(vartheta_local[s, ] * mu[i, ]) + varepsilon[i, s]
      }
    }
    
    J <- i_max - 1
    w <- rep(0, J)
    
    F_treated <- F[1, 1:time]
    F_control <- F[2:(J+1), 1:time]
    Z_treated <- Z[1, ]
    Z_control <- Z[2:(J+1), ]
    X_treated <- X[1, ]
    X_control <- X[2:(J+1), ]
    mu_treated <- mu[1, ]
    mu_control <- mu[2:(J+1), ]
    
    if (alg) {
      
      result_X <- optimize_w_ipop(F_treated, F_control, X, Z, time, dr, dt, i_max, target = "X")
      wX <- result_X$weights
      NSE_X_baseline <- NSE_x(wX, F[, 1:time], X, Z, time, dr, dt, i_max, target = "X")
      
      result_Z <- optimize_w_ipop(F_treated, F_control, X, Z, time, dr, dt, i_max, target = "Z")
      wZ <- result_Z$weights
      NSE_Z_baseline <- NSE_x(wZ, F[, 1:time], X, Z, time, dr, dt, i_max, target = "Z")
      
      # ============ TIMEOUT USING mcparallel (5 min = 300 sec) ============
      job <- mcparallel({
        find_best_B(B = B, F_treated = F_treated, F_control = F_control, 
                    X_treated = X_treated, X_control = X_control,
                    Z_treated = Z_treated, Z_control = Z_control,
                    t_max = time, dr = dr, dt = dt, i_max = i_max, 
                    NSE_Z_baseline = NSE_Z_baseline, 
                    NSE_X_baseline = NSE_X_baseline, 
                    eta_Z = eta_Z, eta_X = eta_X)
      })
      
      # Wait for result with timeout
      B_result <- mccollect(job, wait = TRUE, timeout = timeout_seconds)
      
      # Check if timed out (result will be NULL if timeout)
      if (is.null(B_result) || length(B_result) == 0) {
        cat(sprintf("Timeout at iteration S=%d. Using zero weights.\n", S))
        # Kill the process if still running
        try(tools::pskill(job$pid), silent = TRUE)
        w <- rep(1/J, J)
        b_F <- initial_B[1]
        b_Z <- initial_B[2]
        b_X <- initial_B[3]
      } else {
        # Extract result from list
        B_result <- B_result[[1]]
        
        # Check if result is an error or NULL
        if (is.null(B_result) || inherits(B_result, "try-error")) {
          cat(sprintf("Error at iteration S=%d. Using zero weights.\n", S))
          w <- rep(1/J, J)
          b_F <- initial_B[1]
          b_Z <- initial_B[2]
          b_X <- initial_B[3]
        } else {
          optimal_par <- B_result$best_B_list
          w <- B_result$best_w_star
          
          if (any(is.na(optimal_par))) {
            cat("Optimization failed. Using initial values.\n")
            b_F <- initial_B[1]
            b_Z <- initial_B[2]
            b_X <- initial_B[3]
          } else {
            b_F <- optimal_par$b_F
            b_Z <- optimal_par$b_Z
            b_X <- optimal_par$b_X
          }
        }
      }
      # ============ END TIMEOUT ============
      
      print(S)
    }  
    
    mu_nse_values[S] <- sum((mu_treated - t(w) %*% mu_control)^2)/du
    Z_nse_values[S] <- sum((Z_treated - t(w) %*% Z_control)^2)/dr
    X_nse_values[S] <- sum((X_treated - t(w) %*% X_control)^2)/dt
    F_nse_values[S] <- sum((F_treated - t(w) %*% F_control)^2)/time
    
    synthetic_treated_Y <- t(w) %*% Y[2:(J+1), ]
    
    Y_null <- mean(Y[1, ])
    Y_1 <- mean(Y[1, ] + causal)
    Y_0 <- mean(Y[2:(J+1), ])
    F_1 <- mean(F_treated)
    F_0 <- mean(F_control)
    
    groud_truth <- mean(causal)
    estimate_linear <- (Y_1 - Y_0)-(F_1 - F_0)
    estimate_log <- Y_1 - (F_1/F_0) * Y_0
    
    simulation_bias_sc[S] <- mean(Y[1, ] - synthetic_treated_Y)
    simulation_bias_linear[S] <- groud_truth - estimate_linear
    simulation_bias_log[S] <- groud_truth - estimate_log
  }
  
  ci_bias_sc <- simulation_bias_sc
  ci_bias_linear <- simulation_bias_linear
  ci_bias_log <- simulation_bias_log
  ci_mu_nse <- mu_nse_values
  ci_Z_nse <- Z_nse_values
  ci_X_nse <- X_nse_values
  ci_F_nse <- F_nse_values
  
  return(list(mean_bias_log = simulation_bias_log, ci_bias_log = ci_bias_log, 
              mean_mu_nse = mu_nse_values, ci_mu_nse = ci_mu_nse,
              mean_Z_nse = Z_nse_values, ci_Z_nse = ci_Z_nse,
              mean_bias_sc = simulation_bias_sc, ci_bias_sc = ci_bias_sc,
              mean_bias_linear = simulation_bias_linear, ci_bias_linear = ci_bias_linear,
              mean_X_nse = X_nse_values, ci_X_nse = ci_X_nse,
              mean_F_nse = F_nse_values, ci_F_nse = ci_F_nse
  ))
}


optimize_w_b_gurobi_B_list <- function(par = stop("B missing"),  
                                            F_treated = stop("F_treated missing"), F_control = stop("F_control missing"),
                                       X_treated = stop("X missing"), X_control = stop("X missing"),
                                       Z_treated = stop("Z missing"), Z_control = stop("Z missing"), 
                                       t_max = stop("t_max missing"), 
                                       dr = stop("dr missing"), dt = stop("dt missing"), i_max = stop("i_max missing"), 
                                       NSE_Z_baseline = stop("NSE_Z_baseline missing"), 
                                       NSE_X_baseline = stop("NSE_X_baseline missing"), 
                                       eta_Z = 0.1, eta_X = 0.1) {
  
  J <- i_max - 1
  
  # Extract treated and control matrices
  F_treated <- F_treated[1:t_max]
  F_control <- F_control[, 1:t_max]
  
  # Normalize B so that sum(B) = 1
  b_F <- par$b_F
  b_Z <- par$b_Z
  b_X <- par$b_X
  
  # Quadratic matrix for the objective function
  H <- b_F/t_max * (F_control %*% t(F_control)) +
    b_Z/dr * (Z_control %*% t(Z_control)) +
    b_X/dt * (X_control %*% t(X_control))
  
  # Linear term for the objective function
  c <- -2 * (b_F/t_max * t(F_treated) %*% t(F_control) +
               b_Z/dr * t(Z_treated) %*% t(Z_control) +
               b_X/dt * t(X_treated) %*% t(X_control))
  
  ## -------- constraints -------- ##
  scale_n = 1
  
  # Quadratic constraints for Z and X
  Q_Z <- Z_control %*% t(Z_control) * scale_n
  L_Z <- -2 * t(matrix(Z_treated)) %*% t(Z_control) * scale_n
  #rhs_Z <- dr *  eta_Z * NSE_Z_baseline - sum(Z_treated^2)
  rhs_Z <- dr * ((1 + eta_Z) * (1 + NSE_Z_baseline * scale_n) - 1) - sum(Z_treated^2) * scale_n
  
  Q_X <- X_control %*% t(X_control) * scale_n
  L_X <- -2 * t(matrix(X_treated)) %*% t(X_control) * scale_n
  #rhs_X <- dt * eta_X * NSE_X_baseline - sum(X_treated^2)
  rhs_X <- dt * ((1 + eta_X) * (1 + NSE_X_baseline * scale_n) - 1) - sum(X_treated^2) * scale_n
  
  ## -------- Gurobi Model -------- ##
  model <- list()
  
  # Objective function (quadratic)
  model$Q <- H
  model$obj <- as.vector(c)
  model$modelsense <- "min"
  
  # Quadratic constraints (Z and X)
  model$quadcon <- list(
    list(Qc = Q_Z, q = as.vector(L_Z), rhs = rhs_Z),
    list(Qc = Q_X, q = as.vector(L_X), rhs = rhs_X)
  )
  
  # Linear constraint: sum(w) = 1
  model$A <- matrix(1, nrow = 1, ncol = J)  # Row of ones for the sum constraint
  model$rhs <- c(1)
  model$sense <- c("=")
  
  
  # Variable bounds
  model$lb <- rep(0, J)  # Non-negativity bounds
  model$ub <- rep(1, J)  # No upper bound
  
  # Solve the model using Gurobi
  params <- list(OutputFlag = 0, TimeLimit = 100) # debug model false

  # Extract weights
  result <- tryCatch({
    gurobi(model, params = params)
  }, error = function(e) {
    return(NA)
  })
  
  # Check if result is NA due to timeout or other issues
  if (is.null(result) || !is.list(result) || is.null(result$status) || result$status != "OPTIMAL") {
    return(list(loss.B = NA, solution.w = NA))
  }
  
  # Extract weights
  solution.w <- as.numeric(result$x)
  if (is.null(solution.w) || any(is.na(solution.w))) {
    return(list(loss.B = NA, solution.w = NA))
  }
  
  solution.w[solution.w < 0.005] <- 0
  
  # Calculate the loss for w
  if (any(is.na(solution.w))) {
    loss.w <- NA
    loss.B <- NA
  } else {
    loss.w <- as.numeric(
      b_F * (sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max) +
        b_Z * (sum((matrix(Z_treated) - t(Z_control) %*% solution.w)^2) / dr) +
        b_X * (sum((matrix(X_treated) - t(X_control) %*% solution.w)^2) / dt) 
    )
    
    # Return the loss for B
    loss.B <- as.numeric(sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max)
  }

  return(list(loss.B = loss.B, solution.w = solution.w))
}



### --------------- Backup algorithm ------------------- ###
# Data-driven weight optimization algorithm
## ------------------------- QCQP with gurobi -------------------------------##
# Objective function 




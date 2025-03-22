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


### --------------- Backup algorithm ------------------- ###
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

optimize_w_b_gurobi <- function(par = stop("B missing"), F = stop("F missing"), X = stop("X missing"), 
                                Z = stop("Z missing"), t_max = stop("t_max missing"), 
                                dr = stop("dr missing"), dt = stop("dt missing"), i_max = stop("i_max missing"), 
                                NSE_Z_baseline = stop("NSE_Z_baseline missing"), 
                                NSE_X_baseline = stop("NSE_X_baseline missing"), 
                                eta_Z = 10, eta_X = 10) {
  
  J <- i_max - 1
  
  # Extract treated and control matrices
  F_treated <- F[1, 1:t_max]
  F_control <- F[2:(J+1), 1:t_max]
  Z_treated <- Z[1, ]
  Z_control <- Z[2:(J+1), ]
  X_treated <- X[1, ]
  X_control <- X[2:(J+1), ]
  
  # Normalize B so that sum(B) = 1
  par <- abs(par) / sum(abs(par))
  B <- par
  b_F <- B[1]
  b_Z <- B[2]
  b_X <- B[3]
  
  # Quadratic matrix for the objective function
  H <- b_F/t_max * (F_control %*% t(F_control)) +
    b_Z/dr * (Z_control %*% t(Z_control)) +
    b_X/dt * (X_control %*% t(X_control))
  
  # Linear term for the objective function
  c <- -2 * (b_F/t_max * t(F_treated) %*% t(F_control) +
               b_Z/dr * t(Z_treated) %*% t(Z_control) +
               b_X/dt * t(X_treated) %*% t(X_control))
  
  ## -------- constraints -------- ##
  # Quadratic constraints for Z and X
  Q_Z <- Z_control %*% t(Z_control)
  L_Z <- -2 * t(matrix(Z_treated)) %*% t(Z_control)
  #rhs_Z <- dr *  eta_Z * NSE_Z_baseline - sum(Z_treated^2)
  rhs_Z <- dr * ((1 + eta_Z) * (1 + NSE_Z_baseline) - 1) - sum(Z_treated^2)
  
  Q_X <- X_control %*% t(X_control)
  L_X <- -2 * t(matrix(X_treated)) %*% t(X_control)
  #rhs_X <- dt * eta_X * NSE_X_baseline - sum(X_treated^2)
  rhs_X <- dt * ((1 + eta_X) * (1 + NSE_X_baseline) - 1) - sum(X_treated^2)
  
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
  params <- list(OutputFlag = 0) # debug model false
  #params$InfUnbdInfo <- 1
  result <- gurobi(model, params = params)
  
  # Extract weights
  solution.w <- as.numeric(result$x)
  
  solution.w[solution.w < 0.005] <- 0  # Threshold for small weights
  #print(solution.w)
  
  # Calculate the loss for w
  loss.w <- as.numeric(
    b_F * (sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max) +
      b_Z * (sum((matrix(Z_treated) - t(Z_control) %*% solution.w)^2) / dr) +
      b_X * (sum((matrix(X_treated) - t(X_control) %*% solution.w)^2) / dt)
  )
  
  print("The optimal value for NSEs")
  print(loss.w)
  print("The optimal value for NSE_Z")
  print(sum((matrix(Z_treated) - t(Z_control) %*% solution.w)^2))
  print("The optimal value for NSE_X")
  print(sum((matrix(X_treated) - t(X_control) %*% solution.w)^2))
  # Return the loss for B
  loss.B <- as.numeric(sum((matrix(F_treated) - t(F_control) %*% solution.w)^2) / t_max)
  
  return(invisible(loss.B))
}




# Solving the parameters using optimization
optimize_parameters <- function(SV1, optimize_fn, F, X, Z, t_max, dr, dt, i_max, 
                                NSE_Z_baseline, NSE_X_baseline, eta_Z = 10, eta_X = 10) {
  # Run the optimization using optimx with non-negativity constraints
  rgV.optim <- optimx(
    par = SV1,                     # Initial guess for parameters
    fn = optimize_fn,              # Objective function c("BFGS", "Nelder-Mead"),  # Choose optimization methods
    method = "L-BFGS-B",           # Constrained optimization method
    lower = rep(0, length(SV1)),   # Lower bounds (non-negative)
    upper = rep(1, length(SV1)), # Upper bounds (unbounded)
    control = list(
      kkt = FALSE,                 # Disable KKT checks
      starttests = FALSE,          # Skip start tests
      dowarn = FALSE,               # Suppress warnings
      trace = 2,                    # Enable debugging output
      factr = 1e4,                 # stricter convergence 1e5
      pgtol = 1e-3                 # gradients
    ),
    # Pass additional arguments required by optimize_fn
    F = F,                         # Matrix for treated and control outcomes (F)
    X = X,                         # Matrix for covariates (X)
    Z = Z,                         # Matrix for covariates (Z)
    t_max = t_max,                 # Maximum time (t_max)
    dr = dr,                       # Dimensionality for Z
    dt = dt,                       # Dimensionality for X
    i_max = i_max,                 # Maximum index for control
    NSE_Z_baseline = NSE_Z_baseline, # Baseline NSE for Z
    NSE_X_baseline = NSE_X_baseline, # Baseline NSE for X
    eta_Z = eta_Z,                 # Constraint parameter for Z
    eta_X = eta_X                  # Constraint parameter for X
  )
  
  # Find the row index with the smallest value
  best_row <- which.min(rgV.optim$value)
  
  # Extract the optimal parameters as a numeric vector
  optimal_par <- as.numeric(rgV.optim[best_row, c("p1", "p2", "p3")])
  
  # Clamp the parameters to be non-negative (additional safeguard)
  optimal_par <- pmax(optimal_par, 0)
  
  # Normalize the parameters to ensure they sum to 1
  optimal_par <- optimal_par / sum(optimal_par)
  
  # Debugging output for validation
  #cat("Optimal parameters after optimization:", optimal_par, "\n")
  
  # Return the optimal parameters
  return(optimal_par)
}


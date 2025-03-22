##########################################
# Data Fusion Simulation Study function
# Author: Zoey Yang
# Date: 
#########################################

# Optimization for W:

solve_w <- function(F, X, Z, b_F, b_X, b_Z, i_max, t_max, dr, dt){
  
  ## ---- test for w and see the visualization -----  ##
  J <- i_max - 1
  w <- rep(0, J)
  
  # Without loss of generation, we choose i = 1 as the treated unit and j = 1 as the control unit. 
  F_treated <- F[1, 1:t_max]
  F_control <- F[2:(J+1), 1:t_max]
  Z_treated <- Z[1, ]
  Z_control <- Z[2:(J+1), ]
  X_treated <- X[1, ]
  X_control <- X[2:(J+1), ]

  weight_Z = sqrt((t_max + dr + dt)/dr * b_Z) 
  weight_X = sqrt((t_max + dr + dt)/dt * b_X) 
  weight_F = sqrt((t_max + dr + dt)/t_max * b_F) 
  
  treated_all <- c(weight_F*F_treated, weight_Z*Z_treated, weight_X*X_treated)
  
  # set up the matrix for the system of equations 
  A_F <- t(F_control)  
  A_Z <- t(Z_control)  
  A_X <- t(X_control)
  
  print(b_F)
  print(b_Z)
  print(b_X)
  
  A_all <- rbind(weight_F*A_F, weight_Z*A_Z, weight_X*A_X)
  
  # Add a row of ones to ensure the sum of weights is 1 
  E <- rep(1, J)
  f <- c(1)
  
  # Inequality constraint matrix for positivity
  G <- diag(J)
  h <- rep(0, J)  # All elements must be >= 0
  
  ## testing ------------------------
  if (any(is.na(A_all)) || any(is.nan(A_all)) || any(is.infinite(A_all))) {
    stop("A_all contains invalid values.")
  }
  if (any(is.na(treated_all)) || any(is.nan(treated_all)) || any(is.infinite(treated_all))) {
    stop("treated_all contains invalid values.")
  }
  if (any(is.na(E)) || any(is.nan(E)) || any(is.infinite(E))) {
    stop("E contains invalid values.")
  }
  if (any(is.na(f)) || any(is.nan(f)) || any(is.infinite(f))) {
    stop("f contains invalid values.")
  }
  if (any(is.na(G)) || any(is.nan(G)) || any(is.infinite(G))) {
    stop("G contains invalid values.")
  }
  if (any(is.na(h)) || any(is.nan(h)) || any(is.infinite(h))) {
    stop("h contains invalid values.")
  }
  
  # Solve the least squares problem with constraints
  result <- lsei(A = A_all, B = treated_all, E = E, F = f, G = G, H = h)
  
  return(result)
}


## --------------------------------------------------------  ##
##                         Weighting Matrix
## --------------------------------------------------------  ##
# Abadies 2010 method
# Using synth method for solving w.
solve_w_matrix <- function(F, X, Z){
    
    F_treated <- F[1, ]
    F_control <- F[2:(J+1), ]
    Z_treated <- Z[1, ]
    Z_control <- Z[2:(J+1), ]
    X_treated <- X[1, ]
    X_control <- X[2:(J+1), ]
    
    #predictor 
    X1 = matrix(c(X_treated, Z_treated))
    X0 = rbind(t(X_control), t(Z_control))
    
    # outcome from pre-intervention 
    Z1 = matrix(F[1, ])
    Z0 = t(F[2:(J+1), ])
    
    w <- as.vector(synth(X1 = X1, X0 = X0, Z1 = Z1, Z0 = Z0)$solution.w)
    
    return(w)
  
}


## --------------------- Calculate the NSEs ------------------------------- ##
NSE_x <- function(w, F, X, Z, t_max, dr, dt, i_max, target = "X") {
  J <- i_max - 1
  if (target == "F") {
    return(sum((F[1, 1:t_max] - t(w) %*% F[2:(J+1), 1:t_max])^2) / t_max)
  } else if (target == "X") {
    return(sum((X[1, ] - t(w) %*% X[2:(J+1), ])^2) / dt)
  } else if (target == "Z") {
    return(sum((Z[1, ] - t(w) %*% Z[2:(J+1), ])^2) / dr)
  } else {
    stop("Invalid target. Choose 'F', 'X', or 'Z'.")
  }
}


## --------------------- Solving for the  ------------------------------- ##
optimize_w_ipop <- function(F_treated, F_control, X, Z, t_max, dr, dt, i_max, target = "X", 
                            margin.ipop = 1e-4, sigf.ipop = 5, bound.ipop = 10) {
  
  J <- i_max - 1
  
  # Extract treated and control data
  F_treated <- F_treated[1:t_max]
  F_control <- F_control[, 1:t_max]
  Z_treated <- Z[1, ]
  Z_control <- Z[2:(J+1), ]
  X_treated <- X[1, ]
  X_control <- X[2:(J+1), ]
  
  # Define objective components based on target
  if (target == "F") {
    X_treated <- matrix(F_treated)
    X_control <- t(F_control)
    V <- diag(rep(1 / t_max, t_max))  # Scaling factor for F
  } else if (target == "X") {
    X_treated <- matrix(X_treated)
    X_control <- t(X_control)
    V <- diag(rep(1 / dt, dt))  # Scaling factor for X
  } else if (target == "Z") {
    X_treated <- matrix(Z_treated)
    X_control <- t(Z_control)
    V <- diag(rep(1 / dr, dr))  # Scaling factor for Z
  } else {
    stop("Invalid target. Choose 'F', 'X', or 'Z'.")
  }
  
  # Define quadratic programming components
  H <- t(X_control) %*% V %*% X_control  # Quadratic term
  
  c <- -1 * t(X_treated) %*% V %*% X_control  # Linear term
  A <- t(rep(1, J))  # Equality constraint (sum w = 1)
  b <- 1
  l <- rep(0, J)  # Lower bound (w >= 0)
  u <- rep(1, J)  # Upper bound (optional, large value)
  r <- 0
  
  # Solve the quadratic programming problem
  res <- ipop(c = c, H = H, A = A, b = b, l = l, u = u, r = r, 
              bound = bound.ipop, margin = margin.ipop, sigf = sigf.ipop)
  
  # Extract weights and compute losses
  solution.w <- as.numeric(primal(res))
  # Set the threshold for w
  solution.w[solution.w < 0.005] <- 0
  loss.w <- sum((X_treated - X_control %*% solution.w)^2)
  
  # Return the optimized weights and loss
  return(list(weights = solution.w, loss = loss.w))
}


##-------------------------------------------##
# Data Fusion Simulation Study 
# Author: Zoey Yang
# Date: 
##-------------------------------------------##

# -----------------------------------------------------------------
# Clean up work space and load or install necessary packages if necessary
rm(list=ls())
want <- c("dplyr", "ggplot2", "limSolve", "scales", "reshape2", "pbapply", "nloptr", "Synth", "gurobi", "kernlab", "optimx")
need <- want[!(want %in% installed.packages()[,"Package"])]
if(length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only = TRUE))
rm(want, need)

# Working directories
dir <- list()
dir$root <- getwd()
dir$output <- paste(dir$root, "/output", sep= "")
dir$code <- paste(dir$root, "/code", sep= "")


##---------------------   Initialization   ----------------------##
##---------------------   Synthetic Control Data Generating Process   ----------------------##
set.seed(45) #45


# Initialization parameters
i_max <- 31  # Number of units (donor pool + treated unit)
t_max <- 20  # Reference domain (time points)
s_max <- 5   # Target domain (time points)

# Dimensions for covariates and latent factors
dr <- 3  # Number of observed covariates in the reference domain
dt <- 3  # Number of observed covariates in the target domain
du <- 3  # Number of unobserved latent factors

# Generate latent factors and observed covariates
mu <- matrix(runif(i_max * du, min=0, max=1), nrow = i_max, ncol = du) # Latent confounder
Z <- matrix(runif(i_max * dr, min=0, max=1), nrow = i_max, ncol = dr)   # Reference domain covariates
X <- matrix(runif(i_max * dt, min=0, max=1), nrow = i_max, ncol = dt)   # Target domain covariates

# Generate domain-specific parameters
dist_max <- 10  # Max value for uniform distribution

causal <- sort(runif(s_max,min=2,max=5))

# Parameters for reference domain (F)
rho <- sort(runif(t_max, 0, 20))
#phi <- matrix(runif(dr * t_max, min=0, max=dist_max), nrow = t_max, ncol = dr)
theta <- matrix(runif(du * t_max, min=0, max=dist_max), nrow = t_max, ncol = du)
epsilon <- matrix(rnorm(i_max * t_max, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = t_max)

# Parameters for target domain (Y)
varrho <- sort(runif(s_max, 0, 10))
#varphi <- matrix(runif(dt * s_max, min=0, max=10), nrow = s_max, ncol = dt)
vartheta <- matrix(runif(du * s_max, min=0, max=10), nrow = s_max, ncol = du)
varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)

# Calculate the delta 
delta_vec = colMeans(theta) - colMeans(vartheta)
vartheta <- sweep(vartheta, 2, delta_vec, "+")

print("vartheta == ", mean(vartheta))
print(colMeans(vartheta))
print("theta == ", mean(theta))
print(colMeans(theta))

# Initialize potential outcomes
F <- matrix(0, nrow = i_max, ncol = t_max)
Y <- matrix(0, nrow = i_max, ncol = s_max)

# Calculate F_it^N for each combination of i and t
for (t in 1:t_max) {
  for (i in 1:i_max) {
    F[i, t] <- rho[t]  + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
  }
}

for (t in 1:s_max) {
  for (i in 1:i_max) {
    Y[i, t] <- varrho[t]  + sum(vartheta[t, ] * mu[i, ]) + varepsilon[i, t]
  }
}

# Calibration step to reduce the difference between core components for the treated unit and donor pool
# Compute the mean differences



mean(F[1, ]) - mean(F[2:i_max, ])
mean(Y[1, ]) - mean(Y[2:i_max, ])

linear_equi_conf_plot_DGP <- function(F, Y){
  
  # normalized plot: 
  
  Y_0 = mean(Y[2:i_max, ])
  Y_1 = mean(Y[1, ] + causal) # observed target unit
  Y_0_t = mean(Y[1, ])
  print(Y_0_t - Y_0)
  
  F_1 = mean(F[1, ])
  F_0 = mean(F[2:i_max,])
  print(F_1 - F_0)
  
  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(Y_0, F_0, Y_1, F_1, Y_0_t),
    group = factor(c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'), levels = c('Y_0', 'F_0', 'Y_1', 'F_1', 'Y_0_t'))
  )
  
  # Map labels for plotting
  labels <- c('Y_0' = "E * '[' * Y[0] * ']'", 'F_0' = "E * '[' * F[0] * ']'", 'Y_1' = "E * '[' * Y[1] * ']'", 'F_1' = "E * '[' * F[1] * ']'", 'Y_0_t' = "E * '[' * Y[1]^(0) * ']'")
  
  # Plot with modified expression
  p <- ggplot(data, aes(x = time, y = value)) +
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
  
  print(p)
  #ggsave(filename = paste0(dir$output, "/linear_eq.pdf"), width = 8, height = 6, dpi = 300)
}

##---------------   Synthetic Control Data Fusion   -------------##
# number of control unit 
J <- i_max - 1
w <- rep(0, J)

source(paste0(dir$code, "/0_function_simulation.R"))
source(paste0(dir$code, "/0_algorithm.R"))
source(paste0(dir$code, "/1_function_plot.R"))
source(paste0(dir$code, "/2_function_optim_w.R"))

#------------------------------------------------------------------------------#
###.                             Setup for alg                               ###
#------------------------------------------------------------------------------#
# Initialize an empty list to store the results
B <- list()

# Initialize counter
c <- 1

# Iterate over b_F, b_Z, and b_X ensuring they sum to 1
for (b_F in seq(0.01, 0.99, by = 0.02)) {
  for (b_Z in seq(0.01, 1 - b_F, by = 0.02)) {
    b_X <- 1 - b_F - b_Z  # Directly calculate b_X to satisfy the constraint
    if (b_X >= 0.02) {    # Ensure b_X meets the minimum step size
      # Add the current values of b_F, b_Z, and b_X as a list to result_list
      B[[c]] <- list(b_F = b_F, b_Z = b_Z, b_X = b_X)
      c <- c + 1
    } 
  }
}

# Calculate the baseline for NSE_Z and NSE_X
# number of control unit 

# Without loss of generation, we choose i = 1 as the treated unit and j = 1 as the control unit. 
F_treated <- F[1, ]
F_control <- F[2:(J+1), ]
Z_treated <- Z[1, ]
Z_control <- Z[2:(J+1), ]
X_treated <- X[1, ]
X_control <- X[2:(J+1), ]
mu_treated <- mu[1, ]
mu_control <- mu[2:(J+1), ]


result_X <- optimize_w_ipop(F_treated, F_control, X, Z, t_max, dr, dt, i_max, target = "X")
wX <- result_X$weights
NSE_X_baseline <- NSE_x(wX, F, X, Z, t_max, dr, dt, i_max, target = "X")
print(NSE_X_baseline)


result_Z <- optimize_w_ipop(F_treated, F_control, X, Z, t_max, dr, dt, i_max, target = "Z")
wZ <- result_Z$weights
NSE_Z_baseline <- NSE_x(wZ, F, X, Z, t_max, dr, dt, i_max, target = "Z")
print(NSE_Z_baseline)


# One test for the algorithm 
b_list <- find_best_B(B, F_treated, F_control, 
                      X_treated, X_control, Z_treated, Z_control, 
                      t_max, dr, dt, i_max, 
                      NSE_Z_baseline = NSE_Z_baseline, NSE_X_baseline = NSE_X_baseline, eta_Z=0.1, eta_X=0.1) 

b_list


##---------------   Synthetic Control Data Fusion   -------------##

# Solving the w
w <- b_list$best_w_star

print(w)
print(sum(w))
cat("w greater than 0:", sum(w>0), "\n")




##-------------------- Data Visualization ---------------------------##
# Calculate the synthetic control target unit for the reference domain
##-------------------------------------------------------------------##
synthetic_reference_plot(F = F, w = w, t_max = t_max, i_max = i_max)



##-------------------- Data Visualization --------------------------##
# Calculate the synthetic control target unit for the target domain
##------------------------------------------------------------------##
synthetic_target_plot(Y = Y, w = w, s_max = s_max, i_max = i_max)



#----------------------    Placebo Test   ---------------------------#
#--------------------------------------------------------------------#
#Placebo_test(F = F, Y = Y, J = J, t_max = t_max, s_max = s_max, w = w)

##---------------           Data Analysis         --------------------##
##---------------      Equi-confounding Data Fusion      -------------##
# plot for linear equi-confounding(Figure 1a)
linear_equi_conf_plot_DGP(F = F, Y = Y)



#------------------------------------------------------------#
##                       Simulation                        ###
#------------------------------------------------------------#

# number of simulation 
n = 300 #1000


initial <- 10
end_point <- 100
step_size <- 10

# Reinitialize outcome variables and parameters for simulation
F <- matrix(0, nrow = i_max, ncol = end_point)
Y <- matrix(0, nrow = i_max, ncol = s_max)

# Calculate the core components for each i and t or s
rho <- sort(runif(end_point, 0, 20))

phi <- matrix(runif(dr * end_point, min=0, max= dist_max), nrow = end_point, ncol = dr)
theta <- matrix(runif(du * end_point, min =0, max = dist_max), nrow = end_point, ncol = du)


# Define the range of t_max values
t_max_values <- seq(initial, end_point, by = step_size)



results_list <- pblapply(t_max_values, function(t){
  print(t)
  simulation_linear_equi(
    time = t)
})



save(results_list, file = paste0(dir$output, "results_list_DGP_linear0119.RData"))
#save(results_list, file = paste0(dir$output, "results_list0331.RData"))

load(paste0(dir$output, "results_list_DGP_linear0401_20.RData"))
results_list20 <- results_list

#load(paste0(dir$output, "results_list0331.RData"))
#results_list10 <- results_list

results_list10 <- results_list
combined_results <- c(results_list10, results_list20)
#------------------------------------------------------------------------#
##                          Graph Plot                         ###
#------------------------------------------------------------------------#
# --------------------------  Bias plot ------------------------------- #
bias_plot(results_list)


# --------------------------     NSE plot      ------------------------ #
NSE_plot_ZXF(combined_results)

NSE_plot_ZX(combined_results)



# simulation function 
simulation_B_list_DGP_equi <- function(time, B, eta_Z, eta_X, alg = FALSE) {
  
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
    
    for (t in 1:s_max) {
      for (i in 1:i_max) {
        Y[i, t] <- varrho[t] + sum(varphi[t, ] * X[i, ]) + sum(vartheta[t, ] * mu[i, ]) + varepsilon[i, t]
      }
    }
    
    # Compute the mean differences
    diff_F <- colMeans(F[1, , drop = FALSE] - colMeans(F[2:i_max, ]))
    diff_Y <- colMeans(Y[1, , drop = FALSE] - colMeans(Y[2:i_max, ]))
    
    # Compute the offset needed to match the average difference
    offset <- mean(diff_F) - mean(diff_Y)
    
    # Adjust core components to make them closer
    Y[1, ] <- Y[1, ] + offset
    
    mean(F[1, ]) - mean(F[2:i_max, ])
    mean(Y[1, ]) - mean(Y[2:i_max, ])
    
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
  ))
}

simulation_linear_equi <- function(time, alg = FALSE) {
  
  print(paste0("Function testing:", time))
  
  # Initialized parameter:
  # n, rhi, phi, theta
  
  # Initialize the simulation vector to store MSE values
  simulation_bias_linear <- numeric(n)
  simulation_bias_log <- numeric(n)
  
  theta_local <- theta[1:time, ]
  delta_vec = colMeans(theta_local) - colMeans(vartheta)
  vartheta_local <- sweep(vartheta, 2, delta_vec, "+")
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
        F[i, t] <- rho[t] + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
      }
    }
    
    for (s in 1:s_max) {
      for (i in 1:i_max) {
        Y[i, s] <- varrho[s] + sum(vartheta_local[s, ] * mu[i, ]) + varepsilon[i, s]
      }
    }
    
    # Calculate the weights
    J <- i_max - 1
    w <- rep(0, J)
    
    # Without loss of generality, choose i = 1 as the treated unit and j = 1 as the control unit
    F_treated <- F[1, 1:time]
    F_control <- F[2:(J+1), 1:time]
    
    
    Y_null <- mean(Y[1, ])
    Y_1 <- mean(Y[1, ] + causal)
    Y_0 <- mean(Y[2:(J+1), ])
    F_1 <- mean(F_treated)
    F_0 <- mean(F_control)
    
    groud_truth <- mean(causal)
    estimate_linear <- (Y_1 - Y_0)-(F_1 - F_0)
    estimate_log <- Y_1 - (F_1/F_0) * Y_0
    
    # Calculate average NSE for this simulation
    simulation_bias_linear[S] <- groud_truth - estimate_linear
    simulation_bias_log[S] <- groud_truth - estimate_log
    
  }
  
  # Calculate the mean MSE across all simulations
  mean_bias_linear <- mean(simulation_bias_linear)
  mean_bias_log <- mean(simulation_bias_log)
  
  
  ci_bias_linear <- simulation_bias_linear
  ci_bias_log <- simulation_bias_log
  
  # Note: change back;  
  return(list(mean_bias_log = simulation_bias_log, ci_bias_log = ci_bias_log, 
              mean_bias_linear = simulation_bias_linear, ci_bias_linear = ci_bias_linear
              #              F = F, Y = Y
  ))
}


library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------
# Extract mean and CI from results
# ------------------------------

df_summary <- data.frame(
  time = t_max_values,
  mean_linear = sapply(results_list, function(x) mean(x$mean_bias_linear)),
  mean_log = sapply(results_list, function(x) mean(x$mean_bias_log)),
  ci_lower_linear = sapply(results_list, function(x) quantile(x$ci_bias_linear, 0.025)),
  ci_upper_linear = sapply(results_list, function(x) quantile(x$ci_bias_linear, 0.975)),
  ci_lower_log = sapply(results_list, function(x) quantile(x$ci_bias_log, 0.025)),
  ci_upper_log = sapply(results_list, function(x) quantile(x$ci_bias_log, 0.975))
)

# ------------------------------
# Plot: Mean Bias with 95% CI
# ------------------------------

ggplot(df_summary, aes(x = time)) +
  # Linear estimator
  geom_ribbon(aes(ymin = ci_lower_linear, ymax = ci_upper_linear), 
              fill = "blue", alpha = 0.2) +
  geom_line(aes(y = mean_linear, color = "Linear"), linewidth = 1) +
  geom_point(aes(y = mean_linear, color = "Linear"), size = 2) +
  # Log estimator
  geom_ribbon(aes(ymin = ci_lower_log, ymax = ci_upper_log), 
              fill = "red", alpha = 0.2) +
  geom_line(aes(y = mean_log, color = "Log"), linewidth = 1) +
  geom_point(aes(y = mean_log, color = "Log"), size = 2) +
  # Reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Bias Comparison: Linear vs Log Estimator",
    x = "Time (T)",
    y = "Bias",
    color = "Estimator"
  ) +
  scale_color_manual(values = c("Linear" = "blue", "Log" = "red")) +
  theme_minimal() +
  theme(legend.position = "bottom")



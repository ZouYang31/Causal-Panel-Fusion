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
library(R.utils)
# Working directories
dir <- list()
dir$root <- getwd()
dir$output <- paste(dir$root, "/output", sep= "")
dir$code <- paste(dir$root, "/code", sep= "")


source(paste0(dir$code, "/0_function_simulation.R"))
source(paste0(dir$code, "/0_algorithm.R"))
source(paste0(dir$code, "/1_function_plot.R"))
Sys.setenv(GRB_LICENSE_FILE = "/Users/zoeyy/Desktop/bu/work_research/gurobi.lic")

##---------------------   Initialization   ----------------------##
##---------------------   Data Generating Process   ----------------------##

set.seed(45) 
# Try : 510, 420, 894, 649,
#45
#80 #99 #86
#95 good one. 
#45

# initialization 
i_max <- 31 #donor pool + 1 treated unit
t_max <- 20 # reference domain #20
s_max <- 5 # target domain

initial <- 10
end_point <- 100
step_size <- 10


dr <- 3 #3 #5  # Number of observed covariates reference domain
dt <- 3        # Number of observed covariates target domain
#Note: dr is not necessarilly equal to dt. 
du <- 3  #3  # Number of unobserved factor loadings d_u

mu <- matrix(runif(i_max * du, min=0, max=1), nrow = i_max, ncol = du) # latent confounder
Z <- matrix(runif(i_max * dr, min=0, max=0), nrow = i_max, ncol = dr) # reference domain
X <- matrix(runif(i_max * dt, min=0, max=0), nrow = i_max, ncol = dt) # target domain

# Generate parameters
dist_max <- 10 #max value for the uniform distribution

# ------------ reference domain -----------------#
rho <- sort(runif(t_max, 0, 20))

phi <- matrix(runif(dr * t_max, min=0, max= dist_max), nrow = t_max, ncol = dr)
theta <- matrix(runif(du * t_max, min =0, max = dist_max), nrow = t_max, ncol = du)
epsilon <- matrix(rnorm(i_max * t_max, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = t_max)

# #
# # ------------ target domain -------------------# 
varrho <- sort(runif(s_max, 0, 10))

#adjust the varrho. 
C <- mean(rho) / mean(varrho)

varrho <- C * varrho
mean(rho)
mean(varrho)

varphi <- matrix(runif(dt * s_max, min=0, max =10), nrow = s_max, ncol = dt)
vartheta <- matrix(runif(du * s_max, min =0, max = 10), nrow = s_max, ncol = du)
varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)

causal <- sort(runif(s_max,min=2,max=5))


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

log(mean(F[1, ])) - log(mean(F[2:i_max, ]))
log(mean(Y[1, ])) - log(mean(Y[2:i_max, ]))


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
logorithmic_equi_conf_plot_DGP <- function(F, Y){
  
  
  log_Y_1 = log(mean(Y[1,] + causal-varepsilon[1,])) # observed target unit
  log_Y_0_t = log(mean(Y[1,]))
  
  log_Y_0 = log(mean(Y[2:(J+1),] - varepsilon[2:(J+1), ]))
  
  log_F_1 = log(mean(F[1,] - epsilon[1,]))
  log_F_0 = log(mean(F[2:(J+1),] - epsilon[2:(J+1),]))
  
  data <- data.frame(
    time = c(1, 0, 1, 0, 1),
    value = c(log_Y_0, log_F_0, log_Y_1, log_F_1, log_Y_0_t),
    group = factor(c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'), levels = c('log_Y_0', 'log_F_0', 'log_Y_1', 'log_F_1', 'log_Y_0_t'))
  )
  
  # Map labels for plotting
  labels <- c('Y_0' = "log(E * '[' * Y[0] * ']')", 'F_0' = "log(E * '[' * F[0] * ']')", 'Y_1' = "log(E * '[' * Y[1] * ']')", 'F_1' = "log(E * '[' * F[1] * ']')", 'Y_0_t' = "log(E * '[' * Y[1]^(0) * ']')")
  
  # Plot with modified expression
  p <- ggplot(data, aes(x = time, y = value)) +
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
    scale_y_continuous(limits = c(1.8, 2.8)) +
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
  
  print(p)
  #ggsave(filename = paste0(dir$output, "/log_eq.pdf"), width = 8, height = 6, dpi = 300)
}
#------------------------------------------------------------------------------#
###.                             Setup for alg                               ###
#------------------------------------------------------------------------------#
# Initialize an empty list to store the B results
B <- list()

# Initialize counter
c <- 1

# B =
# Iterate over b_F, b_Z, and b_X ensuring they sum to 1
for (b_F in seq(0.01, 0.99, by = 0.01)) {
  for (b_Z in seq(0.01, 1 - b_F, by = 0.01)) {
    b_X <- 1 - b_F - b_Z  # Directly calculate b_X to satisfy the constraint
    if (b_X >= 0.01) {    # Ensure b_X meets the minimum step size
      # Add the current values of b_F, b_Z, and b_X as a list to result_list
      B[[c]] <- list(b_F = b_F, b_Z = b_Z, b_X = b_X)
      c <- c + 1
    }
  }
}

# number of control unit 
J <- i_max - 1
w <- rep(0, J)

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

## plot for Logrithmic Equi-confounding(Figure 1b)
logorithmic_equi_conf_plot_DGP(F= F, Y = Y)


#------------------------------------------------------------#
##                       Simulation                        ###
#------------------------------------------------------------#
# number of simulation 
n = 300 #300


initial <- 10
end_point <- 100
step_size <- 10

# Reinitialize outcome variables and parameters for simulation
F <- matrix(0, nrow = i_max, ncol = end_point)
Y <- matrix(0, nrow = i_max, ncol = s_max)


rho <- sort(runif(end_point, 0, 20))

phi <- matrix(runif(dr * end_point, min=0, max= dist_max), nrow = end_point, ncol = dr)
theta <- matrix(runif(du * end_point, min =0, max = dist_max), nrow = end_point, ncol = du)


# Define the range of t_max values
t_max_values <- seq(initial, end_point, by = step_size)


# This may take 24h!
results_list <- pblapply(t_max_values, function(t){
  print(t) 
  simulation_B_list_3DGP(
    time = t,
    B = B,
    eta_Z = 0.1, #0.02
    eta_X = 0.1, #0.05
    alg = TRUE
  )
})



#------------------------------------------------------------------------#
##                          Graph Plot                         ###
#------------------------------------------------------------------------#
# --------------------------  Bias plot ------------------------------- #
bias_plot(results_list)


# --------------------------     NSE plot      ------------------------ #
NSE_plot_ZXF(results_list)

NSE_plot_ZX(results_list)


#------------------------------------------------------------------------#
##                          Graph Plot                         ###
#------------------------------------------------------------------------#
# --------------------------  Bias plot ------------------------------- #
load(paste0(dir$output, "results_list0228.RData"))
results_list20 <- results_list

load(paste0(dir$output, "results_list0301.RData"))
results_list10 <- results_list

sc_result_list <- c(results_list10, results_list20)

load(paste0(dir$output, "results_list_DGP_linear0119.RData"))
linear_result_list <- results_list

load(paste0(dir$output, "results_list_DGP_log0119.RData"))
log_result_list <- results_list

combined_results_list <- lapply(1:length(t_max_values), function(i) {
  list(
    mean_bias_linear = linear_result_list[[i]]$mean_bias_linear,
    mean_bias_log = log_result_list[[i]]$mean_bias_log,
    mean_bias_sc = sc_result_list[[i]]$mean_bias_sc 
  )
})




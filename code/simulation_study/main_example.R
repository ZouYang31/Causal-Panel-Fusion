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
##---------------------   Data Generating Process   ----------------------##
set.seed(45) #45
#80 #99 #86
#95 good one. 
#45

# initialization 
i_max <- 31 #donor pool + 1 treated unit
t_max <- 20 # reference domain #20
s_max <- 5 # target domain


dr <- 3 #3 #5  # Number of observed covariates reference domain
dt <- 3        # Number of observed covariates target domain
#Note: dr is not necessarilly equal to dt. 
du <- 3  #3  # Number of unobserved factor loadings d_u


mu <- matrix(runif(i_max * du, min=0, max=1), nrow = i_max, ncol = du) # latent confounder
Z <- matrix(runif(i_max * dr, min=0, max=1), nrow = i_max, ncol = dr) # reference domain
X <- matrix(runif(i_max * dt, min=0, max=1), nrow = i_max, ncol = dt) # target domain

# Generate parameters
dist_max <- 10 #max value for the uniform distribution

# ------------ reference domain -----------------#
rho <- sort(runif(t_max, 0, 20))

phi <- matrix(runif(dr * t_max, min=0, max= dist_max), nrow = t_max, ncol = dr)
theta <- matrix(runif(du * t_max, min =0, max = dist_max), nrow = t_max, ncol = du)
epsilon <- matrix(rnorm(i_max * t_max, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = t_max)

# ------------ target domain -------------------# 
varrho <- sort(runif(s_max, 0, 10))

varphi <- matrix(runif(dt * s_max, min=0, max =10), nrow = s_max, ncol = dt)
vartheta <- matrix(runif(du * s_max, min =0, max = 10), nrow = s_max, ncol = du)
varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)

causal <- sort(runif(s_max,min=2,max=5))

# Potential outcomes
F <- matrix(0, nrow = i_max, ncol = t_max)
Y <- matrix(0, nrow = i_max, ncol = s_max)

# Calculate F_it^N for each combination of i and t
for (t in 1:t_max) {
  for (i in 1:i_max) {
    F[i, t] <- rho[t] + sum(phi[t, ] * Z[i, ]) + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
  }
}

for (t in 1:s_max) {
  for (i in 1:i_max) {
    Y[i, t] <- varrho[t] + sum(varphi[t, ] * X[i, ]) + sum(vartheta[t, ] * mu[i, ]) + varepsilon[i, t]
  }
}

# Print the resulting matrices
#print(F)
#print(Y)


#source(paste0(dir$code, "/0_function_simulation.R"))
source(paste0(dir$code, "/0_algorithm.R"))
source(paste0(dir$code, "/1_function_plot.R"))
#source(paste0(dir$code, "/2_function_optim_w.R"))

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
                      NSE_Z_baseline = NSE_Z_baseline, NSE_X_baseline = NSE_X_baseline, eta_Z=0.02, eta_X=0.02) 

b_list


# Print the resulting matrices
#print(F)
#print(Y)

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
Placebo_test(F = F, Y = Y, J = J, t_max = t_max, s_max = s_max, w = w)










##---------------           Data Analysis         --------------------##
##---------------      Equi-confounding Data Fusion      -------------##
# plot for linear equi-confounding(Figure 1a)
linear_equi_conf_plot(F = F, Y = Y)


## plot for Logrithmic Equi-confounding(Figure 1b)
logorithmic_equi_conf_plot(F= F, Y = Y)


#----------------------
## Linear Equi-confounding 
Y_1 = mean(Y[1, ] + causal - varepsilon[1, ])
Y_0_t = mean(Y[1, ] - varepsilon[1, ])
Y_0 = mean(Y[2:(J+1), ] - varepsilon[2:(J+1), ] )

F_1 = mean(F[1, ] - epsilon[1, ])
F_0 = mean(F[2:(J+1), ] - epsilon[2:(J+1), ])


estimate <- (Y_1 - F_1) - (Y_0 - F_0)
bias = abs(mean(causal) - estimate) 

print(paste("The value of the ETT is",estimate))
print(paste("The value of the bias is",bias))

#----------------------
## Logarithm Equi-confounding 

estimate <- Y_1 - (F_1/F_0) * Y_0
bias = abs(mean(causal) - estimate)

print(paste("The value of the ETT is",estimate))
print(paste("The value of the bias is",bias))




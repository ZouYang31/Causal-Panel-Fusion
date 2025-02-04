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


dr <- 3   # Number of observed covariates in reference domain
dt <- 3   # Number of observed covariates in target domain
#Note: dr is not necessarily equal to dt. 
du <- 3  #3  # Number of unobserved factor loadings

mu <- matrix(runif(i_max * du, min=0, max=1), nrow = i_max, ncol = du) # latent confounder
Z <- matrix(runif(i_max * dr, min=0, max=1), nrow = i_max, ncol = dr) # covariates in reference domain
X <- matrix(runif(i_max * dt, min=0, max=1), nrow = i_max, ncol = dt) # covariates in target domain

# Generate parameters
dist_max <- 10 # max value for the uniform distribution

# ------------ reference domain -----------------#
# rho <- sort(runif(t_max, 0, 20))

# phi <- matrix(runif(dr * t_max, min=0, max= dist_max), nrow = t_max, ncol = dr)
# theta <- matrix(runif(du * t_max, min =0, max = dist_max), nrow = t_max, ncol = du)
# epsilon <- matrix(rnorm(i_max * t_max, mean = 0, sd = sqrt(2)), nrow = i_max, ncol = t_max)

# ------------ target domain -------------------# 
varrho <- sort(runif(s_max, 0, 10))

varphi <- matrix(runif(dt * s_max, min=0, max =10), nrow = s_max, ncol = dt)
vartheta <- matrix(runif(du * s_max, min =0, max = 10), nrow = s_max, ncol = du)
varepsilon <- matrix(rnorm(i_max * s_max, mean = 0, sd = sqrt(1/2)), nrow = i_max, ncol = s_max)

causal <- sort(runif(s_max,min=2,max=5))

# Potential outcomes
# F <- matrix(0, nrow = i_max, ncol = t_max)
# Y <- matrix(0, nrow = i_max, ncol = s_max)
# 
# # Calculate F_it for each combination of i and t
# for (t in 1:t_max) {
#   for (i in 1:i_max) {
#     F[i, t] <- rho[t] + sum(phi[t, ] * Z[i, ]) + sum(theta[t, ] * mu[i, ]) + epsilon[i, t]
#   }
# }
# 
# for (t in 1:s_max) {
#   for (i in 1:i_max) {
#     Y[i, t] <- varrho[t] + sum(varphi[t, ] * X[i, ]) + sum(vartheta[t, ] * mu[i, ]) + varepsilon[i, t]
#   }
# }
# 
# # Treated Unit and Contro Unit
# F_treated <- F[1,]
# F_control <- F[-1,]
#print(Y)


source(paste0(dir$code, "/0_function_simulation.R"))
source(paste0(dir$code, "/0_algorithm.R"))
source(paste0(dir$code, "/1_function_plot.R"))


#------------------------------------------------------------------------------#
###.                             Setup for alg                               ###
#------------------------------------------------------------------------------#
# Initialize an empty list to store the B results
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


#------------------------------------------------------------#
##                       Simulation                        ###
#------------------------------------------------------------#
# number of simulation 
n = 200 #1000

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

results_list <- pblapply(t_max_values, function(t){
  print(t)
  simulation_B_list(
    time = t,
    B = B,
    eta_Z = 0.1,
    eta_X = 0.1,
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







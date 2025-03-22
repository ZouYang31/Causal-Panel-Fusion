########################################################################
# Causal Data Fusion: Application to COVID-19 Vaccination Rate Analysis
# 
#######################################################################

# -----------------------------------------------------------------
# Clean up work space and load or install necessary packages if necessary
rm(list=ls())
want <- c("dplyr", "ggplot2", "scales", "reshape2", "tidyr", "gurobi", "kernlab", "optimx", "knitr") # "limSolve"
need <- want[!(want %in% installed.packages()[,"Package"])]
if(length(need)) install.packages(need)
lapply(want, function(i) require(i, character.only = TRUE))
rm(want, need)

library(knitr)


# Load Required Libraries
library(dplyr)
library(ggplot2)
library(limSolve)
library(tidyr)
library(kernlab)

# Working directories
dir <- list()
dir$root <- getwd()
dir$output <- paste(dir$root, "/output", sep= "")
dir$code <- paste(dir$root, "/code", sep= "")
dir$data <- paste(dir$root, "/data", sep= "")


# Import the scripts
source(paste0(dir$code,"/function_data_preprocess.R"))
source(paste0(dir$code,"/0_algorithm.R"))
source(paste0(dir$code,"/1_function_plot.R"))
source(paste0(dir$code,"/2_function_optim_w.R"))
source(paste0(dir$code,"/function_graphs.R"))
source(paste0(dir$code,"/function_etts.R"))


#Load Data
covid_demo_yr <- read.csv(paste0(dir$data, "/covariates_20cities.csv"))
covid_vac_mo <- read.csv(paste0(dir$data, "/monthly_covid_vaccination_rate_20cities.csv"))

# Filter covariates data only for 2021
covid_demo_yr <- covid_demo_yr %>% filter(year >= 2021)
# Filter data for periods March 2021 to June 2022
covid_vac_mo <- covid_vac_mo %>% filter(year_month < "2022-07")


# -----------------------------
# Data Processing and Filtering
# -----------------------------

# Reorder Data Based on target unit.
df_reordered_hispanic <- reorder_unit_first(covid_demo_yr, "Chelsea-hispanic")
df_reordered_black <- reorder_unit_first(covid_demo_yr, "Chelsea-black")

# Extract Subsets for targeted sub-population.
df_black <- df_reordered_black %>% filter(grepl("-black", unit))
df_hispanic <- df_reordered_hispanic  %>% filter(grepl("-hispanic", unit))

# Normalize the reference domain
df_black <- normalize_columns(
  df_black, 
  c("median_income","median_age", "proportion", "X65p_proportion")
)

# Normalize the target domain
df_hispanic <- normalize_columns(
  df_hispanic, 
  c("median_income","median_age", "proportion", "X65p_proportion")
)


# -------------------------------------
# Process Covariates for Reference (Z) and Target (X) Domains
# -------------------------------------

# Covariates in referece domain
Z <- df_black[df_black$year == 2021, c("median_income", "proportion", "median_age", "X65p_proportion"), drop = FALSE]
colnames(Z) <- NULL
Z <- as.matrix(sapply(Z, as.numeric))


# Covariates in target domain
X <- df_hispanic[df_hispanic$year == 2021, c("median_income", "proportion", "median_age", "X65p_proportion"), drop = FALSE]
colnames(X) <- NULL
X <- as.matrix(sapply(X, as.numeric))

num_city <- 20   # Total number of cities
J <- num_city - 1 # Control units
w <- rep(0, J)    # Initialize weight vector

# -------------------------------------
# Reshape COVID Data by Month
# -------------------------------------
covid_reshape_mo <- covid_vac_mo %>%
  pivot_wider(names_from = year_month, values_from = fully_vac_rate)

# -------------------------------------
# Extract & Reorder Data for Black and Hispanic Populations
# -------------------------------------
df_black_mo <- covid_reshape_mo %>%
  filter(grepl("-black", unit)) %>%
  reorder_unit_first("Chelsea-black")  # Uses pre-defined function

df_hispanic_mo <- covid_reshape_mo %>%
  filter(grepl("-hispanic", unit)) %>%
  reorder_unit_first("Chelsea-hispanic")  # Uses pre-defined function

F <- convert_to_matrix(df_black_mo)
Y <- convert_to_matrix(df_hispanic_mo)

# -------------------------------------
# Extract Treated and Control Groups Returning Values
# -------------------------------------
Y_treated <- as.numeric(Y[1, ])
Y_control <- apply(Y[2:(J + 1), ], 2, as.numeric)

F_treated <- as.numeric(F[1, ])
F_control <- apply(F[2:(J + 1), ], 2, as.numeric)

Z_treated <- as.numeric(Z[1, ])
Z_control <- apply(Z[2:(J + 1), ], 2, as.numeric)

X_treated <- as.numeric(X[1, ])
X_control <- apply(X[2:(J + 1), ], 2, as.numeric)

# -------------------------------------
# Define Constants
# -------------------------------------

# Number of unique time points in vaccination data
t_max <- ncol(df_black_mo) -1 
s_max <- ncol(df_hispanic_mo) -1 


# Number of covariates in reference and target domains
dr <- length(Z_treated)   # Reference domain covariates
dt <- length(X_treated)  # Target domain covariates

#------------------------------------------------------------------------------#
###.                             Setup for alg                               ###
#------------------------------------------------------------------------------#
# Create the B list.
B <- generate_b_list(step = 0.01, min_value = 0.01)


# Calculate the baseline for NSE_Z and NSE_X
result_X <- optimize_w_ipop(
  F_treated = F_treated,
  F_control = F_control,
  X = X,
  Z = Z,
  t_max = t_max,
  dr = dr,
  dt = dt,
  i_max = num_city,
  target = "X",
)

wX <- result_X$weights
NSE_X_baseline <- NSE_x(wX, F= F, X= X, Z= Z, t_max, dr, dt, i_max=num_city, target = "X")
print(NSE_X_baseline)


result_Z <- optimize_w_ipop(
  F_treated = F_treated,
  F_control = F_control,
  X = X,
  Z = Z,
  t_max = t_max,
  dr = dr,
  dt = dt,
  i_max = num_city,
  target = "Z",
)

wZ <- result_Z$weights
NSE_Z_baseline <- NSE_x(wZ,  F = F, X = X, Z = Z, t_max, dr, dt, i_max=num_city, target = "Z")
print(NSE_Z_baseline)


# Optimize Best B 
b_list <- find_best_B(
  B = B,
  F_treated = F_treated,
  F_control = F_control,
  X_treated = X_treated, 
  X_control = X_control, 
  Z_treated = Z_treated, 
  Z_control = Z_control, 
  t_max = t_max,
  dr = dr,
  dt = dt,
  i_max = num_city,
  NSE_Z_baseline = NSE_Z_baseline,
  NSE_X_baseline = NSE_X_baseline,
  eta_Z = 0.1, 
  eta_X = 0.1  
)


# Solve for weight vector (w) using the optimal B values
w <- b_list$best_w_star

# Display results
cat("Optimized weight vector (w):\n", w, "\n")
cat("Sum of weights:", sum(w), "\n")
cat("Number of positive weights (w > 0):", sum(w > 0), "\n")


# ------------------------------------------------------------------------------#
# Generate Graphs for Reference and Target Domains
# ------------------------------------------------------------------------------#
synthetic_reference_plot_real_data(F = F, w, t_max, i_max=num_city)
synthetic_target_plot_real_data(Y = Y, w, s_max = t_max, i_max=num_city)

# ------------------------------------------------------------------------------#
# Compute Estimated Treatment Effect (ETT) for Synthetic Control
# ------------------------------------------------------------------------------#

# Calculate the Estimated Treatment Effect (ETT)
ett_synthetic <- synthetic_ett(Y_treated, Y = Y, w, J)

# Display the result
cat("Estimated Treatment Effect (ETT) using the synthetic control data fusion method:", 
    round(mean(ett_synthetic), 5), "\n")

# ---------------------------------------------------------------------------------#
# Compute Estimated Treatment Effect (ETT) for Equi-Confounding Logarithmic and Linear
# ----------------------------------------------------------------------------------#
# Calculate the Estimated Treatment Effect (ETT)
ett_linear <- linear_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
cat("Estimated Treatment Effect (ETT) using linear equi-confounding:", round(ett_linear, 5), "\n")

ett_log <- log_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
cat("Estimated Treatment Effect (ETT) using logarithm equi-confounding:", round(ett_log, 5), "\n")

# ------------------------------------------------------------------------------#
# Generate Graphs for Linear Equi Confounding
# ------------------------------------------------------------------------------#
linear_equi_conf_plot_real_data(Y_treated, Y_control, F_treated, F_control)

# ------------------------------------------------------------------------------#
# Generate Covid Vaccination Rate Predictor Means 
# ------------------------------------------------------------------------------#
# Vaccination Rate Predictor Means
# Reference domain
synth_median_income   <- sum(df_black[df_black$year == 2021, "median_income"][-1] * w)
synth_proportion      <- sum(df_black[df_black$year == 2021, "proportion"][-1] * w)
synth_median_age      <- sum(df_black[df_black$year == 2021, "median_age"][-1] * w)
synth_65proportion    <- sum(df_black[df_black$year == 2021, "X65p_proportion"][-1] * w)


Chelsea_median_income <- df_black[df_black$year == 2021, "median_income"][1]
Chelsea_proportion <- df_black[df_black$year == 2021, "proportion"][1]
Chelsea_median_age <- df_black[df_black$year == 2021, "median_age"][1]
Chelsea_65proportion <- df_black[df_black$year == 2021, "X65p_proportion"][1]

avg_median_income   <- mean(df_black[df_black$year == 2021, "median_income"][-1])
avg_proportion      <- mean(df_black[df_black$year == 2021, "proportion"][-1])
avg_median_age      <- mean(df_black[df_black$year == 2021, "median_age"][-1])
avg_65proportion    <- mean(df_black[df_black$year == 2021, "X65p_proportion"][-1])

# Target domain
synth_median_income_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "median_income"][-1] * w)
synth_proportion_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "proportion"][-1] * w)
synth_median_age_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "median_age"][-1] * w)
synth_65proportion_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "X65p_proportion"][-1] * w)


Chelsea_median_income_hispanic <- df_hispanic[df_hispanic$year == 2021, "median_income"][1]
Chelsea_proportion_hispanic <- df_hispanic[df_hispanic$year == 2021, "proportion"][1]
Chelsea_median_age_hispanic <- df_hispanic[df_hispanic$year == 2021, "median_age"][1]
Chelsea_65proportion_hispanic <- df_hispanic[df_hispanic$year == 2021, "X65p_proportion"][1]


avg_median_income_hispanic <- mean(df_hispanic[df_hispanic$year == 2021, "median_income"][-1])
avg_proportion_hispanic <- mean(df_hispanic[df_hispanic$year == 2021, "proportion"][-1])
avg_median_age_hispanic <- mean(df_hispanic[df_hispanic$year == 2021, "median_age"][-1])
avg_65proportion_hispanic <- mean(df_hispanic[df_hispanic$year == 2021, "X65p_proportion"][-1])

# ------------------------------------------------------------------------------#
# Generate Graph for Placebo Test
# ------------------------------------------------------------------------------#

Placebo_test_data(F = F, Y = Y, J = J, t_max = t_max, s_max = s_max, 
             w = w, X = X, Z = Z,
             dr = dr,
             dt = dt,
             i_max = num_city,
             eta_Z = 0.1, eta_X = 0.1)


# ------------------------------------------------------------------------------#
# Sensitive Analysis
#
# ------------------------------------------------------------------------------#
eta_values <- c(0.05, 0.1, 0.15, 0.2)
eta_combinations <- expand.grid(eta_X = eta_values, eta_Z = eta_values)

# Initialize a dataframe to store results
results <- eta_combinations %>%
  mutate(
    diff_median_income_black = NA,
    diff_proportion_black = NA,
    diff_median_age_black = NA,
    diff_65proportion_black = NA,
    diff_vac_black = NA, 
    
    diff_median_income_hispanic = NA,
    diff_proportion_hispanic = NA,
    diff_median_age_hispanic = NA,
    diff_65proportion_hispanic = NA,
    diff_vac_hispanic = NA,
    
    ett_synthetic = NA
  )



# Loop over all 16 combinations of eta_X and eta_Z
for (i in seq_len(nrow(results))) {
  eta_X <- results$eta_X[i]
  eta_Z <- results$eta_Z[i]
  
  # Call the function with different eta_X and eta_Z values
  b_list <- find_best_B(
    B,
    F_treated = F_treated,
    F_control = F_control,
    X_treated = as.matrix(X_treated),
    X_control = as.matrix(X_control),
    Z_treated = as.matrix(Z_treated),
    Z_control = as.matrix(Z_control),
    t_max = t_max,
    dr = dr,
    dt = dt,
    i_max = num_city,
    NSE_Z_baseline = NSE_Z_baseline,
    NSE_X_baseline = NSE_X_baseline,
    eta_Z = eta_Z,
    eta_X = eta_X
  )
  
  # Solve for weight vector (w) using the optimal B values
  w <- b_list$best_w_star
  
  # Sensitivity matrix
  synth_median_income <- sum(df_black[df_black$year == 2021, "median_income"][-1] * w)
  synth_proportion <- sum(df_black[df_black$year == 2021, "proportion"][-1] * w)
  synth_median_age <- sum(df_black[df_black$year == 2021, "median_age"][-1] * w)
  synth_65proportion <- sum(df_black[df_black$year == 2021, "X65p_proportion"][-1] * w)
  synth_vac_black <- sum(df_black[df_black$year == 2021, "covid_fully_vac_rate"][-1] * w)
  
  
  Chelsea_median_income <- df_black[df_black$year == 2021, "median_income"][1]
  Chelsea_proportion <- df_black[df_black$year == 2021, "proportion"][1]
  Chelsea_median_age <- df_black[df_black$year == 2021, "median_age"][1]
  Chelsea_65proportion <- df_black[df_black$year == 2021, "X65p_proportion"][1]
  Chelsea_vac_black <- df_black[df_black$year == 2021, "covid_fully_vac_rate"][1]
  
  diff_median_income_black <- abs(Chelsea_median_income - synth_median_income)
  diff_proportion_black <- abs(Chelsea_proportion - synth_proportion)
  diff_median_age_black <- abs(Chelsea_median_age - synth_median_age)
  diff_65proportion_black <- abs(Chelsea_65proportion - synth_65proportion)
  diff_vac_black <- abs(Chelsea_vac_black - synth_vac_black)
  
  
  synth_median_income_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "median_income"][-1] * w)
  synth_proportion_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "proportion"][-1] * w)
  synth_median_age_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "median_age"][-1] * w)
  synth_65proportion_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "X65p_proportion"][-1] * w)
  synth_vac_hispanic <- sum(df_hispanic[df_hispanic$year == 2021, "covid_fully_vac_rate"][-1] * w)
  
  Chelsea_median_income_hispanic <- df_hispanic[df_hispanic$year == 2021, "median_income"][1]
  Chelsea_proportion_hispanic <- df_hispanic[df_hispanic$year == 2021, "proportion"][1]
  Chelsea_median_age_hispanic <- df_hispanic[df_hispanic$year == 2021, "median_age"][1]
  Chelsea_65proportion_hispanic <- df_hispanic[df_hispanic$year == 2021, "X65p_proportion"][1]
  Chelsea_vac_hispanic <- df_hispanic[df_hispanic$year == 2021, "covid_fully_vac_rate"][1]
  
  diff_median_income_hispanic <- abs(Chelsea_median_income_hispanic - synth_median_income_hispanic)
  diff_proportion_hispanic <- abs(Chelsea_proportion_hispanic - synth_proportion_hispanic)
  diff_median_age_hispanic <- abs(Chelsea_median_age_hispanic - synth_median_age_hispanic)
  diff_65proportion_hispanic <- abs(Chelsea_65proportion_hispanic - synth_65proportion_hispanic)
  diff_vac_hispanic <- abs(Chelsea_vac_hispanic - synth_vac_hispanic)
  
  ett_synthetic <- synthetic_ett(Y_treated, Y = combined_Y, w, J)
  
  # Store results
  results$diff_median_income_black[i] <- diff_median_income_black
  results$diff_proportion_black[i] <- diff_proportion_black
  results$diff_median_age_black[i] <- diff_median_age_black
  results$diff_65proportion_black[i] <- diff_65proportion_black
  results$diff_vac_black[i] <- round(diff_vac_black, 3)
  
  results$diff_median_income_hispanic[i] <- diff_median_income_hispanic
  results$diff_proportion_hispanic[i] <- diff_proportion_hispanic
  results$diff_median_age_hispanic[i] <- diff_median_age_hispanic
  results$diff_65proportion_hispanic[i] <- diff_65proportion_hispanic
  results$diff_vac_hispanic[i] <- round(diff_vac_hispanic, 3)
  
  results$ett_synthetic[i] <- round(ett_synthetic, 3)
}



# Reshape data for heatmap plotting
results_long <- melt(results, id.vars = c("eta_X", "eta_Z"))

# List of metrics to plot
metrics <- c("diff_median_income_black", "diff_proportion_black", "diff_median_age_black", 
             "diff_65proportion_black", "diff_vac_black",
             "diff_median_income_hispanic", "diff_proportion_hispanic", 
             "diff_median_age_hispanic", "diff_65proportion_hispanic", "diff_vac_hispanic",
             "ett_synthetic")

# Loop over each metric and plot individually
for (metric in metrics) {
  p <- ggplot(results_long %>% filter(variable == metric), aes(x = eta_X, y = eta_Z, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 3)), size = 5, color = "black") +  
    scale_fill_gradient2(#low = "#D6EAF8", mid = "skyblue", high = "blue", 
                         low = "#D6EAF8", mid ="#ADD8E6", high = "skyblue", 
                         midpoint = median(results_long$value, na.rm = TRUE),
                         name = "Causal effect") +
    labs(title = paste(""),
         x = expression(eta[X]),
         y = expression(eta[Z])) +
    theme_minimal()
  
  print(p)
  ggsave(filename = paste0(dir$output, metric, "_heatmap.png"), width = 6, height = 5) 
}



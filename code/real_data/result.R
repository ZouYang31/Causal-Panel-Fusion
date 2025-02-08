# Result

source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/function_preprocess.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/0_algorithm.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/1_function_plot.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/2_function_optim_w.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/function_graphs.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/function_etts.R")
source("/Users/leeseunghee/Dropbox/Vaccine_Causal/Code/function_sensitive_analysis.R")

library(knitr)
# -----------------------------
# Data Processing and Filtering
# -----------------------------
covid_demo_yr <- covid_demo_yr %>% filter(year >= 2021)
covid_vac_mo <- covid_vac_mo %>% filter(year_month < "2022-07")

#For the future usage
covid_demo_copy <- covid_demo_yr

# -----------------------------
# Data Processing and Filtering
# -----------------------------

# Apply Normalization to Selected Columns
covid_demo_yr <- normalize_columns(
  covid_demo_yr, 
  c("median_income","median_age", "proportion", "X65p_proportion")
)

# Reorder Data Based on Specific Units
df_covid_reordered <- reorder_unit_first(covid_demo_yr, "Chelsea-hispanic")
df_covid_reordered_black <- reorder_unit_first(covid_demo_yr, "Chelsea-black")

# Extract Subsets for Specific Demographic Groups
df_black <- df_covid_reordered_black %>% filter(grepl("-black", unit))
df_hispanic <- df_covid_reordered %>% filter(grepl("-hispanic", unit))

# Filter Out Specific Groups from Backup Dataset
df_black_filtered<- covid_demo_copy %>%
  filter(unit != "Chelsea-black" & !grepl("-asian|-white|-hispanic", unit))

df_hispanic_filtered <- covid_demo_copy %>%
  filter(unit != "Chelsea-hispanic" & !grepl("-asian|-white|-black", unit))

# -------------------------------------
# Aggregate Vaccination Rates for Black and Hispanic Populations
# -------------------------------------
combined_black <- combine_yearly_data(df_black, 2021, 2022, "covid_fully_vac_rate")
combined_hispanic <- combine_yearly_data(df_hispanic, 2021, 2022, "covid_fully_vac_rate")

# -------------------------------------
# Process Covariates for Reference (Z) and Target (X) Domains
# -------------------------------------
combined_Z <- combine_yearly_data(df_black, 2021, 2022, c("proportion", "median_income"))
combined_X <- combine_yearly_data(df_hispanic, 2021, 2022, c("proportion", "median_income"))


num_city <- 20   # Total number of cities
J <- num_city - 1 # Number of comparisons (excluding reference city)
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

combined_F <- convert_to_numeric_matrix(df_black_mo)
combined_Y <- convert_to_numeric_matrix(df_hispanic_mo)

# -------------------------------------
# Extract Treated and Control Groups Returning Values
# -------------------------------------
Y_treated <- as.numeric(combined_Y[1, ])
Y_control <- apply(combined_Y[2:(J + 1), ], 2, as.numeric)

F_treated <- as.numeric(combined_F[1, ])
F_control <- apply(combined_F[2:(J + 1), ], 2, as.numeric)

Z_treated <- as.numeric(combined_Z[1, ])
Z_control <- apply(combined_Z[2:(J + 1), ], 2, as.numeric)

X_treated <- as.numeric(combined_X[1, ])
X_control <- apply(combined_X[2:(J + 1), ], 2, as.numeric)

# -------------------------------------
# Define Constants
# -------------------------------------

# Number of unique time points in vaccination data
t_max <- ncol(df_black_mo) -1 
s_max <- ncol(df_hispanic_mo) -1 


# Number of covariates in reference and target domains
dr <- 4  # Reference domain covariates
dt <- 4  # Target domain covariates

#------------------------------------------------------------------------------#
###.                             Setup for alg                               ###
#------------------------------------------------------------------------------#
# Initialize an empty list to store the results
B <- list()

# Initialize counter
c <- 1

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

# Convert all columns to numeric matrices
combined_Z_num <- as.matrix(sapply(combined_Z, as.numeric))
combined_X_num <- as.matrix(sapply(combined_X, as.numeric))

# Calculate the baseline for NSE_Z and NSE_X
result_X <- optimize_w_ipop(
  F_treated = F_treated,
  F_control = F_control,
  X = combined_X_num,
  Z = combined_Z_num,
  t_max = t_max,
  dr = dr,
  dt = dt,
  i_max = num_city,
  target = "X",
)

wX <- result_X$weights
NSE_wX <- NSE_x(wX, F=combined_F, X=combined_X_num, Z=combined_Z_num, t_max, dr, dt, i_max=num_city, target = "X")
print(NSE_wX)
NSE_X_baseline = NSE_wX

result_Z <- optimize_w_ipop(
  F_treated = F_treated,
  F_control = F_control,
  X = combined_X_num,
  Z = combined_Z_num,
  t_max = t_max,
  dr = dr,
  dt = dt,
  i_max = num_city,
  target = "Z",
)

wZ <- result_Z$weights
NSE_wZ <- NSE_x(wZ,  F=combined_F, X=combined_X_num, Z=combined_Z_num, t_max, dr, dt, i_max=num_city, target = "Z")
NSE_Z_baseline <- NSE_wZ
print(NSE_wZ)

# Optimize Best B 
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
  eta_Z = 0.02,
  eta_X = 0.05
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
synthetic_reference_plot_real_data(F=combined_F, w, t_max, i_max=num_city)
synthetic_target_plot_real_data(Y=combined_Y, w, s_max, i_max=num_city)

# ------------------------------------------------------------------------------#
# Compute Estimated Treatment Effect (ETT) for Synthetic Control
# ------------------------------------------------------------------------------#

# Calculate the Estimated Treatment Effect (ETT)
ett_synthetic <-synthetic_ett(Y_treated, Y=combined_Y, w, J)

# Display the result
cat("Estimated Treatment Effect (ETT) using the synthetic control data fusion method:", 
    round(mean(ett_synthetic), 5), "\n")

# ---------------------------------------------------------------------------------#
# Compute Estimated Treatment Effect (ETT) for Equi-Confounding Logarithmic and Linear
# ----------------------------------------------------------------------------------#
# Calculate the Estimated Treatment Effect (ETT)
ett_linear <- linear_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)
# Display the result
cat("Estimated Treatment Effect (ETT) using linear equi-confounding:", round(ett_linear, 5), "\n")

ett_log <- log_equi_confounding_ett(Y_treated, Y_control, F_treated, F_control)

# Display the result
cat("Estimated Treatment Effect (ETT) using logarithm equi-confounding:", round(ett_log, 5), "\n")

# ------------------------------------------------------------------------------#
# Generate Graphs for Linear Equi Confounding
# ------------------------------------------------------------------------------#
linear_equi_conf_plot_real_data(Y_treated, Y_control, F_treated, F_control)

# ------------------------------------------------------------------------------#
# Generate Covid Vaccination Rate Predictor Means 
# ------------------------------------------------------------------------------#

covariate_columns <- c("median_income", "proportion", "median_age", "X65p_proportion")

# Process data for Black population
black_results <- process_group_data(df_black_filtered, "Chelsea-black", w)

# Process data for Hispanic population
hispanic_results <- process_group_data(df_hispanic_filtered, "Chelsea-hispanic", w)

# ------------------------------------------------------------------------------#
# Generate Graph for Placebo Test
# ------------------------------------------------------------------------------#
# Placebo_test
placebo_result <- Placebo_test_v2(F = combined_F, Y = combined_Y, J = J, t_max = t_max, s_max = s_max, w = w, X = combined_X_num,
                                  Z = combined_Z_num,
                                  dr = dr,
                                  dt = dt,
                                  i_max = num_city,
                                  B = B)

# ------------------------------------------------------------------------------#
# Generate Sensitive Anlysis
# ------------------------------------------------------------------------------#
sensitivity_results <- sensitivity_analysis(
  B = B, 
  F_treated = F_treated, 
  F_control = F_control, 
  X_treated = X_treated, 
  X_control = X_control, 
  Z_treated = Z_treated, 
  Z_control = Z_control, 
  F = combined_F, 
  X = combined_X_num, 
  Z = combined_Z_num, 
  t_max = t_max, 
  dr = dr, 
  dt = dt, 
  num_city = num_city, 
  NSE_Z_baseline = NSE_Z_baseline, 
  NSE_X_baseline = NSE_X_baseline
)

# Display results
knitr::kable(sensitivity_results$EST_NSEs, "latex")
knitr::kable(sensitivity_results$B_lists, "latex")


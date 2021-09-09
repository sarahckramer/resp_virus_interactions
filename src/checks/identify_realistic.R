# ---------------------------------------------------------------------------------------------------------------------
# Explore which combinations of parameters lead to "realistic" outbreaks, and how this changes with rho1/rho2 values
# Do LHS and plot pairs plot according to whether outbreaks are realistic or not
# Also see how realistic values vary based on theta_lambda1/theta_lambda2
# Finally, what about when we hold R10/R20 to "realistic" values (can allow R120 to vary)?
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
# library(tidyverse)
# library(pomp)
library(tgp)
# library(gridExtra)

# Set seed:
set.seed(1902574195)

# Set global parameters:
n_lhs <- 500
n_sim <- 5

# Load functions:
source('src/functions/functions_sim_dat.R')

# ---------------------------------------------------------------------------------------------------------------------

# Prep pomp object

# Set necessary parameter values:
vir1 <- 'flu_A' # 'flu_A', 'flu_B'
vir2 <- 'rsv'
yr <- 2006

debug_bool <- FALSE

Ri_max1 <- 2.0
Ri_max2 <- 2.0
delta_min <- 7 / 60.0

# Load pomp object:
source('src/resp_interaction_model.R')

# ---------------------------------------------------------------------------------------------------------------------

# First: Get realistic parameter values with no interaction and "original" rho values

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 0.75, 0.75, 0.75))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = TRUE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
parms_matrix_ORIG <- parms_matrix
rm(parm_sets)

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]])

# Clean up:
rm(sim_res)

# ---------------------------------------------------------------------------------------------------------------------

# Run similar analysis at various values of rho1/rho2

# Set rho values of interest:
rho1_vals <- c(0.25, 0.50)
rho2_vals <- c(0.05, 0.15, 0.35)

# Loop through reporting rates and run:
for (r1 in rho1_vals) {
  for (r2 in rho2_vals) {
    
    # Which rho values currently?:
    print(paste0(r1, '_', r2))
    
    # Update parms_matrix:
    parms_matrix <- parms_matrix_ORIG
    parms_matrix['rho1', ] <- r1
    parms_matrix['rho2', ] <- r2
    
    # Run simulations and identify realistic:
    sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)
    sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)
    
    # Plot results:
    plot_res(sim_res[[2]], plot_title = paste0(r1, '_', r2))
    
    # Clean up:
    rm(sim_res)
    
  }
}

# ---------------------------------------------------------------------------------------------------------------------

# Update rho2:
parms_matrix_ORIG['rho2', ] <- 0.15

# ---------------------------------------------------------------------------------------------------------------------

# Run similar analysis at various levels of theta_lambda1/theta_lambda2

# Set interaction strengths:
int_vals <- c(0, 0.25, 0.5, 0.75)

# Loop through theta_lambda2 values and run:
for (int in int_vals) {
  
  # Which int values currently?:
  print(paste0('theta_lambda2=', int))
  
  # Update parms_matrix:
  parms_matrix <- parms_matrix_ORIG
  parms_matrix['theta_lambda2', ] <- int
  
  # Run simulations and identify realistic:
  sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)
  sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)
  
  # Plot results:
  plot_res(sim_res[[2]], plot_title = paste0('theta_lambda2=', int))
  
  # Clean up:
  rm(sim_res)
  
}

# Loop through theta_lambda1 values and run:
for (int in int_vals) {
  
  # Which int values currently?:
  print(paste0('theta_lambda1=', int))
  
  # Update parms_matrix:
  parms_matrix <- parms_matrix_ORIG
  parms_matrix['theta_lambda1', ] <- int
  
  # Run simulations and identify realistic:
  sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)
  sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)
  
  # Plot results:
  plot_res(sim_res[[2]], plot_title = paste0('theta_lambda1=', int))
  
  # Clean up:
  rm(sim_res)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Run similar analysis using "realistic" values of R10/R20/R120

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Adjust values of R20/R120:
parms_matrix['R10', ] <- 0.50
parms_matrix['R20', ] <- 0.35

# Adjust rho2:
parms_matrix['rho2', ] <- 0.15

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]], plot_title = 'Realistic R', incl_R2 = FALSE)

# Clean up:
rm(sim_res)

# Again, but allow rho2 to "fit":

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'rho2')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 1.0))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Adjust values of R20/R120:
parms_matrix['R10', ] <- 0.50
parms_matrix['R20', ] <- 0.35

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]], plot_title = 'Realistic R', incl_R2 = FALSE)

# Clean up:
rm(sim_res)


# ---------------------------------------------------------------------------------------------------------------------

# Run similar analysis setting R20/R120 to zero

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 0.75))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Adjust rho2:
parms_matrix['rho2', ] <- 0.15

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]], plot_title = 'R20=R120=0', incl_R2 = FALSE)

# Clean up:
rm(sim_res)

# Again, but allow rho2 to "fit":

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'rho2')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 0.75, 1.0))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]], plot_title = 'R20=R120=0', incl_R2 = FALSE)

# Clean up:
rm(sim_res)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

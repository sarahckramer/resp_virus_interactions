# ---------------------------------------------------------------------------------------------------------------------
# Explore how changing the values of interaction parameters influences model likelihoods/simulations
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(gridExtra)
library(testthat)

# Set whether plots should be stored:
save_plots <- TRUE

# Set model parameters:
prof_lik <- FALSE

# ---------------------------------------------------------------------------------------------------------------------

# Functions

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, true_estpars, run_name) {
  
  # # Compile estpars:
  # true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full = list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get parameter estimates and log-likelihoods:
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
  no_best <- max(no_best, 50)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Return formatted results:
  return(pars_top)
  
}

# Function to generate slices/grids and calculate log-likelihoods:
calculate_slices_and_grid <- function(pars_top, vir1, shared_estpars) {
  
  # Set estpars:
  estpars <- names(pars_top %>%
                     select(-loglik))
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get MLE:
  mle <- setNames(object = as.numeric(pars_top[1, 1:length(estpars)]),
                  nm = estpars)
  
  # Get slices and calculate ll:
  slices <- slice_design(center = mle,
                         theta_lambda1 = seq(from = 0, to = 1.0, by = 0.025),
                         theta_lambda2 = seq(from = 0, to = 1.0, by = 0.025),
                         delta = seq(from = 0.01, to = 1.4, by = 0.025)) %>%
    mutate(ll = NA)
  
  for (ix in 1:nrow(slices)) {
    x0 <- slices[ix, 1:length(estpars)]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    slices$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Get grid and calculate ll:
  theta_lambda1_vals <- seq(from = 0, to = 1.0, by = 0.025)
  delta_vals <- seq(from = 0.01, to = 1.4, by = 0.025)
  
  grid <- expand_grid(theta_lambda1_vals, delta_vals) %>%
    rename('theta_lambda1' = 'theta_lambda1_vals',
           'delta' = 'delta_vals') %>%
    mutate(ll = NA)
  
  x0 <- mle
  for (ix in 1:nrow(grid)) {
    x0[c('theta_lambda1', 'delta')] <- grid[ix, ]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    grid$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Return data frames:
  return(list(slices, grid))
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1', 'theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in results (round 2):
res_h1 <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci/',
                                       true_estpars = true_estpars,
                                       run_name = 'flu_H1') %>%
  mutate(theta_lambda2 = 1.0) %>%
  select(rho1:theta_lambda1, theta_lambda2, `s13-14_Ri1`:loglik)
res_b <- load_and_format_mega_results(filepath = 'results/round2_flu_b_round1ci/',
                                      true_estpars = true_estpars,
                                      run_name = 'flu_B') %>%
  mutate(theta_lambda2 = 1.0) %>%
  select(rho1:theta_lambda1, theta_lambda2, `s13-14_Ri1`:loglik)
res_h1_BOTH <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci_thetalambda1_thetalambda2/',
                                            true_estpars = true_estpars,
                                            run_name = 'flu_H1_BOTH')

res_list <- list(res_h1, res_b, res_h1_BOTH)
names(res_list) <- c('H1', 'B', 'H1_both')
rm(res_h1, res_b, res_h1_BOTH)

# ---------------------------------------------------------------------------------------------------------------------

# Slice/grid likelihoods
# Slices: theta_lambda1, theta_lambda2, delta
# Grid: theta_lambda1, delta

# Loop through results and get slices:
slice_list = grid_list = vector('list', length = length(res_list))
for (i in 1:length(res_list)) {
  
  # Get vir1:
  if (str_detect(names(res_list)[i], 'H1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Calculate ll over slices:
  res_temp <- calculate_slices_and_grid(res_list[[i]], vir1, shared_estpars)
  
  # Store slices:
  slice_list[[i]] <- res_temp[[1]]
  grid_list[[i]] <- res_temp[[2]]
  
}
rm(i, vir1, res_temp, dat_pomp, hk_dat, obj_fun_list, po_list, seasons, yr, yr_index)

# Save plots to file?:
if (save_plots) {
  pdf('results/plots/explore_interaction_param_impact_HK_likelihood.pdf', width = 12, height = 9)
}

# Plot slices:
par(mfrow = c(3, 3))
for (i in 1:length(res_list)) {
  slice_temp <- slice_list[[i]]
  
  for (par in c('theta_lambda1', 'theta_lambda2', 'delta')) {
    slices_cur <- filter(slice_temp, slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = names(res_list)[i])
  }
}
rm(i, slice_temp, par, slices_cur)

# Plot grids:
for (i in 1:length(res_list)) {
  grid_temp <- grid_list[[i]]
  
  p_temp <- ggplot(data = grid_temp, aes(x = theta_lambda1, y = delta, fill = ll)) +
    geom_tile() + theme_classic() +
    scale_fill_viridis() +
    scale_x_continuous(n.breaks = 10, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 20, expand = c(0, 0)) +
    labs(fill = 'LL')
  print(p_temp)
}

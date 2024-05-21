# ---------------------------------------------------------------------------------------------------------------------
# Format and plot results from "round 2" of trajectory matching (using mega-likelihood)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Data from Canada?:
fit_canada <- FALSE

# Data from Germany?:
fit_germany <- FALSE

# Sensitivity analysis?:
if (fit_canada | fit_germany) {
  sens <- 'sinusoidal_forcing'
}

# Set directory where final results from round2 fits are stored:
if (fit_canada) {
  res_dir <- 'results/round2_fit/sens/canada/round2_3_flu/'
  res_dir_round1 <- 'results/round2_fit/sens/canada/round1_fitsharedFALSE/'
} else if (fit_germany) {
  res_dir <- 'results/round2_fit/sens/germany/round2_2_flu/'
  res_dir_round1 <- 'results/round2_fit/sens/germany/round1_fitsharedFALSE/'
} else {
  res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'
  res_dir_round1 <- 'results/round1_fitsharedFALSE/'
}

# Check that directory for storing plots exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/plots/')) {
  dir.create('results/plots')
}

# Get/format date (for saving results):
date <- format(Sys.Date(), '%d%m%y')

# ---------------------------------------------------------------------------------------------------------------------

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, shared_estpars, unit_estpars, run_name) {
  
  # Compile estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full = list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get parameter estimates and log-likelihoods:
  if (str_detect(res_files[1], 'PARALLEL')) {
    
    res_full <- do.call('c', res_full)
    num_errors <- length(which(res_full == 'error'))
    
    if (num_errors > 0) {
      res_full <- res_full[-which(res_full == 'error')]
    }
    
  }
  
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  if (str_detect(res_files[1], 'PARALLEL')) {
    expect_true(nrow(pars_df) == (length(res_files) * 25) - num_errors)
  } else {
    expect_true(nrow(pars_df) == length(res_files))
  }
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Check whether estimations stopped due to tolerance or time:
  print(table(pars_df$message))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  # expect_equal(df_use, 54)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  # no_best <- max(no_best, 20)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Format for combining with other results:
  pars_top_LONG <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 6),
           param = str_sub(param, 8, str_length(param))) %>%
    pivot_wider(names_from = param, values_from = 'value') %>%
    pivot_longer(-c(loglik, year), names_to = 'param') %>%
    select(year, param:value, loglik) %>%
    mutate(method = run_name)
  
  # Explore correlations:
  pars_corr <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 6),
           param = str_sub(param, 8, str_length(param))) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    # select(-year) %>%
    select(all_of(true_estpars), loglik) %>%
    mutate(method = run_name)
  
  # Return formatted results:
  return(list(pars_top, pars_top_LONG, pars_corr))
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format results for all runs

# Set shared and unit parameters:
if (fit_canada | fit_germany) {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
} else {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
}
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Round 2 results:
res_r2 <- load_and_format_mega_results(filepath = res_dir,
                                       shared_estpars = shared_estpars,
                                       unit_estpars = unit_estpars,
                                       run_name = 'round2')

# ---------------------------------------------------------------------------------------------------------------------

# Get results from round 1

# Read in results:
res_r1 <- read_rds(paste0(res_dir_round1, 'traj_match_round1_byvirseas_TOP.rds'))

# Compile to data frames:
res_r1 <- bind_rows(res_r1)

# Format:
res_r1 <- res_r1 %>%
  mutate(theta_lambda1 = NA,
         theta_lambda2 = NA,
         delta1 = NA,
         d2 = NA,
         alpha = NA,
         phi = NA,
         eta_temp1 = NA,
         eta_temp2 = NA,
         eta_ah1 = NA,
         eta_ah2 = NA,
         b1 = NA,
         b2 = NA,
         phi1 = NA,
         phi2 = NA,
         method = 'round1_noInt') %>%
  select(year:method) %>%
  pivot_longer(-c(year, loglik, method),
               names_to = 'param') %>%
  select(year, param:value, loglik:method)

# ---------------------------------------------------------------------------------------------------------------------

# Visualize results

# Save pairs plots/estimate comparisons as pdf:
if (fit_canada) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_CANADA.pdf'),
      width = 22, height = 12)
} else if (fit_germany) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_GERMANY.pdf'),
      width = 22, height = 12)
} else {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2.pdf'),
      width = 22, height = 12)
}

# Pairs plots:
pairs(res_r2[[3]][, c(shared_estpars, unit_estpars, 'loglik')], pch = 20, main = unique(res_r2[[3]]$method))

# Compare estimates across methods/subtypes:
res <- res_r2[[2]] %>%
  bind_rows(res_r1) %>%
  filter(param %in% c(shared_estpars, unit_estpars))
rm(res_r1)
res$param <- factor(res$param)
res$param <- factor(res$param, levels = c(shared_estpars, unit_estpars))
res$method <- factor(res$method)

p1 <- ggplot(data = res, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~ param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p1)

dev.off()

# Plot correlations between global parameters and Ris/initial conditions:
if (fit_canada) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_corr_CANADA.pdf'),
      width = 18, height = 10)
} else if (fit_germany) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_corr_GERMANY.pdf'),
      width = 18, height = 10)
} else {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_corr.pdf'),
      width = 18, height = 10)
}

for (param in unit_estpars) {
  res_r2[[1]] %>%
    select(all_of(shared_estpars),
           contains(param)) %>%
    plot(pch = 20,
         main = param)
}

# And calculate correlation between estimates of eta_temp and eta_ah:
if (!fit_canada & !fit_germany) {
  par(mfrow = c(2, 1), mar = c(4, 4, 1, 0.5))
  pars_temp <- res_r2[[1]] %>% select(eta_temp1:eta_ah2)
  plot(pars_temp$eta_temp1, pars_temp$eta_ah1, pch = 20)
  plot(pars_temp$eta_temp2, pars_temp$eta_ah2, pch = 20)
  cor.test(pars_temp$eta_temp1, pars_temp$eta_ah1) %>% print()
  cor.test(pars_temp$eta_temp2, pars_temp$eta_ah2) %>% print()
}

dev.off()

# Plot slices:
if (fit_canada) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_slices_CANADA.pdf'),
      width = 20, height = 20)
} else if (fit_germany) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_slices_GERMANY.pdf'),
      width = 20, height = 20)
} else {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_slices.pdf'),
      width = 20, height = 20) 
}

true_estpars <- c(shared_estpars, unit_estpars)
prof_lik <- FALSE

# Set estpars:
estpars <- names(res_r2[[1]])[1:(length(names(res_r2[[1]])) - 1)]

# Set vir1:
if (fit_canada | fit_germany) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

# Read in pomp models:
source('src/functions/setup_global_likelilhood.R')

# Loop through top 5 parameter sets and calculate/plot slices over global params:
par(mfrow = c(10, 6), bty = 'l')

for (j in 1:5) {
  mle <- setNames(object = as.numeric(res_r2[[1]][j, 1:(length(names(res_r2[[1]])) - 1)]),
                  nm = estpars)
  
  if (fit_canada | fit_germany) {
    
    slices <- slice_design(center = mle,
                           rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                           rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20),
                           theta_lambda1 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda1'], to = 1.1 * mle['theta_lambda1'], length.out = 20),
                           theta_lambda2 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda2'], to = 1.1 * mle['theta_lambda2'], length.out = 20),
                           delta1 = seq(from = 0.9 * mle['delta1'], to = 1.1 * mle['delta1'], length.out = 20),
                           d2 = seq(from = 0.9 * mle['d2'], to = 1.1 * mle['d2'], length.out = 20),
                           alpha = seq(from = 0.9 * mle['alpha'], to = 1.1 * mle['alpha'], length.out = 20),
                           phi = seq(from = 0, to = 52, length.out = 20),#seq(from = 0.9 * mle['phi'], to = 1.1 * mle['phi'], length.out = 20),
                           b1 = seq(from = 0.9 * mle['b1'], to = 1.1 * mle['b1'], length.out = 20),
                           b2 = seq(from = 0.9 * mle['b2'], to = 1.1 * mle['b2'], length.out = 20),
                           phi1 = seq(from = 0.9 * mle['phi1'], to = 1.1 * mle['phi1'], length.out = 20),
                           phi2 = seq(from = 0.9 * mle['phi2'], to = 1.1 * mle['phi2'], length.out = 20)) %>%
      mutate(ll = NA)
    
  } else {
    
    slices <- slice_design(center = mle,
                           rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                           rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20),
                           theta_lambda1 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda1'], to = 1.1 * mle['theta_lambda1'], length.out = 20),
                           theta_lambda2 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda2'], to = 1.1 * mle['theta_lambda2'], length.out = 20),
                           delta1 = seq(from = 0.9 * mle['delta1'], to = 1.1 * mle['delta1'], length.out = 20),
                           d2 = seq(from = 0.9 * mle['d2'], to = 1.1 * mle['d2'], length.out = 20),
                           alpha = seq(from = 0.9 * mle['alpha'], to = 1.1 * mle['alpha'], length.out = 20),
                           phi = seq(from = 0, to = 52, length.out = 20),#seq(from = 0.9 * mle['phi'], to = 1.1 * mle['phi'], length.out = 20),
                           eta_temp1 = seq(from = 0.9 * mle['eta_temp1'], to = 1.1 * mle['eta_temp1'], length.out = 20),
                           eta_temp2 = seq(from = 0.9 * mle['eta_temp2'], to = 1.1 * mle['eta_temp2'], length.out = 20),
                           eta_ah1 = seq(from = 0.9 * mle['eta_ah1'], to = 1.1 * mle['eta_ah2'], length.out = 20),
                           eta_ah2 = seq(from = 0.9 * mle['eta_ah2'], to = 1.1 * mle['eta_ah2'], length.out = 20)) %>%
      mutate(ll = NA) 
    
  }
  
  for (k in 1:nrow(slices)) {
    x0 <- slices[k, 1:(length(names(res_r2[[1]])) - 1)]
    expect_true(all(names(x0) == estpars))
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    slices$ll[k] <- -1 * calculate_global_loglik(x0_trans)
  }
  rm(k, x0, x0_trans)
  
  for (par in shared_estpars) {
    slices_cur <- filter(slices, slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = par)
    
  }
  rm(par, slices_cur)
  
}
rm(j, mle, slices)

par(mfrow = c(1, 1), bty = 'l')

for (j in 1:5) {
  mle <- setNames(object = as.numeric(res_r2[[1]][j, 1:(length(names(res_r2[[1]])) - 1)]),
                  nm = estpars)
  
  param_vals <- rbind(c(mle[3], mle[5]),
                      expand.grid(theta_lambda1 = seq(0, 1, by = 0.05), delta1 = 7 / c(15, seq(30, 180, by = 15)))) %>%
    mutate(ll = NA)
  
  for (k in 1:nrow(param_vals)) {
    x0 <- mle
    x0[names(x0) == 'theta_lambda1'] <- param_vals[k, 1]
    x0[names(x0) == 'delta1'] <- param_vals[k, 2]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    param_vals$ll[k] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  p1 <- ggplot(data = as_tibble(param_vals[2:148, ])) + 
    geom_tile(aes(x = theta_lambda1, y = 7 / delta1, fill = ll)) +
    theme_classic() +
    scale_fill_viridis() +
    labs(x = 'theta_lambda1', y = '7 / delta1', fill = 'LL')
  p2 <- ggplot(data = as_tibble(param_vals[2:148, ])) +
    geom_point(aes(x = theta_lambda1, y = ll, col = 7 / delta1)) +
    geom_line(aes(x = theta_lambda1, y = ll, col = 7 / delta1, group = delta1)) +
    geom_point(data = param_vals[1, ], aes(x = theta_lambda1, y = ll, col = 7 / delta1), col = 'red', shape = 8, size = 4) +
    theme_classic() +
    scale_color_viridis() +
    labs(x = 'theta_lambda1', y = 'LL', col = '7 / delta1')
  
  grid.arrange(p1, p2, ncol = 1)
}
rm(j, mle, param_vals)

dev.off()

# Plot simulations:
if (fit_canada) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_simulations_CANADA.pdf'),
      width = 18, height = 10)
} else if (fit_germany) {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_simulations_GERMANY.pdf'),
      width = 18, height = 10)
} else {
  pdf(paste0('results/plots/', date, '_trajectory_matching_round2_simulations.pdf'),
      width = 18, height = 10)
}

# # Read in pomp models:
# source('src/functions/setup_global_likelilhood.R')

# Create list to store plots:
plot_list <- list()

# Simulate/plot each season:
for (j in 1:length(seasons)) {
  
  # Get year:
  yr <- seasons[j]
  
  # Get pomp object:
  resp_mod <- po_list[[j]]
  
  # Get parameter values:
  pars_temp <- res_r2[[1]] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  # Plot top 5 parameter sets:
  for (k in 1:5) {
    coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp[k, ]
    
    sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
    # traj_temp <- trajectory(resp_mod, format = 'data.frame')
    
    sim_temp <- sim_temp %>%
      as_tibble() %>%
      select(time:.id, n_P1:n_P2) %>%
      arrange(.id) %>%
      cbind(t(resp_mod@data))
    names(sim_temp)[5:6] <- c('obs1', 'obs2')
    
    p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
      geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') + 
      geom_point(aes(x = time, y = obs1, group = .id)) + geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
      theme_classic() +
      labs(x = 'Time', y = '# Positive Tests', title = paste('sim', k, sep = '_'))
    plot_list[[j * 5 - 4 + k - 1]] <- p_temp
    
  }
}

# Print plots:
do.call('grid.arrange', c(plot_list, ncol = 5))

dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

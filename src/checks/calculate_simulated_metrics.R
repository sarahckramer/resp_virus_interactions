# ---------------------------------------------------------------------------------------------------------------------
# Code to compare predicted/observed outbreak metrics for influenza and RSV when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Load necessary functions:
source('src/functions/functions_evalutate_res.R')

# Get names of fitted parameters:
shared_estpars_hk <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                       'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
shared_estpars_can <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                        'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars_hk <- c(shared_estpars_hk, unit_estpars)
true_estpars_can <- c(shared_estpars_can, unit_estpars)

# Set parameter values necessary for loading models:
prof_lik <- FALSE

# ---------------------------------------------------------------------------------------------------------------------

# Read in MLEs:
mle_hk <- read_rds('results/MLEs_flu_h1_plus_b.rds')
mle_can <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')

# Simulate several outbreaks for each season (HK):
fit_canada <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk

source('src/functions/setup_global_likelilhood.R')
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$'h1_plus_b_rsv'

res_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Set seed:
  set.seed(1078543)
  
  # Run simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'stochastic', n_sim = 100)
  
  # Calculate metrics (simulated):
  sim_metrics <- calculate_metrics(sim_temp)
  
  # Calculate metrics (observed):
  obs_metrics <- calculate_metrics(hk_dat %>% filter(season == seasons[i]))
  
  # Also calculate correlation coefficients for full outbreaks:
  sim_corr <- sim_temp %>%
    full_join(hk_dat %>% filter(season == seasons[i]),
              by = 'time') %>%
    group_by(.id) %>%
    summarise(corr1 = cor(n_P1.x, n_P1.y, use = 'pairwise.complete.obs'),
              corr2 = cor(n_P2.x, n_P2.y, use = 'pairwise.complete.obs'))
  
  # Determine which simulations within 2wk/25% of observed values:
  sim_metrics <- check_accuracy_metrics(sim_metrics, obs_metrics, seasons[i], pt_acc = 2, pi_acc = 25)
  
  # Store results:
  res_list[[i]] <- sim_metrics
  
}

res_hk <- do.call('rbind', res_list) %>%
  as_tibble() %>%
  mutate(loc = 'hk')

# Simulate several outbreaks for each season (Canada):
fit_canada <- TRUE
vir1 <- 'flu'
true_estpars <- true_estpars_can

source('src/functions/setup_global_likelilhood.R')
can_dat <- read_csv('data/formatted/dat_canada.csv')

res_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Set seed:
  set.seed(1078543)
  
  # Run simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'stochastic', n_sim = 100)
  
  p_temp <- ggplot() +
    geom_line(data = sim_temp, aes(x = time, y = n_P2, group = .id), col = 'gray80', alpha = 0.5) +
    geom_point(data = can_dat %>% filter(season == seasons[i]), aes(x = time, y = n_P2)) +
    theme_classic()
  print(p_temp)
  
  # Calculate metrics (simulated):
  sim_metrics <- calculate_metrics(sim_temp)
  
  # Calculate metrics (observed):
  obs_metrics <- calculate_metrics(can_dat %>% filter(season == seasons[i]))
  
  # Also calculate correlation coefficients for full outbreaks:
  sim_corr <- sim_temp %>%
    full_join(can_dat %>% filter(season == seasons[i]),
              by = 'time') %>%
    group_by(.id) %>%
    summarise(corr1 = cor(n_P1.x, n_P1.y, use = 'pairwise.complete.obs'),
              corr2 = cor(n_P2.x, n_P2.y, use = 'pairwise.complete.obs'))
  
  # Determine which simulations within 2wk/25% of observed values:
  sim_metrics <- check_accuracy_metrics(sim_metrics, obs_metrics, seasons[i], pt_acc = 2, pi_acc = 25)
  
  # Store results:
  res_list[[i]] <- sim_metrics
  
}

res_can <- do.call('rbind', res_list) %>%
  as_tibble() %>%
  mutate(loc = 'canada')

# Combine results from both locations:
res <- bind_rows(res_hk, res_can)
print(res)

# Clean up:
rm(list = ls())

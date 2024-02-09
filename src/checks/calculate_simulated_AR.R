# ---------------------------------------------------------------------------------------------------------------------
# Code to calculate the total attack rate of influenza and RSV each season when parameters are set to best-fit values
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

# Get mean observed cases (HK):
fit_canada <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk

source('src/functions/setup_global_likelilhood.R')

ar_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  res_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'deterministic', obs_only = TRUE, analysis = 'basic')
  
  ar_temp <- res_temp %>%
    group_by(season) %>%
    summarise(H1 = sum(H1),
              H2 = sum(H2),
              obs1 = sum(obs1, na.rm = TRUE),
              obs2 = sum(obs2, na.rm = TRUE))
  
  ar_list[[i]] <- ar_temp
  
}

ar_hk <- bind_rows(ar_list) %>%
  mutate(loc = 'hk')

# Get mean observed cases (Canada):
fit_canada <- TRUE
vir1 <- 'flu'
true_estpars <- true_estpars_can

source('src/functions/setup_global_likelilhood.R')

ar_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  res_temp <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'deterministic', obs_only = TRUE, analysis = 'basic')
  
  ar_temp <- res_temp %>%
    group_by(season) %>%
    summarise(H1 = sum(H1),
              H2 = sum(H2),
              obs1 = sum(obs1, na.rm = TRUE),
              obs2 = sum(obs2, na.rm = TRUE))
  
  ar_list[[i]] <- ar_temp
  
}

ar_can <- bind_rows(ar_list) %>%
  mutate(loc = 'canada')

# Clean up:
rm(dat_pomp, hk_dat, obj_fun_list, po_list, res_temp, ar_temp,
   debug_bool, vir1, vir2, Ri_max1, Ri_max2, d2_max, prof_lik,
   fit_canada, age_structured, yr, yr_index, nrow_check)

# Save results for supplementary figure:
ar_df <- bind_rows(ar_hk, ar_can) %>%
  select(loc, season:H2) %>%
  pivot_longer(H1:H2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'H1', 'Influenza', 'RSV'))
write_rds(ar_df, 'results/simulated_ar.rds')

# Print attack rate ranges:
ar_df %>% filter(loc == 'hk' & virus == 'Influenza') %>% pull(attack_rate) %>% summary() %>% print()
ar_df %>% filter(loc == 'hk' & virus == 'RSV') %>% pull(attack_rate) %>% summary() %>% print()

ar_df %>% filter(loc == 'canada' & virus == 'Influenza') %>% pull(attack_rate) %>% summary() %>% print()
ar_df %>% filter(loc == 'canada' & virus == 'RSV') %>% pull(attack_rate) %>% summary() %>% print()

# Clean up:
rm(list = ls())

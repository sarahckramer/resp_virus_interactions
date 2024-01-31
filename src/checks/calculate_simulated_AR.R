# ---------------------------------------------------------------------------------------------------------------------
# Code to calculate the total attack rate of influenza and RSV each season when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in MLEs:
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')

# Set necessary parameters:
prof_lik <- FALSE
fit_canada <- FALSE

if (fit_canada) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in pomp models:
source('src/functions/setup_global_likelilhood.R')

# Get data:
dat_temp <- hk_dat$h1_plus_b_rsv

# Loop through seasons and get trajectories:
ar_list <- vector('list', length = length(seasons))
for (j in 1:length(seasons)) {
  
  # Get year:
  yr <- seasons[j]
  
  # Get pomp object:
  resp_mod <- po_list[[j]]
  
  # Get parameter values:
  pars_temp <- mle[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[13:19] <- unit_estpars
  
  # Get trajectory at MLE:
  coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp
  traj_temp <- trajectory(resp_mod, format = 'data.frame') %>%
    mutate(season = yr) %>%
    left_join(dat_temp, by = c('time', 'season')) %>%
    rename('i_ILI' = 'GOPC') %>%
    mutate(i_ILI = i_ILI / 1000)
  
  # Calculate attack rates (total):
  ar_tot <- traj_temp %>%
    select(time, H1:H2) %>%
    summarise(H1 = sum(H1),
              H2 = sum(H2)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
  # Calculate mean number of observed cases:
  rho1 <- as.numeric(pars_temp['rho1'])
  rho2 <- as.numeric(pars_temp['rho2'])
  alpha <- as.numeric(pars_temp['alpha'])
  phi <- as.numeric(pars_temp['phi'])
  
  rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H1 / traj_temp$i_ILI
  rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H2 / traj_temp$i_ILI
  
  rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
  rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
  
  expect_equal(nrow(traj_temp), length(rho1_w))
  expect_equal(nrow(traj_temp), length(rho2_w))
  
  traj_temp$rho1_w <- rho1_w
  traj_temp$rho2_w <- rho2_w
  
  ar_obs <- traj_temp %>%
    mutate(obs1 = rho1_w * n_T,
           obs2 = rho2_w * n_T) %>%
    summarise(obs1 = sum(obs1, na.rm = TRUE),
              obs2 = sum(obs2, na.rm = TRUE)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
  # Store ARs in list:
  ar_list[[j]] <- inner_join(ar_tot, ar_obs, by = c('virus1', 'season')) %>%
    select(virus1:season, H1:H2, obs1:obs2)
  
}

# Clean up:
rm(dat_pomp, hk_dat, obj_fun_list, po_list, resp_mod, traj_temp, debug_bool,
   vir1, vir2, Ri_max1, Ri_max2, d2_max, prof_lik, age_structured,
   j, yr, yr_index)

# Save results for supplementary figure:
ar_df <- bind_rows(ar_list) %>%
  select(season:H2) %>%
  pivot_longer(H1:H2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'H1', 'Influenza', 'RSV'))
write_rds(ar_df, 'results/simulated_ar.rds')

# Print attack rate ranges:
ar_df %>% filter(virus == 'Influenza') %>% pull(attack_rate) %>% summary() %>% print()
ar_df %>% filter(virus == 'RSV') %>% pull(attack_rate) %>% summary() %>% print()

# Clean up:
rm(list = ls())

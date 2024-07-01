# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates, peak timings, and outbreak concentration for flu and RSV outbreaks, and compare
# to metrics calculated for simulations when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Load necessary functions:
source('src/functions/functions_evalutate_res.R')

# ---------------------------------------------------------------------------------------------------------------------

# Calculate observed metrics

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv
can_dat <- read_csv('data/formatted/dat_canada.csv')

# Format:
hk_dat <- hk_dat %>%
  select(time, season, pop, n_T:n_P2) %>%
  mutate(n_T1 = n_T,
         n_T2 = n_T) %>%
  select(time:pop, n_T1:n_T2, n_P1:n_P2) %>%
  mutate(loc = 'hk')

can_dat <- can_dat %>%
  select(time, season, pop, n_T1:n_T2, n_P1, n_P2) %>%
  mutate(loc = 'can')

# Loop through all seasons and calculate metrics (HK):
metrics_hk <- vector('list', length = length(unique(hk_dat$season)))

for (i in 1:length(unique(hk_dat$season))) {
  
  # Get relevant season's data:
  dat_temp <- hk_dat %>% filter(season == unique(hk_dat$season)[i])
  
  # Calculate peak timing, peak intensity, and attack rate:
  obs_metrics <- calculate_metrics(dat_temp)# %>%
  # mutate(pt_diff = pt1 - pt2)
  
  # Calculate peak timing based on proportion (rather than number) positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             mutate(prop1 = n_P1 / n_T1, prop2 = n_P2 / n_T2) %>%
                                             summarise(pt_prop1 = which.max(prop1), pt_prop2 = which.max(prop2))) %>%
    mutate(pt_diff = pt_prop1 - pt_prop2)
  
  # Calculate attack rate as total proportion positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             summarise(ar_prop1 = sum(n_P1, na.rm = TRUE) / sum(n_T1, na.rm = TRUE) * 100,
                                                       ar_prop2 = sum(n_P2, na.rm = TRUE) / sum(n_T2, na.rm = TRUE) * 100))
  
  # Calculate duration and concentration (number of weeks containing 75% of reported cases):
  obs_metrics <- obs_metrics %>% bind_cols(calculate_duration_and_concentration(dat_temp))
  
  # Save in list:
  metrics_hk[[i]] <- obs_metrics %>% mutate(season = unique(hk_dat$season)[i])
  
}

metrics_hk <- bind_rows(metrics_hk)
rm(i, hk_dat, dat_temp, obs_metrics)

# Loop through all seasons and calculate metrics (Canada):
metrics_can <- vector('list', length = length(unique(can_dat$season)))

for (i in 1:length(unique(can_dat$season))) {
  
  # Get relevant season's data:
  dat_temp <- can_dat %>% filter(season == unique(can_dat$season)[i])
  
  # Calculate peak timing, peak intensity, and attack rate:
  obs_metrics <- calculate_metrics(dat_temp)
  
  # Calculate peak timing based on proportion (rather than number) positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             mutate(prop1 = n_P1 / n_T1, prop2 = n_P2 / n_T2) %>%
                                             summarise(pt_prop1 = which.max(prop1), pt_prop2 = which.max(prop2))) %>%
    mutate(pt_diff = pt_prop1 - pt_prop2)
  
  # Calculate attack rate as total proportion positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             summarise(ar_prop1 = sum(n_P1, na.rm = TRUE) / sum(n_T1, na.rm = TRUE) * 100,
                                                       ar_prop2 = sum(n_P2, na.rm = TRUE) / sum(n_T2, na.rm = TRUE) * 100))
  
  # Calculate duration and concentration (number of weeks containing 75% of reported cases):
  obs_metrics <- obs_metrics %>% bind_cols(calculate_duration_and_concentration(dat_temp))
  
  # Save in list:
  metrics_can[[i]] <- obs_metrics %>% mutate(season = unique(can_dat$season)[i])
  
}

metrics_can <- bind_rows(metrics_can)
rm(i, can_dat, dat_temp, obs_metrics)

# Print:
metrics_hk %>%
  select(season, contains('1'), pt_diff) %>%
  mutate(pt1 = if_else(season == 's16-17', pt1 + 45 - 53, pt1 + 45 - 52),
         pt_prop1 = if_else(season == 's16-17', pt_prop1 + 45 - 53, pt_prop1 + 45 - 52)) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

metrics_hk %>%
  select(season, contains('2')) %>%
  mutate(pt2 = if_else(season == 's16-17', pt2 + 45 - 53, pt2 + 45 - 52),
         pt_prop2 = if_else(season == 's16-17', pt_prop2 + 45 - 53, pt_prop2 + 45 - 52)) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

metrics_can %>%
  select(season, contains('1'), pt_diff) %>%
  mutate(pt1 = pt1 + 34 - 52,
         pt_prop1 = pt_prop1 + 34 - 52) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

metrics_can %>%
  select(season, contains('2')) %>%
  mutate(pt2 = pt2 + 34 - 52,
         pt_prop2 = pt_prop2 + 34 - 52) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

# Plot "realistic" values:
p1 <- metrics_hk %>% mutate(loc = 'Hong Kong') %>%
  bind_rows(metrics_can %>% mutate(loc = 'Canada')) %>%
  pivot_longer(-c(season:loc), names_to = 'metric', values_to = 'value') %>%
  mutate(vir = if_else(str_detect(metric, '2'), 'RSV', 'Influenza'),
         metric = str_remove(metric, '1'),
         metric = str_remove(metric, '2')) %>%
  filter(metric %in% c('ar_prop', 'dur', 'pt_prop', 'pt_diff')) %>%
  ggplot(aes(x = loc, y = value)) +
  geom_violin(fill = 'gray90') +
  facet_grid(metric ~ vir, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value')
print(p1)

# ---------------------------------------------------------------------------------------------------------------------

# Calculate simulated metrics and compare to observed

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

# Combine results from all locations:
res <- bind_rows(res_hk, res_can)
print(res)

# Clean up
rm(list = ls())

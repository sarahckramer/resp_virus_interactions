# ---------------------------------------------------------------------------------------------------------------------
# Code to compare predicted/observed outbreak metrics for influenza and RSV when parameters are set to best-fit values
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
hk_dat <- hk_dat$h1_plus_b_rsv

# Simulate several outbreaks for each season:
res_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Set seed:
  set.seed(1078543)
  
  # Get year:
  yr <- seasons[i]
  
  # Get pomp object:
  resp_mod <- po_list[[i]]
  
  # Get parameter values:
  pars_temp <- mle[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[13:19] <- unit_estpars
  
  # Simulate at MLE:
  coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame') %>%
    mutate(season = yr) %>%
    select(season, time:.id, n_P1:n_P2) %>%
    as_tibble()
  
  # Calculate peak timing, peak intensity, and attack rate for each run:
  sim_metrics <- sim_temp %>%
    group_by(.id) %>%
    summarise(pt1 = which.max(n_P1),
              pt2 = which.max(n_P2),
              pi1 = max(n_P1, na.rm = TRUE),
              pi2 = max(n_P2, na.rm = TRUE),
              ar1 = sum(n_P1, na.rm = TRUE),
              ar2 = sum(n_P2, na.rm = TRUE))
  
  # Calculate observed peak timing, peak intensity, and attack rate:
  dat_temp <- hk_dat %>%
    filter(season == yr)
  obs_metrics <- dat_temp %>%
    summarise(pt1 = which.max(n_P1),
              pt2 = which.max(n_P2),
              pi1 = max(n_P1, na.rm = TRUE),
              pi2 = max(n_P2, na.rm = TRUE),
              ar1 = sum(n_P1, na.rm = TRUE),
              ar2 = sum(n_P2, na.rm = TRUE))
  
  # Also calculate correlation coefficients for full outbreaks:
  sim_corr <- sim_temp %>%
    full_join(dat_temp, by = 'time') %>%
    group_by(.id) %>%
    summarise(corr1 = cor(n_P1.x, n_P1.y),
              corr2 = cor(n_P2.x, n_P2.y))
  
  # Determine which simulations within 2wk/25% of observed values:
  sim_metrics <- sim_metrics %>%
    mutate(pt1 = (abs(pt1 - obs_metrics$pt1) <= 2),
           pt2 = (abs(pt2 - obs_metrics$pt2) <= 2),
           pi1 = (abs(pi1 - obs_metrics$pi1) + obs_metrics$pi1) <= 1.25 * obs_metrics$pi1,
           pi2 = (abs(pi2 - obs_metrics$pi2) + obs_metrics$pi2) <= 1.25 * obs_metrics$pi2,
           ar1 = (abs(ar1 - obs_metrics$ar1) + obs_metrics$ar1) <= 1.25 * obs_metrics$ar1,
           ar2 = (abs(ar2 - obs_metrics$ar2) + obs_metrics$ar2) <= 1.25 * obs_metrics$ar2) %>%
    select(-.id) %>%
    mutate(season = yr)
  
  # Store results:
  res_list[[i]] <- sim_metrics
  
}

# Clean up:
rm(dat_temp, dat_pomp, hk_dat, obj_fun_list, po_list, resp_mod, pars_temp,
   sim_temp, sim_metrics, obs_metrics, debug_bool, vir1, vir2, Ri_max1, Ri_max2,
   d2_max, prof_lik, nrow_check, sens, age_structured, i, yr, yr_index)

# Determine % of simulations within 2wk/25% of observed values:
res_list <- lapply(res_list, function(ix) {
  c(ix %>% filter(pt1) %>% pull(pt1) %>% length,
    ix %>% filter(pt2) %>% pull(pt1) %>% length,
    ix %>% filter(pi1) %>% pull(pt1) %>% length,
    ix %>% filter(pi2) %>% pull(pt1) %>% length,
    ix %>% filter(ar1) %>% pull(pt1) %>% length,
    ix %>% filter(ar2) %>% pull(pt1) %>% length)
})

res <- do.call('rbind', res_list) %>%
  as_tibble
names(res) <- c('pt1', 'pt2', 'pi1', 'pi2', 'ar1', 'ar2')

print(res)

# Clean up:
rm(list = ls())

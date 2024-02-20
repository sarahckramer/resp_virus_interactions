# ---------------------------------------------------------------------------------------------------------------------
# Calculate population attributable fraction (PAF) of the interaction effect at the MLE
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)

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

# Generate synthetic data at the MLE, with and without interaction

# Read in MLEs:
mle_hk <- read_rds('results/MLEs_flu_h1_plus_b.rds')
mle_can <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')

# Get synthetic data (HK):
fit_canada <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk

source('src/functions/setup_global_likelilhood.R')

traj_list1 = traj_list2 = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  traj_list1[[i]] <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'deterministic', return_obs = FALSE, analysis = 'paf')
  traj_list2[[i]] <- run_sim(po_list[[i]], seasons[i], list(mle_hk, mle_can), shared_estpars_hk, unit_estpars, model_type = 'deterministic', return_obs = FALSE, analysis = 'paf')
  
}

res_hk1 <- bind_rows(traj_list1) %>%
  mutate(loc = 'hk') %>%
  as_tibble()
res_hk2 <- bind_rows(traj_list2) %>%
  mutate(loc = 'hk') %>%
  as_tibble()

# Get synthetic data (Canada):
fit_canada <- TRUE
vir1 <- 'flu'
true_estpars <- true_estpars_can

source('src/functions/setup_global_likelilhood.R')

traj_list1 = traj_list2 = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  traj_list1[[i]] <- run_sim(po_list[[i]], seasons[i], list(mle_can, mle_hk), shared_estpars_can, unit_estpars, model_type = 'deterministic', return_obs = FALSE, analysis = 'paf')
  traj_list2[[i]] <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'deterministic', return_obs = FALSE, analysis = 'paf')
  
}

res_can1 <- bind_rows(traj_list1) %>%
  mutate(loc = 'canada') %>%
  as_tibble()
res_can2 <- bind_rows(traj_list2) %>%
  mutate(loc = 'canada') %>%
  as_tibble()

# Combine results from both locations:
res <- list(bind_rows(res_hk1, res_can1),
            bind_rows(res_hk2, res_can2))

# ---------------------------------------------------------------------------------------------------------------------

# Calculate statistics and plot

# Calculate PAF (or negative of PAF, i.e., % increase in AR expected if interaction were not present):
res_ars <- lapply(res, function(ix) {
  ix %>%
    group_by(loc, season, .id) %>%
    summarise(attack_rate_H1 = sum(H1),
              attack_rate_H2 = sum(H2))
})

for (loc_val in c('hk', 'canada')) {
  print(loc_val)
  
  for (i in 1:length(res_ars)) {
    
    paf_fluOnRSV = paf_RSVonFlu = c()
    res_loc <- res_ars[[i]] %>%
      filter(loc == loc_val)
    
    for (yr in unique(res_loc$season)) {
      
      res_temp <- res_loc %>%
        filter(season == yr)
      
      if (nrow(res_temp) > 0) {
        
        # Impact of flu on RSV:
        ar_pop <- res_temp$attack_rate_H2[res_temp$.id == 1]
        ar_unexposed <- res_temp$attack_rate_H2[res_temp$.id == 3]
        
        # paf_fluOnRSV <- c(paf_fluOnRSV, (ar_pop - ar_unexposed) / ar_pop) # negative
        paf_fluOnRSV <- c(paf_fluOnRSV, (ar_unexposed - ar_pop) / ar_pop)
        
        # Impact of RSV on flu:
        ar_pop <- res_temp$attack_rate_H1[res_temp$.id == 1]
        ar_unexposed <- res_temp$attack_rate_H1[res_temp$.id == 5]
        
        paf_RSVonFlu <- c(paf_RSVonFlu, (ar_unexposed - ar_pop) / ar_pop)
        
      }
      
    }
    
    print(summary(paf_fluOnRSV))
    print(summary(paf_RSVonFlu))
    
  }
  print('')
  
}

# Clean up:
rm(list = ls())

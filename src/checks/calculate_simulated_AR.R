# ---------------------------------------------------------------------------------------------------------------------
# Code to calculate the total attack rate of influenza and RSV each season when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in MLEs:
mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
mle_b <- read_rds('results/MLEs_flu_b.rds')

# Set necessary parameters:
prof_lik <- FALSE

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Loop through results and get trajectories:
ar_list <- vector('list', 2)
for (i in 1:length(ar_list)) {
  
  # Set vir1:
  if (i == 1) {
    vir1 <- 'flu_h1'
    mle <- mle_h1
  } else {
    vir1 <- 'flu_b'
    mle <- mle_b
  }
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get data:
  if (i == 1) {
    dat_temp <- hk_dat$h1_rsv
  } else {
    dat_temp <- hk_dat$b_rsv
  }
  
  # Loop through seasons:
  ar_list_seas <- vector('list', length = length(seasons))
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
    ar_list_seas[[j]] <- inner_join(ar_tot, ar_obs, by = c('virus1', 'season')) %>%
      select(virus1:season, H1:H2, obs1:obs2)
    
  }
  
  # Store data frame of ARs in list:
  ar_list[[i]] <- bind_rows(ar_list_seas)
  
}

# Clean up:
rm(ar_list_seas, dat_pomp, hk_dat, obj_fun_list, po_list, resp_mod, traj_temp,
   debug_bool, vir1, vir2, Ri_max1, Ri_max2, d2_max, prof_lik, lag_val, age_structured,
   i, j, yr, yr_index)

# Plot attack rates by virus/season:
ar_df <- bind_rows(ar_list) %>%
  select(-c(obs1:obs2)) %>%
  pivot_longer(H1:H2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'H1', 'Influenza', 'RSV'),
         virus1 = if_else(virus1 == 'flu_h1', 'H1', 'B'),
         virus1 = factor(virus1),
         virus1 = relevel(virus1, ref = 'H1'),
         attack_rate = attack_rate * 100)
rm(ar_list)

p1 <- ggplot(data = ar_df,
             aes(x = virus, y = attack_rate, group = virus)) +
  geom_violin(fill = 'gray90') +
  theme_classic() +
  facet_wrap(~ virus1) +
  labs(x = 'Virus', y = 'Attack Rate (%)')
print(p1)

# Print attack rate ranges:
ar_df %>% filter(virus1 == 'H1', virus == 'Influenza') %>% pull(attack_rate) %>% summary() %>% print()
ar_df %>% filter(virus1 == 'B', virus == 'Influenza') %>% pull(attack_rate) %>% summary() %>% print()

ar_df %>% filter(virus1 == 'H1', virus == 'RSV') %>% pull(attack_rate) %>% summary() %>% print()
ar_df %>% filter(virus1 == 'B', virus == 'RSV') %>% pull(attack_rate) %>% summary() %>% print()

# Clean up:
rm(list = ls())
dev.off()

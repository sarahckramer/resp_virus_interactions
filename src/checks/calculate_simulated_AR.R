# ---------------------------------------------------------------------------------------------------------------------
# Code to calculate the total attack rate of influenza and RSV each season when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, true_estpars, run_name) {
  
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
  # no_best <- max(no_best, 50)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Return formatted results:
  return(pars_top)
  
}

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1', 'theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in results:
res_h1_round1CI <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci_thetalambda1_thetalambda2/',
                                                true_estpars = true_estpars,
                                                run_name = 'H1_round1CIs')
res_b_round1CI <- load_and_format_mega_results(filepath = 'results/round2_flu_b_round1ci_thetalambda1_thetalambda2/',
                                               true_estpars = true_estpars,
                                               run_name = 'B_round1CIs')

# Compile results:
res_list <- list(res_h1_round1CI, res_b_round1CI)
names(res_list) <- c('flu_h1', 'flu_b')
rm(res_h1_round1CI, res_b_round1CI)

# Set necessary parameters:
prof_lik <- FALSE

# Loop through results and get trajectories:
ar_list <- vector('list', 2)
for (i in 1:length(res_list)) {
  
  # Set vir1:
  if (str_detect(names(res_list)[i], 'h1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Loop through seasons:
  ar_list_seas <- vector('list', length = length(seasons))
  for (j in 1:length(seasons)) {
    
    # Get year:
    yr <- seasons[j]
    
    # Get pomp object:
    resp_mod <- po_list[[j]]
    
    # Get parameter values:
    pars_temp <- res_list[[i]] %>%
      select(all_of(shared_estpars),
             contains(yr))
    names(pars_temp)[6:12] <- unit_estpars
    
    # Get trajectories for top 10 parameter sets:
    ar_list_temp <- vector('list', length = nrow(pars_temp))#min(nrow(pars_temp), 10))
    for (k in 1:nrow(pars_temp)) {#min(nrow(pars_temp), 10)) {
      
      # Set parameter values:
      coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp[k, ]
      
      # Get trajectories:
      traj_temp <- trajectory(resp_mod, format = 'data.frame')
      
      # Calculate attack rates:
      sums_temp <- traj_temp %>%
        select(H1_tot:H2_tot) %>%
        summarise(sum1 = sum(H1_tot),
                  sum2 = sum(H2_tot)) %>%
        mutate(virus1 = vir1,
               season = yr,
               set = k)
      
      # Store in list:
      ar_list_temp[[k]] <- sums_temp
      
    }
    
    # Store data frame of ARs in list:
    ar_list_seas[[j]] <- bind_rows(ar_list_temp)
    
  }
  
  # Store data frame of ARs in list:
  ar_list[[i]] <- bind_rows(ar_list_seas)
  
}

# Clean up:
rm(ar_list_seas, ar_list_temp, dat_pomp, hk_dat, obj_fun_list, po_list, resp_mod, sums_temp, traj_temp,
   debug_bool, vir1, vir2, Ri_max1, Ri_max2, delta_min, early_start_val, prof_lik, i, j, k, yr, yr_index)

# Plot attack rates by virus/season:
ar_df <- bind_rows(ar_list) %>%
  pivot_longer(sum1:sum2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'sum1', 'Flu', 'RSV'),
         virus1 = if_else(virus1 == 'flu_h1', 'H1', 'B'),
         virus1 = factor(virus1),
         virus1 = relevel(virus1, ref = 'H1'),
         attack_rate = attack_rate * 100)
rm(ar_list)

p1 <- ggplot(data = ar_df, aes(x = season, y = attack_rate, fill = virus, group = paste(season, virus))) +
  geom_boxplot() + theme_classic() + facet_wrap(~ virus1) + scale_fill_brewer(palette = 'Set1') +
  labs(x = 'Season', y = 'Attack Rate (%)', fill = 'Virus')
print(p1)

# Print attack rate ranges:
ar_df %>% filter(virus1 == 'H1', virus == 'Flu') %>% group_by(season) %>% summarise(mean = mean(attack_rate), median = median(attack_rate)) %>% print()
ar_df %>% filter(virus1 == 'B', virus == 'Flu') %>% group_by(season) %>% summarise(mean = mean(attack_rate), median = median(attack_rate)) %>% print()

ar_df %>% filter(virus1 == 'H1', virus == 'RSV') %>% group_by(season) %>% summarise(mean = mean(attack_rate), median = median(attack_rate)) %>% print()
ar_df %>% filter(virus1 == 'B', virus == 'RSV') %>% group_by(season) %>% summarise(mean = mean(attack_rate), median = median(attack_rate)) %>% print()

# Clean up:
rm(list = ls())
dev.off()

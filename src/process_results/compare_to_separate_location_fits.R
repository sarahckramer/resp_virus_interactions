# ---------------------------------------------------------------------------------------------------------------------
# Code to compare the log-likelihoods for Hong Kong and Canada when fit simultaneously vs. separately
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set locations of all results:
res_dir_hk <- 'results/round2_fit/round2_3_fluH1_plus_B/'
res_dir_can <- 'results/round2_fit/sens/canada/round2_3_flu/'
res_dir_comb <- 'results/round2_fit/sens/hk_plus_canada/round2_1_/'

# Function to read in and format results:
load_and_format_mega_results <- function(filepath) {
  
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
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  
  # TEMPORARY: If only one top result, take top 25:
  if (no_best == 1) {
    no_best <- 25
  }
  
  print(no_best)
  print(table(pars_df$message))
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # Return formatted results:
  return(pars_top)
  
}

# Load results from individual locations:
res_hk <- load_and_format_mega_results(res_dir_hk)
res_can <- load_and_format_mega_results(res_dir_can)
res_comb <- load_and_format_mega_results(res_dir_comb)

# For combined results, need to calculate log-likelihoods for separate locations
# Load pomp models and relevant functions:
global_estpars <- c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')
shared_estpars <- c('rho1', 'rho2', 'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(global_estpars, shared_estpars, unit_estpars)

sens <- 'sinusoidal_forcing'
source('src/functions/setup_global_likelilhood.R')

# Update function to calculate likelihoods for both locations:
calculate_global_loglik <- function(trans_vals) {
  # Calculates the log-likelihood for each season, and combines to yield a global log-likelihood
  # param trans_vals: Unnamed vector of transformed parameters; fxn only works if this is the only input?
  # returns: The global, negative log-likelihood
  
  # Add names to vector:
  if(is.null(names(trans_vals))) names(trans_vals) <- x0_trans_names
  
  # Split params into global, shared, and unit params:
  shared_hk_in <- as.data.frame(sapply(paste0('^', 'hk_'), grep, x = names(trans_vals)))
  unit_hk_in <- as.data.frame(sapply(paste0('^', seasons_hk, '_hk_'), grep, x = names(trans_vals)))
  names(unit_hk_in) <- seasons_hk
  
  shared_can_in <- as.data.frame(sapply(paste0('^', 'can_'), grep, x = names(trans_vals)))
  unit_can_in <- as.data.frame(sapply(paste0('^', seasons_can, '_can_'), grep, x = names(trans_vals)))
  names(unit_can_in) <- seasons_can
  
  global_params <- trans_vals[-c(unlist(shared_hk_in), unlist(shared_can_in), unique(unlist(unit_hk_in)), unique(unlist(unit_can_in)))]
  
  shared_hk_params <- trans_vals[unlist(shared_hk_in)]
  names(shared_hk_params) <- gsub('hk_', '', names(shared_hk_params))
  shared_hk_params <- c(global_params, shared_hk_params)
  
  shared_can_params <- trans_vals[unlist(shared_can_in)]
  names(shared_can_params) <- gsub('can_', '', names(shared_can_params))
  shared_can_params <- c(global_params, shared_can_params)
  
  unit_hk_params <- list()
  for (i in 1:length(seasons_hk)) {
    unit <- trans_vals[unit_hk_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_hk_params[[i]] <- unit
    }
  }
  
  unit_can_params <- list()
  for (i in 1:length(seasons_can)) {
    unit <- trans_vals[unit_can_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_can_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season (HK):
  units_hk_ll <- rep(NA, length(seasons_hk))
  
  for (i in 1:length(seasons_hk)) {
    if(length(unit_hk_params) > 0) {
      params_temp <- c(shared_hk_params, unit_hk_params[[i]])
    } else {
      params_temp <- c(shared_hk_params)
    }
    
    units_hk_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Get -ll for each season (Canada):
  units_can_ll <- rep(NA, length(seasons_can))
  
  for (i in 1:length(seasons_can)) {
    if (length(unit_can_params) > 0) {
      params_temp <- c(shared_can_params, unit_can_params[[i]])
    } else {
      params_temp <- c(shared_can_params)
    }
    
    units_can_ll[i] <- obj_fun_list[[i + length(seasons_hk)]](params_temp)
  }
  
  # Calculate -ll for each location:
  hk_ll <- sum(units_hk_ll)
  can_ll <- sum(units_can_ll)
  
  # Calculate global -ll:
  glob_ll <- hk_ll + can_ll
  return(c(hk_ll, can_ll, glob_ll))
}

# Calculate log-likelihoods for all best-fit parameter sets:
df_ll <- NULL
for (i in nrow(res_comb)) {
  
  res_comb_red <- res_comb %>%
    select(-loglik)
  
  x0 <- as.numeric(res_comb_red[i, ])
  x0_trans <- transform_params(x0, po_list[[1]], po_list[[7]], seasons_hk, seasons_can, names(res_comb_red), global_estpars, shared_estpars)
  
  ll <- -1 * calculate_global_loglik(x0_trans)
  
  expect_equal(ll[3], res_comb$loglik[i])
  
  df_ll <- rbind(df_ll, ll)
  
}

df_ll <- df_ll %>%
  as_tibble() %>%
  rename('ll_hk' = 'V1',
         'll_can' = 'V2',
         'll_comb' = 'V3') %>%
  pivot_longer(cols = everything()) %>%
  mutate(loc = str_sub(name, 4),
         loglik = value,
         fit = 'combined') %>%
  select(fit, loc:loglik)

# Get log-likelihood values only for individual location runs:
res_hk <- res_hk %>%
  select(loglik) %>%
  mutate(loc = 'hk',
         fit = 'single') %>%
  select(fit, loc, loglik)

res_can <- res_can %>%
  select(loglik) %>%
  mutate(loc = 'can',
         fit = 'single') %>%
  select(fit, loc, loglik)

# Combine all results:
res <- bind_rows(df_ll,
                 res_hk,
                 res_can) %>%
  filter(loc != 'comb')

# Plot:
p1 <- ggplot(data = res, aes(x = fit, y = loglik)) +
  geom_violin(fill = 'gray90') +
  geom_point() +
  facet_wrap(~ loc, scales = 'free_y') +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood')
print(p1)

# Clean up:
rm(list = ls())

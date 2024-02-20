# ---------------------------------------------------------------------------------------------------------------------
# Check whether accounting for stochasticity in simulations leads to substantial variations from deterministic trajectory
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Load necessary functions:
source('src/functions/functions_evalutate_res.R')

# Set directories where final results from round2 fits are stored:
res_dir_hk <- 'results/round2_fit/round2_3_fluH1_plus_B/'
res_dir_can <- 'results/round2_fit/sens/canada/round2_3_flu/'

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

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, true_estpars) {
  
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
  
  # Check whether estimations stopped due to tolerance or time:
  print(table(pars_df$message))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Return formatted results:
  return(pars_top)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Hong Kong

# Read in and format results:
res_hk <- load_and_format_mega_results(filepath = res_dir_hk,
                                       true_estpars = true_estpars_hk)

# Set estpars:
estpars <- names(res_hk)[1:(length(names(res_hk)) - 1)]

# Load pomp models:
fit_canada <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk

source('src/functions/setup_global_likelilhood.R')

# Loop through seasons and get trajectories/simulations:
traj_list = sim_list = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Get deterministic simulations:
  traj_temp <- run_sim(po_list[[i]], seasons[i], res_hk %>% select(-loglik), shared_estpars_hk, unit_estpars, model_type = 'deterministic', return_obs = TRUE) %>%
    select(time:season, obs1:obs2)
  
  # Get stochastic simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], res_hk %>% select(-loglik), shared_estpars_hk, unit_estpars, model_type = 'stochastic', n_sim = 100) %>%
    arrange(.id)
  
  # Store results in lists:
  traj_list[[i]] <- traj_temp
  sim_list[[i]] <- sim_temp
  
}

# Combine and format deterministic trajectories:
dat_traj_hk <- traj_list %>%
  bind_rows() %>%
  pivot_longer(obs1:obs2,
               names_to = 'virus',
               values_to = 'mean') %>%
  mutate(virus = if_else(virus == 'obs1', 'Influenza', 'RSV'),
         loc = 'hk')

# Combine and format stochastic simulations:
dat_sim_hk <- sim_list %>%
  bind_rows() %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'sim') %>%
  mutate(virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'),
         loc = 'hk')

# ---------------------------------------------------------------------------------------------------------------------

# Canada

# Read in and format results:
res_can <- load_and_format_mega_results(filepath = res_dir_can,
                                        true_estpars = true_estpars_can)

# Set estpars:
estpars <- names(res_can)[1:(length(names(res_can)) - 1)]

# Load pomp models:
fit_canada <- TRUE
vir1 <- 'flu'
true_estpars <- true_estpars_can

source('src/functions/setup_global_likelilhood.R')

# Loop through seasons and get trajectories/simulations:
traj_list = sim_list = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Get deterministic simulations:
  traj_temp <- run_sim(po_list[[i]], seasons[i], res_can %>% select(-loglik), shared_estpars_can, unit_estpars, model_type = 'deterministic', return_obs = TRUE) %>%
    select(time:season, obs1:obs2)
  
  # Get stochastic simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], res_can %>% select(-loglik), shared_estpars_can, unit_estpars, model_type = 'stochastic', n_sim = 100) %>%
    arrange(.id)
  
  # Store results in lists:
  traj_list[[i]] <- traj_temp
  sim_list[[i]] <- sim_temp
  
}

# Combine and format deterministic trajectories:
dat_traj_can <- traj_list %>%
  bind_rows() %>%
  pivot_longer(obs1:obs2,
               names_to = 'virus',
               values_to = 'mean') %>%
  mutate(virus = if_else(virus == 'obs1', 'Influenza', 'RSV'),
         loc = 'canada')

# Combine and format stochastic simulations:
dat_sim_can <- sim_list %>%
  bind_rows() %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'sim') %>%
  mutate(virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'),
         loc = 'canada')

# ---------------------------------------------------------------------------------------------------------------------

# Plot

p1 <- ggplot() +
  geom_line(data = dat_sim_hk, aes(x = time, y = sim, group = .id), col = 'gray80', alpha = 0.5) +
  geom_line(data = dat_traj_hk, aes(x = time, y = mean)) +
  facet_grid(virus ~ season, scales = 'free_y') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = 'Time (Weeks)', y = 'Estimated # of Cases')
p2 <- ggplot() +
  geom_line(data = dat_sim_can, aes(x = time, y = sim, group = .id), col = 'gray80', alpha = 0.5) +
  geom_line(data = dat_traj_can, aes(x = time, y = mean)) +
  facet_grid(virus ~ season, scales = 'free_y') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = 'Time (Weeks)', y = 'Estimated # of Cases')

print(p1)
print(p2)

# Clean up:
rm(list = ls())

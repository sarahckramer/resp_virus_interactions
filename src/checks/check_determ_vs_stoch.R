# ---------------------------------------------------------------------------------------------------------------------
# Check whether accounting for stochasticity in simulations leads to substantial variations from deterministic trajectory
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where final results from round2 fits are stored:
res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, shared_estpars, unit_estpars, run_name) {
  
  # Compile estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
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
  expect_equal(df_use, 54)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  # no_best <- max(no_best, 20)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Format for combining with other results:
  pars_top_LONG <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 6),
           param = str_sub(param, 8, str_length(param))) %>%
    pivot_wider(names_from = param, values_from = 'value') %>%
    pivot_longer(-c(loglik, year), names_to = 'param') %>%
    select(year, param:value, loglik) %>%
    mutate(method = run_name)
  
  # Explore correlations:
  pars_corr <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 6),
           param = str_sub(param, 8, str_length(param))) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    # select(-year) %>%
    select(all_of(true_estpars), loglik) %>%
    mutate(method = run_name)
  
  # Return formatted results:
  return(list(pars_top, pars_top_LONG, pars_corr))
  
}

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in and format results:
res_h1_plus_b <- load_and_format_mega_results(filepath = res_dir,
                                              shared_estpars = shared_estpars,
                                              unit_estpars = unit_estpars,
                                              run_name = 'H1_plus_B')

# Set estpars:
estpars <- names(res_h1_plus_b[[1]])[1:(length(names(res_h1_plus_b[[1]])) - 1)]

# Load pomp models:
prof_lik <- FALSE
lag_val <- 0
vir1 <- 'flu_h1_plus_b'
fit_canada <- FALSE

source('src/functions/setup_global_likelilhood.R')

# Loop through seasons and get trajectories/simulations:
traj_list = sim_list = vector('list', length = length(seasons))

for (i in 1:length(seasons)) {
  
  # Get year:
  yr <- seasons[i]
  
  # Get pomp object:
  resp_mod <- po_list[[i]]
  
  # Get parameter values:
  pars_temp <- res_h1_plus_b[[1]][1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  # Set coefficient values:
  coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp
  
  # Get trajectories:
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  
  # Calculate mean expected observed case counts:
  traj_temp <- traj_temp %>%
    as_tibble() %>%
    select(time, H1, H2) %>%
    cbind(resp_mod@covar@table[1, ]) %>%
    cbind(resp_mod@covar@table[2, ])
  names(traj_temp)[4:5] <- c('i_ILI', 'n_T')
  
  rho1 <- as.numeric(pars_temp['rho1'])
  rho2 <- as.numeric(pars_temp['rho2'])
  alpha <- as.numeric(pars_temp['alpha'])
  phi <- as.numeric(pars_temp['phi'])
  
  rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H1 / traj_temp$i_ILI
  rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H2 / traj_temp$i_ILI
  
  rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
  rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
  
  traj_temp$rho1_w <- rho1_w
  traj_temp$rho2_w <- rho2_w
  
  traj_temp <- traj_temp %>%
    mutate(mean1 = rho1_w * n_T,
           mean2 = rho2_w * n_T) %>%
    select(time, mean1:mean2) %>%
    mutate(season = yr) %>%
    as_tibble()
  
  # Get and format simulations:
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    select(time:.id, n_P1:n_P2) %>%
    arrange(.id) %>%
    mutate(season = yr) %>%
    as_tibble()
  
  # Store results in lists:
  traj_list[[i]] <- traj_temp
  sim_list[[i]] <- sim_temp
  
}

# Clean up:
rm(i, traj_temp, sim_temp, yr, resp_mod, pars_temp, rho1, rho2, alpha, phi,
   rho1_w, rho2_w, vir1, vir2, yr_index, Ri_max1, Ri_max2, d2_max, debug_bool,
   age_structured, lag_val, prof_lik, sens, nrow_check)

# Combine and format trajectories:
dat_traj <- traj_list %>%
  bind_rows() %>%
  pivot_longer(mean1:mean2,
               names_to = 'virus',
               values_to = 'mean') %>%
  mutate(virus = if_else(virus == 'mean1', 'Influenza', 'RSV'))

# Combine and format simulations:
dat_sim <- sim_list %>%
  bind_rows() %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'sim') %>%
  mutate(virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'))

# Plot:
p1 <- ggplot() +
  geom_line(data = dat_sim, aes(x = time, y = sim, group = .id), col = 'gray80', alpha = 0.5) +
  geom_line(data = dat_traj, aes(x = time, y = mean)) +
  facet_grid(virus ~ season, scales = 'free_y') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  labs(x = 'Time (Weeks)', y = 'Estimated # of Cases')
print(p1)

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Explore how changing the values of interaction parameters influences model likelihoods/simulations
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(gridExtra)
library(testthat)

# Set whether plots should be stored:
save_plots <- TRUE

# Set model parameters:
prof_lik <- FALSE

# ---------------------------------------------------------------------------------------------------------------------

# Functions

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, true_estpars, run_name) {
  
  # # Compile estpars:
  # true_estpars <- c(shared_estpars, unit_estpars)
  
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
  expect_equal(df_use, 47)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  no_best <- max(no_best, 50)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Return formatted results:
  return(pars_top)
  
}

# Function to generate slices/grids and calculate log-likelihoods:
calculate_slices_and_grid <- function(pars_top, vir1, shared_estpars) {
  
  # Set estpars:
  estpars <- names(pars_top %>%
                     select(-loglik))
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get MLE:
  mle <- setNames(object = as.numeric(pars_top[1, 1:length(estpars)]),
                  nm = estpars)
  
  # Get slices and calculate ll:
  slices <- slice_design(center = mle,
                         theta_lambda1 = seq(from = 0, to = 1.0, by = 0.025),
                         theta_lambda2 = seq(from = 0, to = 1.0, by = 0.025),
                         delta = seq(from = 0.01, to = 1.4, by = 0.025)) %>%
    mutate(ll = NA)
  
  for (ix in 1:nrow(slices)) {
    x0 <- slices[ix, 1:length(estpars)]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    slices$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Get grid and calculate ll:
  theta_lambda1_vals <- seq(from = 0, to = 1.0, by = 0.025)
  delta_vals <- seq(from = 0.01, to = 1.4, by = 0.025)
  
  grid <- expand_grid(theta_lambda1_vals, delta_vals) %>%
    rename('theta_lambda1' = 'theta_lambda1_vals',
           'delta' = 'delta_vals') %>%
    mutate(ll = NA)
  
  x0 <- mle
  for (ix in 1:nrow(grid)) {
    x0[c('theta_lambda1', 'delta')] <- grid[ix, ]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    grid$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Return data frames:
  return(list(slices, grid))
  
}

# Function to simulate for each season and parameter combination:
simulate_grid <- function(pars_top, vir1, shared_estpars, unit_estpars, param1, param2, int_val_ranges, n_sim = 5) {
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get interaction parameter values:
  if (param1 == 'delta') {
    param1_vals <- int_val_ranges[[2]]
  } else {
    param1_vals <- int_val_ranges[[1]]
  }
  
  if (param2 == 'delta') {
    param2_vals <- int_val_ranges[[2]]
  } else {
    param2_vals <- int_val_ranges[[1]]
  }
  
  # Get grid of interaction parameters:
  grid_vals <- expand_grid(param1_vals, param2_vals) %>%
    rename({{param1}} := 'param1_vals',
           {{param2}} := 'param2_vals')
  
  # Get parameters for simulations:
  p_mat_list <- lapply(1:length(seasons), function(ix) {
    pars_temp <- pars_top %>%
      # arrange(desc(loglik)) %>%
      select(all_of(shared_estpars), contains(seasons[ix]))
    names(pars_temp) <- c(shared_estpars, unit_estpars)
    
    coef(po_list[[ix]], names(pars_temp)) <- pars_temp[1, ]
    p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(grid_vals))
    p_mat[param1, ] <- unlist(grid_vals[, param1])
    p_mat[param2, ] <- unlist(grid_vals[, param2])
    
    p_mat
  })
  
  # Simulate for each season and parameter set:
  sim_list_temp <- lapply(1:length(seasons), function(ix) {
    
    simulate(object = po_list[[ix]],
             params = p_mat_list[[ix]],
             nsim = n_sim,
             format = 'data.frame') %>%
      select(time:.id, H1:n_P2) %>%
      pivot_longer(H1:n_P2, names_to = 'output', values_to = 'val') %>%
      mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
             id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2)),
             param1_val = unlist(grid_vals[, param1])[id.parm],
             param2_val = unlist(grid_vals[, param2])[id.parm],
             season = seasons[ix])
    
  })
  
  # Compile simulations:
  sim_df <- bind_rows(sim_list_temp)
  
  # Return simulations:
  return(sim_df)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1', 'theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in results (round 2):
res_h1 <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci/',
                                       true_estpars = true_estpars,
                                       run_name = 'flu_H1') %>%
  mutate(theta_lambda2 = 1.0) %>%
  select(rho1:theta_lambda1, theta_lambda2, `s13-14_Ri1`:loglik)
res_b <- load_and_format_mega_results(filepath = 'results/round2_flu_b_round1ci/',
                                      true_estpars = true_estpars,
                                      run_name = 'flu_B') %>%
  mutate(theta_lambda2 = 1.0) %>%
  select(rho1:theta_lambda1, theta_lambda2, `s13-14_Ri1`:loglik)
res_h1_BOTH <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci_thetalambda1_thetalambda2/',
                                            true_estpars = true_estpars,
                                            run_name = 'flu_H1_BOTH')

res_list <- list(res_h1, res_b, res_h1_BOTH)
names(res_list) <- c('H1', 'B', 'H1_both')
rm(res_h1, res_b, res_h1_BOTH)

# ---------------------------------------------------------------------------------------------------------------------

# Slice/grid likelihoods
# Slices: theta_lambda1, theta_lambda2, delta
# Grid: theta_lambda1, delta

# Loop through results and get slices:
slice_list = grid_list = vector('list', length = length(res_list))
for (i in 1:length(res_list)) {
  
  # Get vir1:
  if (str_detect(names(res_list)[i], 'H1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Calculate ll over slices:
  res_temp <- calculate_slices_and_grid(res_list[[i]], vir1, shared_estpars)
  
  # Store slices:
  slice_list[[i]] <- res_temp[[1]]
  grid_list[[i]] <- res_temp[[2]]
  
}
rm(i, vir1, res_temp, dat_pomp, hk_dat, obj_fun_list, po_list, seasons, yr, yr_index)

# Save plots to file?:
if (save_plots) {
  pdf('results/plots/explore_interaction_param_impact_HK_likelihood.pdf', width = 12, height = 9)
}

# Plot slices:
par(mfrow = c(3, 3))
for (i in 1:length(res_list)) {
  slice_temp <- slice_list[[i]]
  
  for (par in c('theta_lambda1', 'theta_lambda2', 'delta')) {
    slices_cur <- filter(slice_temp, slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = names(res_list)[i])
  }
}
rm(i, slice_temp, par, slices_cur)

# Plot grids:
for (i in 1:length(res_list)) {
  grid_temp <- grid_list[[i]]
  
  p_temp <- ggplot(data = grid_temp, aes(x = theta_lambda1, y = delta, fill = ll)) +
    geom_tile() + theme_classic() +
    scale_fill_viridis() +
    scale_x_continuous(n.breaks = 10, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 20, expand = c(0, 0)) +
    labs(fill = 'LL')
  print(p_temp)
}
rm(i, grid_temp, p_temp)

if (save_plots) {
  dev.off()
}

# ---------------------------------------------------------------------------------------------------------------------

# Simulations at different parameter values

# Set interaction parameter values:
int_vals <- seq(0, 1, by = 0.1)
delta_vals <- c(7/5, 7/14, 7/30, 7/60, 7/200)
int_vals_all <- list(int_vals, delta_vals)
rm(int_vals, delta_vals)

# Loop through results:
sim_list <- vector('list', length = length(res_list))
for (i in 1:length(res_list)) {
  
  # Get vir1:
  if (str_detect(names(res_list)[i], 'H1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Run simulations across interaction parameter values:
  sim_temp <- simulate_grid(res_list[[i]], vir1, shared_estpars, unit_estpars, 'theta_lambda1', 'delta', int_vals_all)
  
  # Store simulations:
  sim_list[[i]] <- sim_temp
  
}
rm(i, int_vals_all, vir1, sim_temp, dat_pomp, hk_dat, obj_fun_list, po_list, seasons, yr, yr_index)

# Generate plots of simulations:
plot_list <- vector('list', length = length(res_list))
seasons_all <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
for (i in 1:length(res_list)) {
  sim_temp <- sim_list[[i]]
  sim_temp <- sim_temp %>%
    mutate(param2_val = 7 / param2_val)
  
  seasons <- res_list[[i]] %>%
    select(contains(seasons_all)) %>%
    pivot_longer(cols = everything(),
                 names_to = 'season') %>%
    mutate(season = str_sub(season, 1, 6)) %>%
    pull(season) %>%
    unique()
  
  plot_list_temp <- vector('list', length = length(seasons))
  
  for (j in 1:length(seasons)) {
    sim_temp_yr <- sim_temp %>%
      filter(season == seasons[j])
    
    p_cases <- ggplot(data = sim_temp_yr %>%
                        filter(output %in% c('H1', 'H2')),
                      aes(x = time, y = val, col = output, group = paste(output, id.sim))) +
      geom_line() + theme_classic() + facet_grid(param1_val ~ param2_val) +
      scale_color_brewer(palette = 'Set1') +
      labs(x = 'Time (Weeks)', y = 'Cases (Total)', col = 'Vir.', title = seasons[j])
    p_obs <- ggplot(data = sim_temp_yr %>%
                      filter(output %in% c('n_P1', 'n_P2')),
                    aes(x = time, y = val, col = output, group = paste(output, id.sim))) +
      geom_line() + theme_classic() + facet_grid(param1_val ~ param2_val) +
      scale_color_brewer(palette = 'Set1') +
      labs(x = 'Time (Weeks)', y = 'Postive Tests', col = 'Vir.', title = seasons[j])
    
    # p_obs_alt1 <- ggplot(data = sim_temp_yr %>%
    #                        filter(output %in% c('n_P1', 'n_P2')),
    #                      aes(x = time, y = val, col = param1_val, lty = output, group = paste(output, id.sim, param1_val))) +
    #   geom_line() + theme_classic() + facet_wrap(~ param2_val) +
    #   scale_color_viridis() +
    #   labs(x = 'Time (Weeks)', y = 'Postive Tests', col = 'theta_lambda1', title = seasons[j])
    # p_obs_alt2 <- ggplot(data = sim_temp_yr %>%
    #                        filter(output %in% c('n_P1', 'n_P2')),
    #                      aes(x = time, y = val, col = param2_val, lty = output, group = paste(output, id.sim, param2_val))) +
    #   geom_line() + theme_classic() + facet_wrap(~ param1_val) +
    #   scale_color_viridis() +
    #   labs(x = 'Time (Weeks)', y = 'Postive Tests', col = 'delta', title = seasons[j])
    # print(p_obs_alt1)
    # print(p_obs_alt2)
    
    plot_list_temp[[j]] <- list(p_cases, p_obs)
  }
  rm(j, sim_temp_yr, p_cases, p_obs)
  
  plot_list[[i]] <- plot_list_temp
}
rm(i, sim_temp, seasons_all, seasons, plot_list_temp)

# Save plots to file:
if (save_plots) {
  pdf('results/plots/explore_interaction_param_impact_HK_simulations.pdf',
      width = 15, height = 10)
  print(plot_list)
  dev.off()
}

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

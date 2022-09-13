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

# Function to generate slices/grids and calculate log-likelihoods:
calculate_slices_and_grid <- function(mle) {
  
  # Set estpars:
  estpars <- names(mle)
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get slices and calculate ll:
  slices <- slice_design(center = unlist(mle),
                         theta_lambda1 = seq(from = 0, to = 1.0, by = 0.025),
                         theta_lambda2 = seq(from = 0, to = 1.0, by = 0.025),
                         delta1 = seq(from = 0.01, to = 1.4, by = 0.025),
                         d2 = seq(from = 0.5, to = 5.0, by = 0.5)) %>%
    mutate(ll = NA)
  
  for (ix in 1:nrow(slices)) {
    x0 <- slices[ix, 1:length(estpars)]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    slices$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Get grid and calculate ll:
  theta_lambda1_vals <- seq(from = 0, to = 1.0, by = 0.025)
  delta1_vals <- seq(from = 0.01, to = 1.4, by = 0.025)
  
  grid <- expand_grid(theta_lambda1_vals, delta1_vals) %>%
    rename('theta_lambda1' = 'theta_lambda1_vals',
           'delta1' = 'delta1_vals') %>%
    mutate(ll = NA)
  
  x0 <- mle
  for (ix in 1:nrow(grid)) {
    x0[c('theta_lambda1', 'delta1')] <- grid[ix, c('theta_lambda1', 'delta1')]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    grid$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  theta_lambda2_vals <- seq(from = 0, to = 1.0, by = 0.025)
  d2_vals <- seq(from = 0.5, to = 5.0, by = 0.5)
  
  grid2 <- expand_grid(theta_lambda2_vals, d2_vals) %>%
    rename('theta_lambda2' = 'theta_lambda2_vals',
           'd2' = 'd2_vals') %>%
    mutate(ll = NA)
  
  x0 <- mle
  for (ix in 1:nrow(grid2)) {
    x0[c('theta_lambda2', 'd2')] <- grid2[ix, c('theta_lambda2', 'd2')]
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    grid2$ll[ix] <- -1 * calculate_global_loglik(x0_trans)
  }
  
  # Return data frames:
  return(list(slices, grid, grid2))
  
}

# Function to simulate for each season and parameter combination:
simulate_grid <- function(mle, param1, int_val_ranges, n_sim = 5) {
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Get interaction parameter values:
  if (param1 == 'theta_lambda1') {
    
    param2 <- 'delta1'
    
    param1_vals <- int_val_ranges[[1]]
    param2_vals <- int_val_ranges[[2]]
    
  } else if (param1 == 'theta_lambda2') {
    
    param2 <- 'd2'
    
    param1_vals <- int_val_ranges[[1]]
    param2_vals <- int_val_ranges[[3]]
    
  }
  
  # Get grid of interaction parameters:
  grid_vals <- expand_grid(param1_vals, param2_vals) %>%
    rename({{param1}} := 'param1_vals',
           {{param2}} := 'param2_vals')
  
  # Get parameters for simulations:
  p_mat_list <- lapply(1:length(seasons), function(ix) {
    
    pars_temp <- mle %>%
      select(all_of(shared_estpars), contains(seasons[ix]))
    names(pars_temp) <- c(shared_estpars, unit_estpars)
    
    coef(po_list[[ix]], names(pars_temp)) <- pars_temp
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
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Read in MLEs:
mle_h1 <- read_rds('results/MLEs_flu_h1.rds')[1, ]
mle_b <- read_rds('results/MLEs_flu_b.rds')[1, ]

# ---------------------------------------------------------------------------------------------------------------------

# Slice/grid likelihoods

# Get slices/grids:
vir1 <- 'flu_h1'
slice_grid_ll_h1 <- calculate_slices_and_grid(mle_h1)

vir1 <- 'flu_b'
slice_grid_ll_b <- calculate_slices_and_grid(mle_b)

res_list <- list(slice_grid_ll_h1, slice_grid_ll_b)
rm(slice_grid_ll_h1, slice_grid_ll_b)

# Save plots to file?:
if (save_plots) {
  pdf('results/plots/explore_interaction_param_impact_HK_likelihood.pdf', width = 12, height = 9)
}

# Plot slices:
par(mfrow = c(2, 4))
for (i in 1:length(res_list)) {
  
  slice_temp <- res_list[[i]][[1]]
  
  for (par in c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')) {
    slices_cur <- filter(slice_temp, slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = names(res_list)[i])
  }
  
}
rm(i, slice_temp, par, slices_cur)

# Plot grids:
for (i in 1:length(res_list)) {
  
  grid_temp1 <- res_list[[i]][[2]]
  
  p_temp1 <- ggplot(data = grid_temp1, aes(x = theta_lambda1, y = delta1, fill = ll)) +
    geom_tile() + theme_classic() +
    scale_fill_viridis() +
    scale_x_continuous(n.breaks = 10, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 20, expand = c(0, 0)) +
    labs(fill = 'LL')
  
  grid_temp2 <- res_list[[i]][[3]]
  p_temp2 <- ggplot(data = grid_temp2, aes(x = theta_lambda2, y = d2, fill = ll)) +
    geom_tile() + theme_classic() +
    scale_fill_viridis() +
    scale_x_continuous(n.breaks = 10, expand = c(0, 0)) +
    scale_y_continuous(n.breaks = 20, expand = c(0, 0)) +
    labs(fill = 'LL')
  
  grid.arrange(p_temp1, p_temp2, nrow = 1)
  
}
rm(i, grid_temp1, grid_temp2, p_temp1, p_temp2)

if (save_plots) {
  dev.off()
}

# ---------------------------------------------------------------------------------------------------------------------

# Simulations at different parameter values

# Set interaction parameter values:
int_vals <- seq(0, 1, by = 0.1)
delta1_vals <- c(7/5, 7/14, 7/30, 7/60, 7/200)
d2_vals <- c(0.5, 1.0, 1.5, 2.0, 5.0)
int_vals_all <- list(int_vals, delta1_vals, d2_vals)
rm(int_vals, delta1_vals, d2_vals)

# Run simulations across interaction parameter values:
vir1 <- 'flu_h1'
sim_grid_h1_1 <- simulate_grid(mle_h1, 'theta_lambda1', int_vals_all)
sim_grid_h1_2 <- simulate_grid(mle_h1, 'theta_lambda2', int_vals_all)

vir1 <- 'flu_b'
sim_grid_b_1 <- simulate_grid(mle_b, 'theta_lambda1', int_vals_all)
sim_grid_b_2 <- simulate_grid(mle_b, 'theta_lambda2', int_vals_all)

res_list <- list(sim_grid_h1_1, sim_grid_h1_2, sim_grid_b_1, sim_grid_b_2)
rm(sim_grid_h1_1, sim_grid_h1_2, sim_grid_b_1, sim_grid_b_2)

# Generate plots of simulations:
plot_list <- vector('list', length = length(res_list))
for (i in 1:length(res_list)) {
  
  sim_temp <- res_list[[i]]
  if (i %in% c(1, 3)) {
    sim_temp <- sim_temp %>%
      mutate(param2_val = 7 / param2_val)
  }
  
  seasons <- unique(sim_temp$season)
  
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
    
    plot_list_temp[[j]] <- list(p_cases, p_obs)
    
  }
  
  rm(j, sim_temp_yr, p_cases, p_obs)
  
  plot_list[[i]] <- plot_list_temp
  
}
rm(i, sim_temp, seasons, plot_list_temp)

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

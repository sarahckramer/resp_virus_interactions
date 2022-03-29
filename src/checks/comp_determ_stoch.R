# ---------------------------------------------------------------------------------------------------------------------
# Code to compare deterministic and stochastic model runs at the MLE, and determine timestep for stochastic model
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(gridExtra)
library(testthat)

# Set necessary parameters for pomp model:
prof_lik <- FALSE
lag_val <- 0

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# List filepaths to read results from:
filepath_list <- list('results/round2_2_fluH1_FULL/', 'results/round2_2_fluB_FULL/')
names(filepath_list) <- c('flu_H1', 'flu_B')

# Loop through filepaths:
for (i in 1:length(filepath_list)) {
  
  # Get list of results files:
  res_files <- list.files(path = filepath_list[[i]], full.names = TRUE)
  
  # Read in results:
  res_full <- list()
  for (j in seq_along(res_files)) {
    res_full[[j]] <- read_rds(res_files[[j]])
  }
  
  # Get parameter estimates and log-likelihoods:
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only MLE:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  mle <- pars_df[1, ]
  
  # Set vir1:
  if (str_detect(names(filepath_list)[i], 'H1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  print(vir1)
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Create list to store plots:
  plot_list <- list()
  
  # Simulate/plot each season:
  for (k in 1:length(seasons)) {
    
    # Get year:
    yr <- seasons[k]
    
    # Get pomp object:
    resp_mod <- po_list[[k]]
    
    # Get parameter values:
    mle_yr <- mle %>%
      select(all_of(shared_estpars),
             contains(yr))
    names(mle_yr)[(length(names(mle_yr)) - 6):length(names(mle_yr))] <- unit_estpars
    
    # Set parameters for simulation:
    coef(resp_mod, true_estpars) <- mle_yr
    
    # Run deterministic and stochastic model:
    traj_temp <- trajectory(resp_mod, format = 'data.frame')
    sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
    
    traj_temp <- traj_temp %>%
      as_tibble() %>%
      select(time:.id, H1:H2) %>%
      mutate(type = 'deterministic',
             .id = as.character(.id))
    sim_temp <- sim_temp %>%
      as_tibble() %>%
      select(time:.id, H1:H2) %>%
      mutate(type = 'stochastic',
             .id = as.character(.id))
    
    sim_temp <- bind_rows(traj_temp, sim_temp)
    rm(traj_temp)
    
    p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = H1, group = paste(type, .id), col = type, size = type)) +
      geom_line(aes(x = time, y = H2, group = paste(type, .id), col = type, size = type), lty = 2) +
      theme_classic() + scale_color_brewer(palette = 'Set1') +
      scale_size_discrete(range = c(1.1, 0.7)) +
      labs(x = 'Time', y = '# of Cases', col = '', title = paste(vir1, yr, sep = '_'))
    plot_list[[k]] <- p_temp
    
  }
  
  # Print plots:
  do.call('grid.arrange', c(plot_list, ncol = 2))
  
}

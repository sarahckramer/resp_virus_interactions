# ---------------------------------------------------------------------------------------------------------------------
# Format and plot results from "round 2" of trajectory matching (using mega-likelihood)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Get/format date (for saving results):
date <- format(Sys.Date(), '%d%m%y')

# ---------------------------------------------------------------------------------------------------------------------

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
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
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

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format results for all runs

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# H1, using round 1 CIs, all estpars:
res_h1_round1CI <- load_and_format_mega_results(filepath = 'results/round2_2_fluH1_FULL/',
                                                shared_estpars = shared_estpars,
                                                unit_estpars = unit_estpars,
                                                run_name = 'H1_round1CIs_FULL')

# B, using round 1 CIs, all estpars:
res_b_round1CI <- load_and_format_mega_results(filepath = 'results/round2_2_fluB_FULL/',
                                               shared_estpars = shared_estpars,
                                               unit_estpars = unit_estpars,
                                               run_name = 'B_round1CIs_FULL')

# Extract results:
res_LIST <- list(res_h1_round1CI, res_b_round1CI)
pars_top_LIST = pars_top_long_LIST = pars_corr_LIST = vector('list', length = length(res_LIST))
for (i in 1:length(res_LIST)) {
  pars_top_LIST[[i]] <- res_LIST[[i]][[1]]
  pars_top_long_LIST[[i]] <- res_LIST[[i]][[2]]
  pars_corr_LIST[[i]] <- res_LIST[[i]][[3]]
}
rm(i)

names(pars_top_LIST) <- c('flu_h1_round1CI_FULL', 'flu_b_round1CI_FULL')
names(pars_top_long_LIST) <- c('flu_h1_round1CI_FULL', 'flu_b_round1CI_FULL')
names(pars_corr_LIST) <- c('flu_h1_round1CI_FULL', 'flu_b_round1CI_FULL')

# Clean up:
try(rm(res_LIST, res_h1_round1CI, res_b_round1CI))

# ---------------------------------------------------------------------------------------------------------------------

# Get results from round 1

# Read in results:
res_r1 <- read_rds('results/round1_fitsharedFALSE/traj_match_round1_byvirseas_TOP.rds')

# Compile to data frames:
res_r1 <- bind_rows(res_r1)

# Format:
res_r1 <- res_r1 %>%
  mutate(theta_lambda1 = NA,
         theta_lambda2 = NA,
         delta1 = NA,
         d2 = NA,
         alpha = NA,
         phi = NA,
         eta_temp1 = NA,
         eta_temp2 = NA,
         eta_ah1 = NA,
         eta_ah2 = NA,
         method = if_else(virus1 == 'flu_b',
                          'B_round1_noInt',
                          'H1_round1_noInt')) %>%
  select(year:method) %>%
  pivot_longer(-c(year, loglik, method),
               names_to = 'param') %>%
  select(year, param:value, loglik:method)

# ---------------------------------------------------------------------------------------------------------------------

# Visualize results

# Save pairs plots/estimate comparisons as pdf:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2.pdf'),
    width = 22, height = 12)

# Pairs plots:
lapply(pars_corr_LIST, function(ix) {
  pairs(ix[, c(shared_estpars, unit_estpars, 'loglik')], pch = 20, main = unique(ix$method))
})

# Compare estimates across methods/subtypes:
res <- bind_rows(pars_top_long_LIST) %>%
  bind_rows(res_r1) %>%
  mutate(vir1 = str_sub(method, 1, 2))
rm(res_r1)
res$param <- factor(res$param)
res$param <- factor(res$param, levels = levels(res$param)[c(14:15, 18:19, 3, 2, 1, 10, 6:7, 4:5, 16:17, 8:9, 11, 13, 12)])

res_h1 <- res %>% filter(vir1 == 'H1')
res_b <- res %>% filter(vir1 == 'B_') %>% mutate(vir1 = 'B')

res_h1$method <- factor(res_h1$method)
res_h1$method <- factor(res_h1$method, levels = levels(res_h1$method)[2:1])

res_b$method <- factor(res_b$method)
res_b$method <- factor(res_b$method, levels = levels(res_b$method)[2:1])

p1 <- ggplot(data = res_h1, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~ param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
p2 <- ggplot(data = res_b, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~ param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')

print(p1)
print(p2)

dev.off()

# Plot correlations between global parameters and Ris/initial conditions:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2_corr.pdf'),
    width = 18, height = 10)

for (i in 1:length(pars_top_LIST)) {
  for (param in unit_estpars) {
    pars_top_LIST[[i]] %>%
      select(all_of(shared_estpars),
             contains(param)) %>%
      plot(pch = 20,
           main = paste(unique(pars_top_long_LIST[[i]]$method),
                        param,
                        sep = '_'))
  }
}
rm(i)

# And calculate correlation between estimates of eta_temp and eta_ah:
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.5))
for (i in 1:length(pars_top_LIST)) {
  pars_temp <- pars_top_LIST[[i]] %>% select(eta_temp1:eta_ah2)
  plot(pars_temp$eta_temp1, pars_temp$eta_ah1, pch = 20)
  plot(pars_temp$eta_temp2, pars_temp$eta_ah2, pch = 20)
  cor.test(pars_temp$eta_temp1, pars_temp$eta_ah1) %>% print()
  cor.test(pars_temp$eta_temp2, pars_temp$eta_ah2) %>% print()
}
dev.off()

# Plot slices:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2_slices.pdf'),
    width = 20, height = 20)

true_estpars <- c(shared_estpars, unit_estpars)
prof_lik <- FALSE
lag_val <- 0

for (i in 1:length(pars_top_LIST)) {
  
  # Set estpars:
  estpars <- names(pars_top_LIST[[i]])[1:(length(shared_estpars) + length(unit_estpars) * 5)]
  
  # Set vir1:
  if (str_detect(names(pars_top_LIST)[i], 'h1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Loop through top 5 parameter sets and calculate/plot slices over global params:
  par(mfrow = c(10, 6), bty = 'l')
  
  for (j in 1:5) {
    mle <- setNames(object = as.numeric(pars_top_LIST[[i]][j, 1:(length(shared_estpars) + length(unit_estpars) * 5)]),
                    nm = estpars)
    # # slices <- slice_design(center = mle,
    # #                        rho1 = seq(from = 0, to = 1.0, by = 0.05),
    # #                        rho2 = seq(from = 0, to = 1.0, by = 0.05),
    # #                        theta_lambda1 = seq(from = 0, to = 1.0, by = 0.05),
    # #                        delta = 7 / seq(from = 1, to = 365, by = 30)) %>%
    # #   mutate(ll = NA)
    # slices <- slice_design(center = mle,
    #                        rho1 = seq(from = 0.01, to = 0.2, by = 0.01),
    #                        rho2 = seq(from = 0.01, to = 0.2, by = 0.01),
    #                        theta_lambda1 = seq(from = 0, to = 0.2, by = 0.01),
    #                        theta_lambda2 = seq(from = 0, to = 0.2, by = 0.01),
    #                        delta = seq(from = 0.01, to = 0.3, by = 0.01)) %>%
    #   mutate(ll = NA)
    slices <- slice_design(center = mle,
                           rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                           rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20),
                           theta_lambda1 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda1'], to = 1.1 * mle['theta_lambda1'], length.out = 20),
                           theta_lambda2 = seq(from = 0, to = 1, by = 0.01), #(from = 0.9 * mle['theta_lambda2'], to = 1.1 * mle['theta_lambda2'], length.out = 20),
                           delta1 = seq(from = 0.9 * mle['delta1'], to = 1.1 * mle['delta1'], length.out = 20),
                           d2 = seq(from = 0.9 * mle['d2'], to = 1.1 * mle['d2'], length.out = 20),
                           alpha = seq(from = 0.9 * mle['alpha'], to = 1.1 * mle['alpha'], length.out = 20),
                           phi = seq(from = 0, to = 52, length.out = 20),#seq(from = 0.9 * mle['phi'], to = 1.1 * mle['phi'], length.out = 20),
                           eta_temp1 = seq(from = 0.9 * mle['eta_temp1'], to = 1.1 * mle['eta_temp1'], length.out = 20),
                           eta_temp2 = seq(from = 0.9 * mle['eta_temp2'], to = 1.1 * mle['eta_temp2'], length.out = 20),
                           eta_ah1 = seq(from = 0.9 * mle['eta_ah1'], to = 1.1 * mle['eta_ah2'], length.out = 20),
                           eta_ah2 = seq(from = 0.9 * mle['eta_ah2'], to = 1.1 * mle['eta_ah2'], length.out = 20)) %>%
      mutate(ll = NA)
    
    for (k in 1:nrow(slices)) {
      x0 <- slices[k, 1:(length(shared_estpars) + length(unit_estpars) * 5)]
      expect_true(all(names(x0) == estpars))
      x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
      slices$ll[k] <- -1 * calculate_global_loglik(x0_trans)
    }
    rm(k, x0, x0_trans)
    
    for (par in shared_estpars) {
      slices_cur <- filter(slices, slice == par)
      plot(slices_cur[[par]], slices_cur$ll, type = 'l',
           xlab = par, ylab = 'Log-Likelihood',
           main = par)
      
    }
    rm(par, slices_cur)
    
  }
  rm(j, mle, slices)
  
}
rm(i)

dev.off()

# Plot simulations:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2_simulations.pdf'),
    width = 18, height = 10)

for (i in 1:length(pars_top_LIST)) {
  
  # Set vir1:
  if (str_detect(names(pars_top_LIST)[i], 'h1')) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Create list to store plots:
  plot_list <- list()
  
  # Simulate/plot each season:
  for (j in 1:length(seasons)) {
    
    # Get year:
    yr <- seasons[j]
    
    # Get pomp object:
    resp_mod <- po_list[[j]]
    
    # Get parameter values:
    pars_temp <- pars_top_LIST[[i]] %>%
      select(all_of(shared_estpars),
             contains(yr))
    names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
    
    # Plot top 5 parameter sets:
    for (k in 1:5) {
      coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp[k, ]
      
      sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
      # traj_temp <- trajectory(resp_mod, format = 'data.frame')
      
      sim_temp <- sim_temp %>%
        as_tibble() %>%
        select(time:.id, n_P1:n_P2) %>%
        arrange(.id) %>%
        cbind(t(resp_mod@data))
      names(sim_temp)[5:6] <- c('obs1', 'obs2')
      
      p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
        geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') + 
        geom_point(aes(x = time, y = obs1, group = .id)) + geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
        theme_classic() +
        labs(x = 'Time', y = '# Positive Tests', title = paste(unique(pars_top_long_LIST[[i]]$method), k, sep = '_'))
      plot_list[[j * 5 - 4 + k - 1]] <- p_temp
      
    }
  }
  
  # Print plots:
  do.call('grid.arrange', c(plot_list, ncol = 5))
}

dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

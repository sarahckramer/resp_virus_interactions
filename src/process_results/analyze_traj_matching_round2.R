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
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
  no_best <- max(no_best, 50)
  
  pars_top <- pars_df[1:no_best, ]
  
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
shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# H1, using round 1 CIs:
res_h1_round1CI <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_round1ci/',
                                                shared_estpars = shared_estpars,
                                                unit_estpars = unit_estpars,
                                                run_name = 'H1_round1CIs')

# H1, broad:
res_h1_broad <- load_and_format_mega_results(filepath = 'results/round2_flu_h1_broad/',
                                             shared_estpars = shared_estpars,
                                             unit_estpars = unit_estpars,
                                             run_name = 'H1_broad')

# B, using round 1 CIs:
res_b_round1CI <- load_and_format_mega_results(filepath = 'results/round2_flu_b_round1ci/',
                                               shared_estpars = shared_estpars,
                                               unit_estpars = unit_estpars,
                                               run_name = 'B_round1CIs')

# B, broad:
res_b_broad <- load_and_format_mega_results(filepath = 'results/round2_flu_b_broad/',
                                            shared_estpars = shared_estpars,
                                            unit_estpars = unit_estpars,
                                            run_name = 'B_broad')

# Extract results:
res_LIST <- list(res_h1_round1CI, res_h1_broad, res_b_round1CI, res_b_broad)
pars_top_LIST = pars_top_long_LIST = pars_corr_LIST = vector('list', length = length(res_LIST))
for (i in 1:length(res_LIST)) {
  pars_top_LIST[[i]] <- res_LIST[[i]][[1]]
  pars_top_long_LIST[[i]] <- res_LIST[[i]][[2]]
  pars_corr_LIST[[i]] <- res_LIST[[i]][[3]]
}
rm(i)

# Clean up:
rm(res_LIST, res_h1_round1CI, res_h1_broad, res_b_round1CI, res_b_broad)

# ---------------------------------------------------------------------------------------------------------------------

# Get results from round 1

# Read in results:
res_r1_noInt <- read_rds('results/traj_match_round1_byvirseas_TOP_noInt.rds')
res_r1_int <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

# Compile to data frames:
res_r1_noInt <- bind_rows(res_r1_noInt)
res_r1_int <- bind_rows(res_r1_int)

# Format:
res_r1_noInt <- res_r1_noInt %>%
  mutate(delta = NA,
         theta_lambda1 = NA,
         method = if_else(virus1 == 'flu_b',
                          'B_round1_noInt',
                          'H1_round1_noInt')) %>%
  select(year:method) %>%
  pivot_longer(-c(year, loglik, method),
               names_to = 'param') %>%
  select(year, param:value, loglik:method)

res_r1_int <- res_r1_int %>%
  mutate(method = if_else(virus1 == 'flu_b',
                          'B_round1_wInt',
                          'H1_round1_wInt')) %>%
  select(year:method) %>%
  pivot_longer(-c(year, loglik, method),
               names_to = 'param') %>%
  select(year, param:value, loglik:method)

# Join:
res_r1 <- bind_rows(res_r1_noInt, res_r1_int)
rm(res_r1_noInt, res_r1_int)

# ---------------------------------------------------------------------------------------------------------------------

# Visualize results

# Save pairs plots/estimate comparisons as pdf:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2.pdf'),
    width = 18, height = 10)

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
res$param <- factor(res$param, levels = levels(res$param)[c(9:10, 2:4, 6, 5, 7:8, 1, 11)])

res_h1 <- res %>% filter(vir1 == 'H1')
res_b <- res %>% filter(vir1 == 'B_') %>% mutate(vir1 = 'B')

res_h1$method <- factor(res_h1$method)
res_h1$method <- factor(res_h1$method, levels = levels(res_h1$method)[c(4, 1, 3, 2)])

res_b$method <- factor(res_b$method)
res_b$method <- factor(res_b$method, levels = levels(res_b$method)[c(4, 1, 3, 2)])

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

dev.off()

# Plot slices:
pdf(paste0('results/plots/', date, '_trajectory_matching_round2_slices.pdf'),
    width = 18, height = 10)

true_estpars <- c(shared_estpars, unit_estpars)
prof_lik <- FALSE

for (i in 1:length(pars_top_LIST)) {
  
  # Set estpars:
  estpars <- names(pars_top_LIST[[i]])[1:39]
  
  # Set vir1:
  if (i %in% 1:2) {
    vir1 <- 'flu_h1'
  } else {
    vir1 <- 'flu_b'
  }
  
  # Read in pomp models:
  source('src/functions/setup_global_likelilhood.R')
  
  # Loop through top 5 parameter sets and calculate/plot slices over global params:
  par(mfrow = c(5, 4), bty = 'l')
  
  for (j in 1:5) {
    mle <- setNames(object = as.numeric(pars_top_LIST[[i]][j, 1:39]),
                    nm = estpars)
    # slices <- slice_design(center = mle,
    #                        rho1 = seq(from = 0, to = 1.0, by = 0.05),
    #                        rho2 = seq(from = 0, to = 1.0, by = 0.05),
    #                        theta_lambda1 = seq(from = 0, to = 1.0, by = 0.05),
    #                        delta = 7 / seq(from = 1, to = 365, by = 30)) %>%
    #   mutate(ll = NA)
    slices <- slice_design(center = mle,
                           rho1 = seq(from = 0.01, to = 0.2, by = 0.01),
                           rho2 = seq(from = 0.01, to = 0.2, by = 0.01),
                           theta_lambda1 = seq(from = 0, to = 0.2, by = 0.01),
                           delta = seq(from = 0.01, to = 0.3, by = 0.01)) %>%
      mutate(ll = NA)
    
    for (k in 1:nrow(slices)) {
      x0 <- slices[k, 1:39]
      x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
      slices$ll[k] <- -1 * calculate_global_loglik(x0_trans)
    }
    rm(k, x0, x0_trans)
    
    for (par in shared_estpars) {
      slices_cur <- filter(slices, slice == par)
      
      if (par == 'delta') {
        plot(slices_cur[[par]], slices_cur$ll, type = 'l',
             xlab = par, ylab = 'Log-Likelihood',
             main = par)
      } else {
        plot(slices_cur[[par]], slices_cur$ll, type = 'l',
             xlab = par, ylab = 'Log-Likelihood',
             main = par)
      }
      
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
  if (i %in% 1:2) {
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
    names(pars_temp)[5:11] <- unit_estpars
    
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

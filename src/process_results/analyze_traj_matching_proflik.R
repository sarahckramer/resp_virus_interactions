# ---------------------------------------------------------------------------------------------------------------------
# Format and plot results from "round 2" of trajectory matching (using mega-likelihood)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Get/format date (for saving results):
date <- format(Sys.Date(), '%d%m%y')

# ---------------------------------------------------------------------------------------------------------------------

# Function to read in and format results:
load_and_format_proflik_results <- function(filepath, prof_par, shared_estpars) {
  
  # Remove prof_par from shared_estpars:
  shared_estpars <- shared_estpars[!(shared_estpars == prof_par)]
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full <- list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get estimated parameter and log-likelihood values:
  res_temp <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    select(all_of(shared_estpars)) %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols(map_chr(str_split(res_files, '_'), 5),
              map_chr(str_split(res_files, '_'), 7),
              map_chr(str_split(map_chr(str_split(res_files, '_'), 8), fixed('.')), 1)) %>%
    rename(vir1 = '...6',
           profpar = '...7',
           run = '...8') %>%
    mutate(profpar = as.numeric(profpar),
           run = as.numeric(run),
           vir1 = if_else(vir1 == 'b', 'B', 'H1'),
           vir1 = factor(vir1),
           vir1 = relevel(vir1, ref = 'H1')) %>%
    arrange(vir1, profpar, run)
  expect_true(nrow(res_temp) == length(res_files))
  expect_true(all(is.finite(res_temp$loglik)))
  
  # Set profpar to correct values:
  if (prof_par == 'delta') {
    res_temp <- res_temp %>%
      mutate(profpar = (7 / seq(5, 255, by = 5))[profpar]) %>%
      mutate(profpar = 7 / profpar)
  } else {
    res_temp <- res_temp %>%
      mutate(profpar = seq(0.0, 1.0, by = 0.02)[profpar])
  }
  
  # # Rename profpar column:
  # res_temp <- res_temp %>%
  #   rename(!!prof_par := 'profpar')
  
  # Return formatted results:
  return(res_temp)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format results for all runs

# Set shared estimated parameters:
shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1', 'theta_lambda2')

# Read in and format results:
res_thetalambda1 <- load_and_format_proflik_results(filepath = 'results/prof_lik_thetalambda1/',
                                                    prof_par = 'theta_lambda1',
                                                    shared_estpars = shared_estpars)
res_thetalambda2 <- load_and_format_proflik_results(filepath = 'results/prof_lik_thetalambda2/',
                                                    prof_par = 'theta_lambda2',
                                                    shared_estpars = shared_estpars)
res_delta <- load_and_format_proflik_results(filepath = 'results/prof_lik_delta/',
                                             prof_par = 'delta',
                                             shared_estpars = shared_estpars)

# Combine all results?:
res_list <- list(res_thetalambda1, res_thetalambda2, res_delta)
names(res_list) <- c('theta_lambda1', 'theta_lambda2', 'delta')
rm(res_thetalambda1, res_thetalambda2, res_delta)

# ---------------------------------------------------------------------------------------------------------------------

# Plot profiles and 99% CIs

# Get maximum log-likelihood values for each virus/profile:
maxloglik <- lapply(res_list, function(ix) {
  ix %>%
    group_by(vir1) %>%
    summarise(loglik = max(loglik)) %>%
    pull(loglik) %>%
    unlist()
})

# Calculate cutoff value for 99% CI and add to tibbles:
ci_cutoff <- lapply(maxloglik, function(ix) {
  ix - 0.5 * qchisq(df = 1, p = 0.99)
})
res_list <- lapply(1:length(res_list), function(ix) {
  res_list[[ix]] %>% mutate(ci = if_else(vir1 == 'B', ci_cutoff[[ix]][2], ci_cutoff[[ix]][1]))
})
names(res_list) <- c('theta_lambda1', 'theta_lambda2', 'delta')

# Get top estimate for each value of profpar:
res_list_top <- lapply(res_list, function(ix) {
  ix %>%
    group_by(vir1, profpar) %>%
    filter(rank(-loglik) == 1) %>%
    ungroup()
})
names(res_list_top) <- c('theta_lambda1', 'theta_lambda2', 'delta')

# Plot profile likelihoods with cutoff for 99% CI:
plot_list <- lapply(1:length(res_list_top), function(ix) {
  ggplot(res_list_top[[ix]], aes(x = profpar, y = loglik)) +
    geom_point() + theme_classic() +
    facet_wrap(~ vir1, scales = 'free_y') +
    geom_smooth(method = 'loess', span = 0.25) +
    geom_hline(color = 'red', aes(yintercept = ci)) +
    labs(x = names(res_list_top)[ix], y = 'Log Likelihood')
})

# ---------------------------------------------------------------------------------------------------------------------

# How do other shared parameter values change as prof_par changes?

for (i in 1:length(res_list)) {
  res_list[[i]] %>%
    filter(vir1 == 'H1') %>%
    select(all_of(shared_estpars[shared_estpars != names(res_list)[i]]), profpar) %>%
    pairs(pch = 20)
  
  res_list[[i]] %>%
    filter(vir1 == 'B') %>%
    select(all_of(shared_estpars[shared_estpars != names(res_list)[i]]), profpar) %>%
    pairs(pch = 20)
}

# ---------------------------------------------------------------------------------------------------------------------

# Save plots to file:
pdf(paste0('results/plots/', date, '_prof_lik_tj.pdf'), width = 10, height = 4.5)
print(plot_list)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

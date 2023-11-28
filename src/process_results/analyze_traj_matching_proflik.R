# ---------------------------------------------------------------------------------------------------------------------
# Format and plot results from "round 2" of trajectory matching (using mega-likelihood)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(patchwork)
library(testthat)

# Set directory where profile likelihood results are stored:
res_dir <- 'results/prof_lik_thetalambda1/'

# Check that directory for storing plots exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/plots/')) {
  dir.create('results/plots')
}

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
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist()) %>%
    bind_cols(map_chr(str_split(res_files, '_'), 10),
              paste0('0.', map_chr(str_split(map_chr(str_split(res_files, '_'), 11), fixed('.')), 2))) %>%
    rename(run = '...14',
           profpar = '...15') %>%
    mutate(run = as.numeric(run),
           profpar = as.numeric(profpar)) %>%
    arrange(profpar, run)
  expect_true(nrow(res_temp) == length(res_files))
  expect_true(all(is.finite(res_temp$loglik)))
  
  res_temp <- res_temp %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # Return formatted results:
  return(res_temp)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format results for all runs

# Set shared estimated parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')

# Read in and format results:
res <- load_and_format_proflik_results(filepath = res_dir,
                                       prof_par = 'theta_lambda1',
                                       shared_estpars = shared_estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Plot profiles and 95% CIs

# Get maximum log-likelihood values for each virus/profile:
maxloglik <- res %>%
  summarise(loglik = max(loglik)) %>%
  pull(loglik)

# Calculate cutoff value for 95% CI and add to tibbles:
ci_cutoff <- maxloglik - 0.5 * qchisq(df = 1, p = 0.95)
res <- res %>%
  mutate(ci = ci_cutoff)

# Get top estimate for each value of profpar:
res <- res %>%
  group_by(profpar) %>%
  filter(rank(-loglik) == 1) %>%
  ungroup()

# Plot profile likelihoods with cutoff for 99% CI:
p1 <- ggplot(res, # filter(loglik > (max(loglik) - 100)), # allows for zooming
             aes(x = profpar, y = loglik)) +
  geom_point() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  geom_smooth(method = 'loess', span = 0.75, color = 'black') +
  geom_hline(color = 'black', aes(yintercept = ci), linewidth = 1, lty = 2) +
  labs(x = bquote(theta[lambda*1]), y = 'Log-Likelihood') +
  scale_x_continuous(n.breaks = 10)
print(p1)

# ---------------------------------------------------------------------------------------------------------------------

# How do other shared parameter values change as prof_par changes?
res %>%
  select(all_of(shared_estpars[shared_estpars != 'theta_lambda1']), profpar) %>%
  pairs(pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

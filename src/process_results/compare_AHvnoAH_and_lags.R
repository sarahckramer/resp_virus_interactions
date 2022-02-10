# ---------------------------------------------------------------------------------------------------------------------
# Format and plot results from "round 2" of trajectory matching (using mega-likelihood)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# ---------------------------------------------------------------------------------------------------------------------

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, cond) {
  
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
  pars_top <- pars_df[1:no_best, ]
  
  # Add label:
  pars_top <- pars_top %>%
    mutate(condition = cond)
  
  # Return formatted results:
  return(pars_top)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Read in results for all runs

# For "full" models, lag 0:
res_h1_full <- load_and_format_mega_results('results/round2_2_fluH1_FULL/', cond = 'full')
res_b_full <- load_and_format_mega_results('results/round2_2_fluB_FULL/', cond = 'full')

# For models w/o AH:
res_h1_noAH <- load_and_format_mega_results('results/round2_3_fluH1_noAH/', cond = 'noAH')
res_b_noAH <- load_and_format_mega_results('results/round2_2_fluB_noAH/', cond = 'noAH')

# For lagged models:
res_h1_lag1 <- load_and_format_mega_results('results/round2_2_fluH1_FULL_lag1/', cond = 'lag1')
res_h1_lag2 <- load_and_format_mega_results('results/round2_2_fluH1_FULL_lag2/', cond = 'lag2')

res_b_lag1 <- load_and_format_mega_results('results/round2_2_fluB_FULL_lag1/', cond = 'lag1')
res_b_lag2 <- load_and_format_mega_results('results/round2_2_fluB_FULL_lag2/', cond = 'lag2')

# ---------------------------------------------------------------------------------------------------------------------

# Compare runs with and without eta_ah1/eta_ah2

# Compare parameter estimates:
summary(res_h1_full %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_h1_noAH %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))

summary(res_b_full %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_b_noAH %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))

# Compare log likelihoods:
res_h1 <- bind_rows(res_h1_full, res_h1_noAH) %>%
  select(rho1:eta_ah2, loglik:condition)
res_b <- bind_rows(res_b_full, res_b_noAH) %>%
  select(rho1:eta_ah2, loglik:condition)

p1 <- ggplot(data = res_h1, aes(x = condition, y = loglik, group = condition)) + geom_jitter() + theme_classic()# + geom_boxplot()
p2 <- ggplot(data = res_b, aes(x = condition, y = loglik, group = condition)) + geom_jitter() + theme_classic()# + geom_boxplot()
grid.arrange(p1, p2, ncol = 2)

# full is significantly better than noAH if 2 * (loglik_full - loglik_noAH) <= qchisq(p = 0.95, df = 2 - 1)
qchisq(p = 0.95, df = 2 - 1)

2 * (min(res_h1$loglik[res_h1$condition == 'full']) - max(res_h1$loglik[res_h1$condition == 'noAH']))
2 * (min(res_b$loglik[res_b$condition == 'full']) - max(res_b$loglik[res_b$condition == 'noAH']))

# ---------------------------------------------------------------------------------------------------------------------

# Compare runs with various lags on climate variables

# Compare parameter estimates:
summary(res_h1_full %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_h1_lag1 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_h1_lag2 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))

summary(res_b_full %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_b_lag1 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))
summary(res_b_lag2 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2')))

# Compare log likelihoods:
res_h1 <- bind_rows(res_h1_full, res_h1_lag1, res_h1_lag2) %>%
  select(rho1:eta_ah2, loglik:condition)
res_b <- bind_rows(res_b_full, res_b_lag1, res_b_lag2) %>%
  select(rho1:eta_ah2, loglik:condition)

p1 <- ggplot(data = res_h1, aes(x = condition, y = loglik, group = condition)) + geom_jitter() + theme_classic()# + geom_boxplot()
p2 <- ggplot(data = res_b, aes(x = condition, y = loglik, group = condition)) + geom_jitter() + theme_classic()# + geom_boxplot()
grid.arrange(p1, p2, ncol = 2)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up
rm(list = ls())

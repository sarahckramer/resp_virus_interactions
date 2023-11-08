# ---------------------------------------------------------------------------------------------------------------------
# Check whether including effect on absolute humidity significantly improves model fit
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Set directory where results are stored:
res_dir_h1 <- 'results/round2_4_fluH1_FULL/'
res_dir_b <- 'results/round2_3_fluB_FULL/'

res_dir_h1_noAH <- 'results/round2_4_fluH1_noAH/'
res_dir_b_noAH <- 'results/round2_3_fluB_noAH/'

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
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  if (cond == 'full') {
    expect_equal(df_use, 54)
  } else if (cond == 'noAH') {
    expect_equal(df_use, 52)
  }

  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # Add label:
  pars_top <- pars_top %>%
    mutate(condition = cond)
  
  # Return formatted results:
  return(pars_top)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Read in results for all runs

# For "full" models:
res_h1_full <- load_and_format_mega_results(res_dir_h1, cond = 'full')
res_b_full <- load_and_format_mega_results(res_dir_b, cond = 'full')

# For models w/o AH:
res_h1_noAH <- load_and_format_mega_results(res_dir_h1_noAH, cond = 'noAH')
res_b_noAH <- load_and_format_mega_results(res_dir_b_noAH, cond = 'noAH')

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

# Check for significance:
# full is significantly better than noAH if 2 * (loglik_full - loglik_noAH) > qchisq(p = 0.95, df = 2)
print(2 * (min(res_h1$loglik[res_h1$condition == 'full']) - max(res_h1$loglik[res_h1$condition == 'noAH'])) > qchisq(p = 0.95, df = 2))
print(2 * (min(res_b$loglik[res_b$condition == 'full']) - max(res_b$loglik[res_b$condition == 'noAH'])) > qchisq(p = 0.95, df = 2))

# ---------------------------------------------------------------------------------------------------------------------

# Clean up
rm(list = ls())

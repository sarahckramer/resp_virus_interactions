# ---------------------------------------------------------------------------------------------------------------------
# Get start ranges for from round2 trajectory matching results (for profile likelihoods, or to rerun round2)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where results from round2 fits are stored:
res_dir_h1 <- 'results/round2_1_fluH1/'
res_dir_b <- 'results/round2_3_fluB_FULL/'

# Are results from a sensitivity analysis?:
sens <- 'main'

# Check that directory for storing results exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/round2_CIs/')) {
  dir.create('results/round2_CIs/')
}

which_round_h1 <- str_split(res_dir_h1, '_')[[1]][3] # which round's results are we using here?
which_round_b <- str_split(res_dir_b, '_')[[1]][3]

if (sens == 'main') {
  
  new_dir_h1 <- paste0('results/round2_CIs/from_2_', which_round_h1, '/')
  if (!dir.exists(new_dir_h1)) {
    dir.create(new_dir_h1)
  }
  new_dir_b <- paste0('results/round2_CIs/from_2_', which_round_b, '/')
  if (!dir.exists(new_dir_b)) {
    dir.create(new_dir_b)
  }
  
} else {
  
  if(!dir.exists('results/round2_CIs/sens/')) {
    dir.create('results/round2_CIs/sens/')
  }
  if(!dir.exists(paste0('results/round2_CIs/sens/', sens, '/'))) {
    dir.create(paste0('results/round2_CIs/sens/', sens, '/'))
  }
  
  new_dir_h1 <- paste0('results/round2_CIs/sens/', sens, '/from_2_', which_round_h1, '/')
  if (!dir.exists(new_dir_h1)) {
    dir.create(new_dir_h1)
  }
  new_dir_b <- paste0('results/round2_CIs/sens/', sens, '/from_2_', which_round_b, '/')
  if (!dir.exists(new_dir_b)) {
    dir.create(new_dir_b)
  }
  
}

# Function to read in and format results:
load_and_format_mega_results <- function(filepath) {
  
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
  # expect_equal(df_use, 54)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  
  # If only one parameter set is in this range, MLE has not yet been reached; take top 5% of fits instead:
  if (no_best == 1) {
    print('MLE not reached!')
    no_best <- 25
  }
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # # Remove where d2 > 10 and theta_lambda2 != 1.0:
  # pars_top <- pars_top %>%
  #   filter(!(d2 > 10.0 & theta_lambda2 < 0.99))
  
  # Set unrealistic values to NA:
  pars_top$delta1[pars_top$delta1 > 7.0] <- NA
  pars_top$d2[pars_top$d2 > 10.0] <- NA
  pars_top$rho1[pars_top$rho1 == 1.0] <- NA
  pars_top$rho2[pars_top$rho2 == 1.0] <- NA
  
  # Since phi=0 is equivalent to phi=52.25, don't use full range; transform so that we can select from only best-supported range:
  hist(pars_top$phi, breaks = 50)
  pars_top <- pars_top %>%
    mutate(phi = if_else(phi < 26, phi + 52.25, phi))
  
  # Return formatted results:
  return(pars_top)
  
}

# Read in results:
res_h1 <- load_and_format_mega_results(res_dir_h1) %>%
  select(-loglik)
res_b <- load_and_format_mega_results(res_dir_b) %>%
  select(-loglik)

# Get minimum and maximum start values:
ci_start_h1 <- as.data.frame(rbind(summarise(res_h1, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                   summarise(res_h1, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))
ci_start_b <- as.data.frame(rbind(summarise(res_b, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                  summarise(res_b, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))

# Possible that d2 ranges are missing if all top fits were > 10; if so, replace:
if (any(ci_start_h1 == Inf)) {
  ci_start_h1$d2 <- c(0, 10)
}
if (any(ci_start_b == Inf)) {
  ci_start_b$d2 <- c(0, 10)
}

# Check that sums of initial conditions can't sum to >1:
init_cond_estpars <- c('I10', 'I20', 'R10', 'R20', 'R120')

sums_h1 <- ci_start_h1 %>%
  mutate(minmax = c('min', 'max')) %>%
  select(contains(init_cond_estpars), minmax) %>%
  pivot_longer(-minmax) %>%
  mutate(season = str_sub(name, 1, 6)) %>%
  group_by(season, minmax) %>%
  summarise(sum = sum(value))
sums_b <- ci_start_b %>%
  mutate(minmax = c('min', 'max')) %>%
  select(contains(init_cond_estpars), minmax) %>%
  pivot_longer(-minmax) %>%
  mutate(season = str_sub(name, 1, 6)) %>%
  group_by(season, minmax) %>%
  summarise(sum = sum(value))

if (any(sums_h1 %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}
if (any(sums_b %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}

if (any(sums_h1 %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  seasons <- sums_h1 %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    pull(season)
  
  for (yr in seasons) {
    
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start_h1[2, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    red_needed <- sums_h1 %>%
      filter(season == yr,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start_h1[1, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start_h1[2, which(str_detect(names(ci_start_h1), yr) &
                           (str_detect(names(ci_start_h1), 'R10') |
                              str_detect(names(ci_start_h1), 'R20') |
                              str_detect(names(ci_start_h1), 'R120')))] <- new_upper_bounds
    expect_lt(ci_start_h1[2, ] %>% select(contains(yr)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(yr)
}

if (any(sums_b %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  seasons <- sums_b %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    pull(season)
  
  for (yr in seasons) {
    
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start_b[2, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    red_needed <- sums_b %>%
      filter(season == yr,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start_b[1, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start_b[2, which(str_detect(names(ci_start_b), yr) &
                          (str_detect(names(ci_start_b), 'R10') |
                             str_detect(names(ci_start_b), 'R20') |
                             str_detect(names(ci_start_b), 'R120')))] <- new_upper_bounds
    expect_lt(ci_start_b[2, ] %>% select(contains(yr)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(yr)
}

# Write start ranges to file:
write_rds(ci_start_h1, file = paste0(new_dir_h1, 'round2CI_startvals_H1.rds'))
write_rds(ci_start_b, file = paste0(new_dir_b, 'round2CI_startvals_B.rds'))

# Also find start ranges for H1+B sensitivity analysis?
res_dir_h1_plus_b <- 'results/round2_fit/sens/round2_1_fluH1_plus_B/'
which_round_h1_plus_b <- str_split(res_dir_h1_plus_b, '_')[[1]][3]

if(!dir.exists('results/round2_CIs/sens/')) {
  dir.create('results/round2_CIs/sens/')
}
if(!dir.exists('results/round2_CIs/sens/flu_h1_plus_b/')) {
  dir.create('results/round2_CIs/sens/flu_h1_plus_b/')
}

new_dir_h1_plus_b <- paste0('results/round2_CIs/sens/flu_h1_plus_b/from_2_', which_round_h1_plus_b, '/')
if (!dir.exists(new_dir_h1_plus_b)) {
  dir.create(new_dir_h1_plus_b)
}

res_h1_plus_b <- load_and_format_mega_results(res_dir_h1_plus_b) %>%
  select(-loglik)
ci_start_h1_plus_b <- as.data.frame(rbind(summarise(res_h1_plus_b, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                          summarise(res_h1_plus_b, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))
if (any(ci_start_h1_plus_b == Inf)) {
  ci_start_h1_plus_b$d2 <- c(0, 10)
}

sums_h1_plus_b <- ci_start_h1_plus_b %>%
  mutate(minmax = c('min', 'max')) %>%
  select(contains(init_cond_estpars), minmax) %>%
  pivot_longer(-minmax) %>%
  mutate(season = str_sub(name, 1, 6)) %>%
  group_by(season, minmax) %>%
  summarise(sum = sum(value))
if (any(sums_h1_plus_b %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}

if (any(sums_h1_plus_b %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  seasons <- sums_h1_plus_b %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    pull(season)
  
  for (yr in seasons) {
    
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start_h1_plus_b[2, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    red_needed <- sums_h1_plus_b %>%
      filter(season == yr,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start_h1_plus_b[1, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start_h1_plus_b[2, which(str_detect(names(ci_start_h1_plus_b), yr) &
                                  (str_detect(names(ci_start_h1_plus_b), 'R10') |
                                     str_detect(names(ci_start_h1_plus_b), 'R20') |
                                     str_detect(names(ci_start_h1_plus_b), 'R120')))] <- new_upper_bounds
    expect_lt(ci_start_h1_plus_b[2, ] %>% select(contains(yr)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(yr)
}

write_rds(ci_start_h1_plus_b, file = paste0(new_dir_h1_plus_b, 'round2CI_startvals_H1_plus_B.rds'))

# Clean up:
rm(list = ls())

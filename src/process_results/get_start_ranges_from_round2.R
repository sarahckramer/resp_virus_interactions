# ---------------------------------------------------------------------------------------------------------------------
# Get start ranges for from round2 trajectory matching results (for profile likelihoods, bootstraps, or to rerun round2)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where results from round2 fits are stored:
res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'

# Which round of fits?:
which_round <- str_split(res_dir, '_')[[1]][which(!is.na(as.numeric(str_split(res_dir, '_')[[1]])))]

# Are results from a sensitivity analysis?:
if (str_detect(res_dir, 'sinusoidal')) {
  sens <- 'sinusoidal_forcing'
} else if (str_detect(res_dir, 'h3_covar')) {
  sens <- 'h3_covar'
} else if (str_detect(res_dir, 'less_circ_h3')) {
  sens <- 'less_circ_h3'
} else if (str_detect(res_dir, 'no_rsv_immune')) {
  sens <- 'no_rsv_immune'
} else if (str_detect(res_dir, 'no_ah')) {
  sens <- 'no_ah'
} else if (str_detect(res_dir, 'rhino_covar')) {
  sens <- 'rhino_covar'
} else {
  sens <- 'main'
}

# Check that directory for storing results exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/round2_CIs/')) {
  dir.create('results/round2_CIs/')
}

if (sens == 'main') {
  
  new_dir <- paste0('results/round2_CIs/from_2_', which_round, '/')
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }
  
  if (str_detect(res_dir, 'age_structured')) {
    
    if(!dir.exists('results/round2_CIs/sens/')) {
      dir.create('results/round2_CIs/sens/')
    }
    if(!dir.exists('results/round2_CIs/sens/age_structured/')) {
      dir.create('results/round2_CIs/sens/age_structured/')
    }
    
    new_dir <- paste0('results/round2_CIs/sens/age_structured/from_2_', which_round, '/')
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
  }
  
} else {
  
  if(!dir.exists('results/round2_CIs/sens/')) {
    dir.create('results/round2_CIs/sens/')
  }
  if(!dir.exists(paste0('results/round2_CIs/sens/', sens, '/'))) {
    dir.create(paste0('results/round2_CIs/sens/', sens, '/'))
  }
  
  new_dir <- paste0('results/round2_CIs/sens/', sens, '/from_2_', which_round, '/')
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
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
  print(table(pars_df$message))
  
  # If only one parameter set is in this range, MLE has not yet been reached; take top 5% of fits instead:
  if (no_best == 1) {
    print('MLE not reached!')
    no_best <- 25
  }
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  print(summary(pars_top$loglik))
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # If none remaining, take top 25 again:
  while(nrow(pars_top) == 0) {
    no_best <- max(no_best, 25)
    
    pars_top <- pars_df[1:no_best, ] %>%
      filter(!str_detect(message, 'maxtime')) %>%
      select(-message)
    
    no_best <- no_best + 25
  }
  no_best <- no_best - 25
  print(no_best)
  print(summary(pars_top$loglik))
  
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
res <- load_and_format_mega_results(res_dir) %>%
  select(-loglik)

# Get minimum and maximum start values:
ci_start <- as.data.frame(rbind(summarise(res, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                summarise(res, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))

# Possible that d2 ranges are missing if all top fits were > 10; if so, replace:
if (any(ci_start == Inf)) {
  ci_start$d2 <- c(0, 10)
}

# Check that sums of initial conditions can't sum to >1:
init_cond_estpars <- c('I10', 'I20', 'R10', 'R20', 'R120')

sums <- ci_start %>%
  mutate(minmax = c('min', 'max')) %>%
  select(contains(init_cond_estpars), minmax) %>%
  pivot_longer(-minmax) %>%
  mutate(season = str_sub(name, 1, 6)) %>%
  group_by(season, minmax) %>%
  summarise(sum = sum(value))

if (any(sums %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}

if (any(sums %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  seasons <- sums %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    pull(season)
  
  for (yr in seasons) {
    
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start[2, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    red_needed <- sums %>%
      filter(season == yr,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start[1, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    
    if (!all(new_upper_bounds > orig_lower_bounds)) {
      new_upper_bounds_try <- orig_upper_bounds
      new_upper_bounds_try[-which(orig_lower_bounds >= new_upper_bounds)] <- orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] - (red_needed * (orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] / sum(orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)])))
      new_upper_bounds <- new_upper_bounds_try
      rm(new_upper_bounds_try)
    }
    
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start[2, which(str_detect(names(ci_start), yr) &
                        (str_detect(names(ci_start), 'R10') |
                           str_detect(names(ci_start), 'R20') |
                           str_detect(names(ci_start), 'R120')))] <- new_upper_bounds
    expect_lt(ci_start[2, ] %>% select(contains(yr)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(yr)
}

# Write start ranges to file:
write_rds(ci_start, file = paste0(new_dir, 'round2CI_startvals_H1_plus_B.rds'))

# Clean up:
rm(list = ls())

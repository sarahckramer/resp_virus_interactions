# ---------------------------------------------------------------------------------------------------------------------
# Get start ranges for from round2 trajectory matching results (for profile likelihoods, or to rerun round2)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Use 95% CI or top 5% of fits?:
method <- 'perc' # 'ci' or 'perc'

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, method) {
  
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
  
  if (method == 'ci') {
    no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
    # no_best <- max(no_best, 50)
  } else if (method == 'perc') {
    no_best <- 25
  }
  
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
res_h1 <- load_and_format_mega_results('results/round2_4_fluH1_FULL/', method) %>%
  select(-loglik)
res_b <- load_and_format_mega_results('results/round2_3_fluB_FULL/', method) %>%
  select(-loglik)

# Get minimum and maximum start values:
ci_start_h1 <- as.data.frame(rbind(summarise(res_h1, across(.cols = everything(), min, na.rm = TRUE)),
                                   summarise(res_h1, across(.cols = everything(), max, na.rm = TRUE))))
ci_start_b <- as.data.frame(rbind(summarise(res_b, across(.cols = everything(), min, na.rm = TRUE)),
                                  summarise(res_b, across(.cols = everything(), max, na.rm = TRUE))))

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
write_rds(ci_start_h1, file = 'results/round2_CIs/round2CI_startvals_H1.rds')
write_rds(ci_start_b, file = 'results/round2_CIs/round2CI_startvals_B.rds')

# Clean up:
rm(list = ls())

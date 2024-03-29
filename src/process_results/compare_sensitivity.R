# ---------------------------------------------------------------------------------------------------------------------
# Check whether including effect on absolute humidity significantly improves model fit
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Set directories where results are stored:
res_dir_main <- 'results/round2_fit/round2_3_fluH1_plus_B/'
res_dir_noAH <- 'results/round2_fit/sens/no_ah/round2_3_fluH1_plus_B/'
res_dir_sinusoidal <- 'results/round2_fit/sens/sinusoidal_forcing/round2_3_fluH1_plus_B/'
res_dir_noint <- 'results/round2_fit/sens/no_int/round2_2_fluH1_plus_B/'
res_dir_noRSVimmune <- 'results/round2_fit/sens/no_rsv_immune/round2_5_fluH1_plus_B/'
res_dir_h3covar <- 'results/round2_fit/sens/h3_covar/round2_3_fluH1_plus_B/'
res_dir_h3covar_lag1 <- 'results/round2_fit/sens/h3_covar/lag1/round2_3_fluH1_plus_B/'
res_dir_h3covar_lag2 <- 'results/round2_fit/sens/h3_covar/lag2/round2_3_fluH1_plus_B/'
res_dir_h3covar_lag15 <- 'results/round2_fit/sens/h3_covar/lag15/round2_3_fluH1_plus_B/'
res_dir_lesscirch3 <- 'results/round2_fit/sens/less_circ_h3/round2_5_fluH1_plus_B/'
res_dir_rhino <- 'results/round2_fit/sens/rhino_covar/round2_3_fluH1_plus_B/'
res_dir_suscplussev <- 'results/round2_fit/sens/susc_plus_sev/round2_4_fluH1_plus_B/'

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
  if (str_detect(res_files[1], 'PARALLEL')) {
    
    res_full <- do.call('c', res_full)
    num_errors <- length(which(res_full == 'error'))
    
    if (num_errors > 0) {
      res_full <- res_full[-which(res_full == 'error')]
    }
    
  }
  
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  
  if (str_detect(res_files[1], 'PARALLEL')) {
    expect_true(nrow(pars_df) == (length(res_files) * 25) - num_errors)
  } else {
    expect_true(nrow(pars_df) == length(res_files))
  }
  
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  
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

res_main <- load_and_format_mega_results(res_dir_main, cond = 'Main')
res_noah <- load_and_format_mega_results(res_dir_noAH, cond = 'No AH')
res_sinusoidal <- load_and_format_mega_results(res_dir_sinusoidal, cond = 'Sinusoidal Forcing')
res_noint <- load_and_format_mega_results(res_dir_noint, cond = 'No Interaction')
res_noRSVimmune <- load_and_format_mega_results(res_dir_noRSVimmune, cond = 'No Immunity to RSV')
res_h3covar <- load_and_format_mega_results(res_dir_h3covar, cond = 'H3 as Covariate')
res_h3covar_lag1 <- load_and_format_mega_results(res_dir_h3covar_lag1, cond = 'H3 as Covariate (Lag 1)')
res_h3covar_lag2 <- load_and_format_mega_results(res_dir_h3covar_lag2, cond = 'H3 as Covariate (Lag 2)')
res_h3covar_lag15 <- load_and_format_mega_results(res_dir_h3covar_lag15, cond = 'H3 as Covariate (Lag 15)')
res_lesscirch3 <- load_and_format_mega_results(res_dir_lesscirch3, cond = 'Low H3 Circulation Seasons')
res_rhino <- load_and_format_mega_results(res_dir_rhino, cond = 'Rhinovirus as Covariate')
res_suscplussev <- load_and_format_mega_results(res_dir_suscplussev, cond = 'Incl. Severity')

# ---------------------------------------------------------------------------------------------------------------------

# Compare runs

# Compare parameter estimates:
summary(res_main %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_noah %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_sinusoidal %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
# summary(res_noint %>%
#           select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_noRSVimmune %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_h3covar %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_h3covar_lag1 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_h3covar_lag2 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_h3covar_lag15 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_lesscirch3 %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_rhino %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))
summary(res_suscplussev %>%
          select(!contains('I10') & !contains('I20') & !contains('Ri') & !contains('R1') & !contains('R2'), -c(loglik, condition)))

# Compare log likelihoods:
res <- bind_rows(res_main,
                 res_noah,
                 res_sinusoidal,
                 res_noint,
                 res_noRSVimmune,
                 res_h3covar,
                 res_h3covar_lag1,
                 res_h3covar_lag2,
                 res_h3covar_lag15,
                 res_rhino,
                 res_suscplussev)

p1 <- ggplot(data = res, aes(x = condition, y = loglik, group = condition)) + geom_jitter() + theme_classic()# + geom_boxplot()
print(p1)

# Check for significance:
# full is significantly better than noAH if 2 * (loglik_full - loglik_noAH) > qchisq(p = 0.95, df = 2)
print(2 * (min(res$loglik[res$condition == 'Main']) - max(res$loglik[res$condition == 'No AH'])) > qchisq(p = 0.95, df = 2))
# print(2 * (min(res$loglik[res$condition == 'Main']) - max(res$loglik[res$condition == 'Sinusoidal Forcing'])) > qchisq(p = 0.95, df = 0))
print(2 * (min(res$loglik[res$condition == 'Main']) - max(res$loglik[res$condition == 'No Interaction'])) > qchisq(p = 0.95, df = 4))
print(2 * (min(res$loglik[res$condition == 'Main']) - max(res$loglik[res$condition == 'No Immunity to RSV'])) > qchisq(p = 0.95, df = 12))
print(2 * (min(res$loglik[res$condition == 'H3 as Covariate']) - max(res$loglik[res$condition == 'Main'])) > qchisq(p = 0.95, df = 1))
print(2 * (min(res$loglik[res$condition == 'Rhinovirus as Covariate']) - max(res$loglik[res$condition == 'Main'])) > qchisq(p = 0.95, df = 1))
print(2 * (min(res$loglik[res$condition == 'Main']) - max(res$loglik[res$condition == 'Incl. Severity'])) > qchisq(p = 0.95, df = 2))

aic_main <- 2 * length(names(res_main %>% select(-c(loglik, condition)))) - 2 * max(res_main$loglik)
aic_noah <- 2 * length(names(res_noah %>% select(-c(loglik, condition)))) - 2 * max(res_noah$loglik)
aic_sinusoidal <- 2 * length(names(res_sinusoidal %>% select(-c(loglik, condition)))) - 2 * max(res_sinusoidal$loglik)
aic_noint <- 2 * length(names(res_noint %>% select(-c(loglik, condition)))) - 2 * max(res_noint$loglik)
aic_noRSVimmune <- 2 * length(names(res_noRSVimmune %>% select(-c(loglik, condition)))) - 2 * max(res_noRSVimmune$loglik)
aic_h3covar <- 2 * length(names(res_h3covar %>% select(-c(loglik, condition)))) - 2 * max(res_h3covar$loglik)
aic_suscplussev <- 2 * length(names(res_suscplussev %>% select(-c(loglik, condition)))) - 2 * max(res_suscplussev$loglik)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up
rm(list = ls())

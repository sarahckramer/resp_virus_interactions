# ---------------------------------------------------------------------------------------------------------------------
# Check adequacy of 500 samples from starting parameter ranges
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Function to load and format results:
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
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(table(pars_df$message))
  
  # If only one parameter set is in this range, MLE has not yet been reached; take top 5% of fits instead:
  if (no_best == 1) {
    print('MLE not reached!')
    no_best <- 25
  }
  print(no_best)
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # If none remaining, print warning:
  if(nrow(pars_top) == 0) {
    print('No convergence among best-fit runs!')
    
    no_best <- 25
    pars_top <- pars_df[1:no_best, ]
    
    pars_top <- pars_top %>%
      filter(!str_detect(message, 'maxtime')) %>%
      select(-message)
  }
  
  # Return formatted results:
  return(pars_top)
  
}

# Read in top parameter sets:
res_dir_500 <- 'results/round2_fit/round2_3_fluH1_plus_B/'
res_dir_2500 <- 'results/round2_fit/sens/increase_sobol_size/round2_4_fluH1_plus_B/'

res_500 <- load_and_format_mega_results(res_dir_500) %>%
  mutate(sobol_size = 500)
res_2500 <- load_and_format_mega_results(res_dir_2500) %>%
  mutate(sobol_size = 2500)

# mle_500 <- read_rds('results/MLEs_flu_h1_plus_b.rds')
# mle_2500 <- read_rds('results/round2_fit/sens/increase_sobol_size/MLEs_flu_h1_plus_b.rds')

# Combine:
res <- bind_rows(res_500, res_2500) %>%
  mutate(sobol_size = factor(sobol_size, levels = c('500', '2500')))
res_mle <- bind_rows(res_500[1, ], res_2500[1, ]) %>%
  mutate(sobol_size = factor(sobol_size, levels = c('500', '2500')))

# Format:
res_ll <- res %>%
  select(loglik:sobol_size)
res_mle <- res_mle %>%
  select(rho1:eta_ah2, sobol_size) %>%
  pivot_longer(rho1:eta_ah2, names_to = 'parameter')
res <- res %>%
  select(rho1:eta_ah2, sobol_size) %>%
  pivot_longer(rho1:eta_ah2, names_to = 'parameter')

# Compare log-likelihoods:
p1 <- ggplot(data = res_ll) +
  geom_violin(aes(x = sobol_size, y = loglik), fill = 'gray95') +
  theme_classic() +
  labs(x = '# of Random Start Sets', y = 'Log-Likelihood')

# Compare shared parameter estimates:
p2 <- ggplot() +
  geom_violin(data = res, aes(x = sobol_size, y = value), fill = 'gray95') +
  geom_point(data = res_mle, aes(x = sobol_size, y = value), shape = '\u2605', size = 4) +
  facet_wrap(~ parameter, scales = 'free_y') +
  theme_classic() +
  labs(x = '# of Random Start Sets', y = 'Parameter Value')

# Plot:
plot(p1)
plot(p2)

# # Also look at round by round results?
# # Read in results (original - 500 starting sets):
# res_dir_500_r1 <- 'results/round2_fit/round2_1_fluH1_plus_B/'
# res_dir_500_r2 <- 'results/round2_fit/round2_2_fluH1_plus_B/'
# res_dir_500_r3 <- 'results/round2_fit/round2_3_fluH1_plus_B/'
# 
# res_500_r1 <- load_and_format_mega_results(res_dir_500_r1) %>%
#   mutate(round = 1)
# res_500_r2 <- load_and_format_mega_results(res_dir_500_r2) %>%
#   mutate(round = 2)
# res_500_r3 <- load_and_format_mega_results(res_dir_500_r3) %>%
#   mutate(round = 3)
# 
# # Read in results (expanded - 2500 starting sets):
# res_dir_2500_r1 <- 'results/round2_fit/sens/increase_sobol_size/round2_1_fluH1_plus_B/'
# res_dir_2500_r2 <- 'results/round2_fit/sens/increase_sobol_size/round2_2_fluH1_plus_B/'
# res_dir_2500_r3 <- 'results/round2_fit/sens/increase_sobol_size/round2_3_fluH1_plus_B/'
# res_dir_2500_r4 <- 'results/round2_fit/sens/increase_sobol_size/round2_4_fluH1_plus_B/'
# 
# res_2500_r1 <- load_and_format_mega_results(res_dir_2500_r1) %>%
#   mutate(round = 1)
# res_2500_r2 <- load_and_format_mega_results(res_dir_2500_r2) %>%
#   mutate(round = 2)
# res_2500_r3 <- load_and_format_mega_results(res_dir_2500_r3) %>%
#   mutate(round = 3)
# res_2500_r4 <- load_and_format_mega_results(res_dir_2500_r4) %>%
#   mutate(round = 4)

# Clean up:
rm(list = ls())

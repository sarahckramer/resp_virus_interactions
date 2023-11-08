# ---------------------------------------------------------------------------------------------------------------------
# Get MLEs from trajectory matching
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where final results from round2 fits are stored:
res_dir_h1 <- 'results/round2_4_fluH1_FULL/'
res_dir_b <- 'results/round2_3_fluB_FULL/'

# Check that directory for storing results exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
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
  
  print(dim(pars_df))
  print(names(pars_df))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  expect_equal(df_use, 54)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  print(table(pars_top$message))
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Return formatted results:
  return(pars_top)
  
}

# Get MLEs:
res_h1 <- load_and_format_mega_results(res_dir_h1) %>%
  select(-loglik)
res_b <- load_and_format_mega_results(res_dir_b) %>%
  select(-loglik)

# Save MLEs:
write_rds(res_h1, file = paste0('results/MLEs_flu_h1.rds'))
write_rds(res_b, file = paste0('results/MLEs_flu_b.rds'))

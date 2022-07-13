# ---------------------------------------------------------------------------------------------------------------------
# Get MLEs from trajectory matching
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

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
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
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
res_h1 <- load_and_format_mega_results('results/round2_4_fluH1_FULL/') %>%
  select(-loglik)
res_b <- load_and_format_mega_results('results/round2_3_fluB_FULL/') %>%
  select(-loglik)

# Save MLEs:
write_rds(res_h1, file = paste0('results/MLEs_flu_h1.rds'))
write_rds(res_b, file = paste0('results/MLEs_flu_b.rds'))

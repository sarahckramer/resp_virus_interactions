# ---------------------------------------------------------------------------------------------------------------------
# Get MLEs from trajectory matching
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where final results from round2 fits are stored:
res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'

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
  # expect_equal(df_use, 54)
  
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
res <- load_and_format_mega_results(res_dir) %>%
  select(-loglik)

# Save MLEs:
if (str_detect(res_dir, 'sens')) {
  
  write_rds(res, file = paste0(paste(str_split(res_dir, '/')[[1]][1:(length(str_split(res_dir, '/')[[1]]) - 2)], collapse = '/'), '/MLEs_flu_h1_plus_b.rds'))
  
} else {
  
  write_rds(res, file = 'results/MLEs_flu_h1_plus_b.rds')
  
}

# ---------------------------------------------------------------------------------------------------------------------
# Code to determine whether any results are missing after running cluster code
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set size of parameter start space:
sobol_size <- 500

# Get list of completed runs:
check_list_flu_h1 <- list.files(path = 'results/', pattern = 'flu_h1')
check_list_flu_b <- list.files(path = 'results/', pattern = 'flu_b')

# Get vector of seasons for which fitting was done:
all_yrs <- unique(str_sub(check_list_flu_h1, 16, 21))
all_yrs_check <- unique(str_sub(check_list_flu_b, 15, 20))
expect_true(all.equal(all_yrs, all_yrs_check))
rm(all_yrs_check)
print(all_yrs)
print(length(all_yrs))

# Check for complete results for flu_A/RSV:
if (length(check_list_flu_h1) != sobol_size * length(all_yrs)) {
  
  for (yr in all_yrs) {
    temp_list <- check_list_flu_h1[str_detect(check_list_flu_h1, pattern = yr)]
    
    if (length(temp_list) != sobol_size) {
      
      completed_runs <- str_split(temp_list, '_') %>%
        lapply(., function(ix) {ix[6]}) %>%
        unlist() %>%
        str_split(fixed('.')) %>%
        lapply(., function(ix) {ix[1]}) %>%
        unlist() %>%
        as.numeric()
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_h1/RSV in ', yr, ', missing run: ', run))
      }
      
    }
  }
}

# Check for complete results for flu_B/RSV:
if (length(check_list_flu_b) != sobol_size * length(all_yrs)) {
  
  for (yr in all_yrs) {
    temp_list <- check_list_flu_b[str_detect(check_list_flu_b, pattern = yr)]
    
    if (length(temp_list) != sobol_size) {
      
      completed_runs <- str_split(temp_list, '_') %>%
        lapply(., function(ix) {ix[6]}) %>%
        unlist() %>%
        str_split(fixed('.')) %>%
        lapply(., function(ix) {ix[1]}) %>%
        unlist() %>%
        as.numeric()
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_b/RSV in ', yr, ', missing run: ', run))
      }
      
    }
  }
}

# Clean up:
rm(list = ls())
print('Done.')

# ---------------------------------------------------------------------------------------------------------------------
# Code to determine whether any results are missing after running cluster code
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Set size of parameter start space:
sobol_size <- 500

# Get list of completed runs:
check_list_fluA <- list.files(path = 'results/', pattern = 'flu_A')
check_list_fluB <- list.files(path = 'results/', pattern = 'flu_B')

# Check for complete results for flu_A/RSV:
if (length(check_list_fluA) != sobol_size * 9) {
  # print(paste0('Missing ', (sobol_size * 5) - length(check_list_fluA), ' record(s).'))
  
  yr_list <- unique(str_sub(check_list_fluA, 15, 18))
  for (yr in yr_list) {
    temp_list <- check_list_fluA[str_detect(check_list_fluA, pattern = yr)]
    
    if (length(temp_list) != sobol_size) {
      # print(paste0('Missing records in ', yr, '!'))
      
      completed_runs <- str_split(temp_list, '_') %>%
        lapply(., function(ix) {ix[6]}) %>%
        unlist() %>%
        str_split(fixed('.')) %>%
        lapply(., function(ix) {ix[1]}) %>%
        unlist() %>%
        as.numeric()
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_A/RSV in ', yr, ', missing run: ', run))
      }
      
    }
  }
}

# Check for complete results for flu_B/RSV:
if (length(check_list_fluB) != sobol_size * 5) {
  # print(paste0('Missing ', (sobol_size * 5) - length(check_list_fluB), ' record(s).'))
  
  yr_list <- unique(str_sub(check_list_fluB, 15, 18))
  for (yr in yr_list) {
    temp_list <- check_list_fluB[str_detect(check_list_fluB, pattern = yr)]
    
    if (length(temp_list) != sobol_size) {
      # print(paste0('Missing records in ', yr, '!'))
      
      completed_runs <- str_split(temp_list, '_') %>%
        lapply(., function(ix) {ix[6]}) %>%
        unlist() %>%
        str_split(fixed('.')) %>%
        lapply(., function(ix) {ix[1]}) %>%
        unlist() %>%
        as.numeric()
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_B/RSV in ', yr, ', missing run: ', run))
      }
      
    }
  }
}

# Clean up:
rm(list = ls())
print('Done.')

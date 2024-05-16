# ---------------------------------------------------------------------------------------------------------------------
# Code to determine whether any results are missing after running cluster code
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Get cluster environmental variables:
fit_canada <- as.logical(Sys.getenv("FITCANADA")); print(fit_canada)
fit_germany <- as.logical(Sys.getenv("FITGERMANY")); print(fit_germany)

# Set size of parameter start space:
sobol_size <- 500

# Get list of completed runs for each virus:
res_files <- list.files(path = 'results', pattern = 'res_', full.names = TRUE)

if (fit_canada | fit_germany) {
  check_list <- list.files(path = 'results/', pattern = 'flu_r')
} else {
  check_list <- list.files(path = 'results/', pattern = 'flu_h1_plus_b_r')
}

# Get vector of seasons for which fitting was done:
if (fit_canada | fit_germany) {
  all_yrs <- unique(str_sub(check_list, 13, 18))
} else {
  all_yrs <- unique(str_sub(check_list, 23, 28))
}
print(all_yrs)
print(length(all_yrs))

# Check for complete results:
if (length(check_list) != sobol_size * length(all_yrs) & length(check_list) > 0) {
  
  for (yr in all_yrs) {
    
    temp_list <- check_list[str_detect(check_list, pattern = yr)]
    
    if (length(temp_list) != sobol_size) {
      
      if (fit_canada | fit_germany) {
        
        completed_runs <- str_split(temp_list, '_') %>%
          lapply(., function(ix) {ix[5]}) %>%
          unlist() %>%
          str_split(fixed('.')) %>%
          lapply(., function(ix) {ix[1]}) %>%
          unlist() %>%
          as.numeric()
        
      } else {
        
        completed_runs <- str_split(temp_list, '_') %>%
          lapply(., function(ix) {ix[8]}) %>%
          unlist() %>%
          str_split(fixed('.')) %>%
          lapply(., function(ix) {ix[1]}) %>%
          unlist() %>%
          as.numeric()
        
      }
      
      missing_runs <- which(!(c(1:sobol_size) %in% completed_runs))
      
      for (run in missing_runs) {
        print(paste0('For flu_h1_plus_b/RSV in ', yr, ', missing run: ', run))
      }
      
    }
    
  }
}

# Clean up:
rm(list = ls())
print('Done.')

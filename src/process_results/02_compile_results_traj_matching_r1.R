# ---------------------------------------------------------------------------------------------------------------------
# Code to compile results from trajectory matching (round 1)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Set estimated parameter names:
# estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2', 'theta_lambda1', 'delta')

# Set parameter values:
vir2 <- 'rsv'
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

# ---------------------------------------------------------------------------------------------------------------------

# Compile results

# Get list of results files:
res_files <- list.files(path = 'results', pattern = 'res_', full.names = TRUE)

# Get virus/season for each result:
which_flu <- str_split(res_files, pattern = '_') %>%
  lapply(function(ix) {
    ix[3]
  }) %>%
  unlist()
which_flu <- paste('flu', which_flu, sep = '_')
yrs <- str_split(res_files, pattern = '_') %>%
  lapply(function(ix) {
    ix[5]
  }) %>%
  unlist()

# Read in all results:
res_full = list()

for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  mutate(virus1 = which_flu,
         year = yrs) %>%
  select(virus1:year, Ri1:loglik)

expect_true(nrow(pars_df) == length(res_files))
expect_true(ncol(pars_df) == (length(estpars) + 3))
expect_true(all(is.finite(pars_df$loglik)))

# Write results to file:
write_rds(res_full, 'traj_match_round1_COMPILED.rds')
write_csv(pars_df, 'res_traj_match_round1.csv')

# ---------------------------------------------------------------------------------------------------------------------

# Sort and store results by virus/season

# Get total number of virus/season pairs:
virus_seasons <- unique(paste(which_flu, yrs, sep = '_'))

# Order unique years correctly:
yrs_unique <- unique(yrs)[order(str_sub(unique(yrs), 2, 3))]
print(yrs_unique)

# Create lists to store results:
res_list_full = res_list = mle_list = slice_list = vector('list', length(virus_seasons))

# Loop through flus/seasons:
counter <- 1
for (vir1 in unique(which_flu)) {
  for (yr in yrs_unique) {
    print(vir1)
    print(yr)
    
    # Check that results exist for combination:
    if (paste(vir1, yr, sep = '_') %in% virus_seasons) {
      
      # Get results for just that flu/season:
      pars_temp <- pars_df %>%
        filter(virus1 == vir1,
               year == yr)
      print(table(pars_temp$virus1))
      print(table(pars_temp$year))
      
      # Sort results and store:
      pars_temp <- pars_temp %>%
        arrange(desc(loglik))
      res_list_full[[counter]] <- pars_temp
      
      # Get only best estimates:
      no_best <- nrow(subset(pars_temp, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = length(estpars))))
      no_best <- max(no_best, 50) # get top 50 if less than 50
      print(no_best)
      
      pars_temp <- pars_temp[1:no_best, ]
      
      # Store results:
      res_list[[counter]] <- pars_temp
      
      # ---------------------------------------------------------------------------------------------------------------
      
      # Get MLE
      
      # Load pomp object:
      source('src/resp_interaction_model.R')
      
      # Create objective function for call to nloptr:
      obj_fun <- traj_objfun(data = resp_mod,
                             est = estpars,
                             partrans = resp_mod@partrans,
                             verbose = TRUE)
      
      # Get MLE and store:
      mle <- setNames(object = as.numeric(pars_temp[1, estpars]),
                      nm = estpars)
      mle_list[[counter]] <- mle
      
      # ---------------------------------------------------------------------------------------------------------------
      
      # Calculate slice likelihoods
      
      # Take slices:
      # slices <- slice_design(center = mle,
      #                        Ri1 = seq(from = 0.9 * mle['Ri1'], to = 1.1 * mle['Ri1'], length.out = 20),
      #                        Ri2 = seq(from = 0.9 * mle['Ri2'], to = 1.1 * mle['Ri2'], length.out = 20),
      #                        I10 = seq(from = 0.9 * mle['I10'], to = 1.1 * mle['I10'], length.out = 20),
      #                        I20 = seq(from = 0.9 * mle['I20'], to = 1.1 * mle['I20'], length.out = 20),
      #                        R10 = seq(from = 0.9 * mle['R10'], to = 1.1 * mle['R10'], length.out = 20),
      #                        R20 = seq(from = 0.9 * mle['R20'], to = 1.1 * mle['R20'], length.out = 20),
      #                        R120 = seq(from = 0.9 * mle['R120'], to = 1.1 * mle['R120'], length.out = 20),
      #                        rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
      #                        rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20)) %>%
      #   mutate(ll = NA)
      slices <- slice_design(center = mle,
                             Ri1 = seq(from = 0.9 * mle['Ri1'], to = 1.1 * mle['Ri1'], length.out = 20),
                             Ri2 = seq(from = 0.9 * mle['Ri2'], to = 1.1 * mle['Ri2'], length.out = 20),
                             I10 = seq(from = 0.9 * mle['I10'], to = 1.1 * mle['I10'], length.out = 20),
                             I20 = seq(from = 0.9 * mle['I20'], to = 1.1 * mle['I20'], length.out = 20),
                             R10 = seq(from = 0.9 * mle['R10'], to = 1.1 * mle['R10'], length.out = 20),
                             R20 = seq(from = 0.9 * mle['R20'], to = 1.1 * mle['R20'], length.out = 20),
                             R120 = seq(from = 0.9 * mle['R120'], to = 1.1 * mle['R120'], length.out = 20),
                             rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                             rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20),
                             theta_lambda1 = seq(from = 0.9 * mle['theta_lambda1'], to = 1.1 * mle['theta_lambda1'], length.out = 20),
                             delta = seq(from = 0.9 * mle['delta'], to = 1.1 * mle['delta'], length.out = 20)) %>%
        mutate(ll = NA)
      
      # Calculate log likelihoods:
      for (i in 1:nrow(slices)) {
        coef(resp_mod, estpars) <- unname(slices[i, estpars])
        x0_trans <- coef(resp_mod, estpars, transform = TRUE)
        slices$ll[i] <- -obj_fun(par = x0_trans)
      }
      rm(i, x0_trans)
      
      # # Check that any NAs are due to initial conditions > 1.0:
      # init_sums <- slices %>%
      #   filter(is.na(ll)) %>%
      #   mutate(init_sum = I10 + I20 + R10 + R20 + R120) %>%
      #   pull(init_sum)
      # expect_true(all(init_sums > 1))
      
      # Check that any NAs are due to initial conditions or reporting rate >= 1.0:
      nas_in_ll <- slices %>%
        filter(is.na(ll)) %>%
        mutate(init_sum = I10 + I20 + R10 + R20 + R120)
      
      # expect_true(all(nas_in_ll$init_sum > 1.0 | nas_in_ll$rho1 > 1.0 | nas_in_ll$rho2 > 1.0))
      expect_true(all(nas_in_ll$init_sum > 1.0 | nas_in_ll$rho1 > 1.0 | nas_in_ll$rho2 > 1.0 | nas_in_ll$theta_lambda1 > 1.0))
      
      # Remove NAs:
      slices <- slices %>%
        filter(!is.na(ll))
      
      # Store:
      slice_list[[counter]] <- slices
      
      # ---------------------------------------------------------------------------------------------------------------
      
      # Iterate to next list position:
      counter <- counter + 1
      
    }
    
  }
}

# Add names to lists:
names(res_list_full) = names(res_list) = names(mle_list) = names(slice_list) = virus_seasons

# Write results to file:
write_rds(res_list_full, 'traj_match_round1_byvirseas_FULL.rds')
write_rds(res_list, 'traj_match_round1_byvirseas_TOP.rds')
write_rds(mle_list, 'traj_match_round1_byvirseas_MLE.rds')
write_rds(slice_list, 'traj_match_round1_byvirseas_SLICE.rds')

print('Done.')

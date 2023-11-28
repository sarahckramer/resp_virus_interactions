# ---------------------------------------------------------------------------------------------------------------------
# Code to compile results from trajectory matching (round 1)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Get cluster environmental variables:
sens <- as.character(Sys.getenv("SENS")); print(sens)
fit_shared <- as.logical(Sys.getenv("FITSHARED")); print(fit_shared)

# Set estimated parameter names:
if (fit_shared) {
  estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2',
               'theta_lambda1', 'theta_lambda2', 'delta1', 'd2', 'alpha', 'phi',
               'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
} else {
  
  if (sens == 'no_rsv_immune') {
    estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'rho1', 'rho2')
  } else {
    estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')
  }
  
}

# Set parameter values:
vir1 <- 'flu_h1_plus_b'
vir2 <- 'rsv'
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
d2_max <- 10.0

# # Fit for synthetic data from age-structured model?:
# age_structured <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Compile results

# Get list of results files:
res_files <- list.files(path = 'results', pattern = 'res_', full.names = TRUE)

# Get season for each result:
yrs <- str_split(res_files, pattern = '_') %>%
  map(~ .x[7]) %>%
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
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  mutate(virus1 = vir1,
         year = yrs) %>%
  select(virus1:year, Ri1:message)

expect_true(nrow(pars_df) == length(res_files))
expect_true(ncol(pars_df) == (length(estpars) + 4))
expect_true(all(is.finite(pars_df$loglik)))

# Write results to file:
write_rds(res_full, 'traj_match_round1_COMPILED.rds')
write_csv(pars_df, 'res_traj_match_round1.csv')

# ---------------------------------------------------------------------------------------------------------------------

# Sort and store results by virus/season

# Get total number of virus/season pairs:
virus_seasons <- unique(paste(vir1, yrs, sep = '_'))

# Order unique years correctly:
yrs_unique <- unique(yrs)[order(str_sub(unique(yrs), 2, 3))]
print(yrs_unique)

# Create lists to store results:
res_list_full = res_list = mle_list = slice_list = vector('list', length(virus_seasons))

# Loop through flus/seasons:
counter <- 1
for (yr in yrs_unique) {
  print(yr)
  
  # Check that results exist for combination:
  if (paste(vir1, yr, sep = '_') %in% virus_seasons) {
    
    # Get results for just that flu/season:
    pars_temp <- pars_df %>%
      filter(virus1 == vir1,
             year == yr)
    print(table(pars_temp$virus1))
    print(table(pars_temp$year))
    
    # Remove fits that didn't converge:
    pars_temp <- pars_temp %>%
      filter(!str_detect(message, 'maxtime')) %>%
      select(-message)
    
    # Check that no model states go below zero:
    source('src/resp_interaction_model.R')
    
    p_mat <- parmat(coef(resp_mod), nrep = nrow(pars_temp))
    for (param in estpars) {
      p_mat[param, ] <- pars_temp %>% pull(param)
    }
    
    rows_to_remove <- trajectory(resp_mod,
                                 params = p_mat,
                                 format = 'data.frame') %>%
      select(!(H1_tot:H2_tot)) %>%
      pivot_longer(X_SS:H2,
                   names_to = 'state') %>%
      filter(value < 0) %>%
      pull(.id) %>%
      unique() %>%
      as.integer()
    print(length(rows_to_remove))
    
    # Remove parameter sets that go negative:
    nrow_check <- nrow(pars_temp) - length(rows_to_remove)
    if (length(rows_to_remove) > 0) {
      pars_temp <- pars_temp[-rows_to_remove, ]
    }
    expect_true(nrow(pars_temp) == nrow_check)
    
    # Sort results and store:
    pars_temp <- pars_temp %>%
      arrange(desc(loglik))
    res_list_full[[counter]] <- pars_temp
    
    # Get only best estimates:
    no_best <- nrow(subset(pars_temp, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = length(estpars))))
    # no_best <- max(no_best, 50) # get top 50 if less than 50
    if (no_best < 50) {
      print('Fewer than 50 top results!')
    }
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
    if (fit_shared) {
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
                             theta_lambda2 = seq(from = 0.9 * mle['theta_lambda2'], to = 1.1 * mle['theta_lambda2'], length.out = 20),
                             delta1 = seq(from = 0.9 * mle['delta1'], to = 1.1 * mle['delta1'], length.out = 20),
                             d2 = seq(from = 0.9 * mle['d2'], to = 1.1 * mle['d2'], length.out = 20),
                             alpha = seq(from = 0.9 * mle['alpha'], to = 1.1 * mle['alpha'], length.out = 20),
                             phi = seq(from = 0.9 * mle['phi'], to = 1.1 * mle['phi'], length.out = 20),
                             eta_temp1 = seq(from = 0.9 * mle['eta_temp1'], to = 1.1 * mle['eta_temp1'], length.out = 20),
                             eta_temp2 = seq(from = 0.9 * mle['eta_temp2'], to = 1.1 * mle['eta_temp2'], length.out = 20),
                             eta_ah1 = seq(from = 0.9 * mle['eta_ah1'], to = 1.1 * mle['eta_ah1'], length.out = 20),
                             eta_ah2 = seq(from = 0.9 * mle['eta_ah2'], to = 1.1 * mle['eta_ah2'], length.out = 20)) %>%
        mutate(ll = NA)
    } else {
      
      if (sens == 'no_rsv_immune') {
        slices <- slice_design(center = mle,
                               Ri1 = seq(from = 0.9 * mle['Ri1'], to = 1.1 * mle['Ri1'], length.out = 20),
                               Ri2 = seq(from = 0.9 * mle['Ri2'], to = 1.1 * mle['Ri2'], length.out = 20),
                               I10 = seq(from = 0.9 * mle['I10'], to = 1.1 * mle['I10'], length.out = 20),
                               I20 = seq(from = 0.9 * mle['I20'], to = 1.1 * mle['I20'], length.out = 20),
                               R10 = seq(from = 0.9 * mle['R10'], to = 1.1 * mle['R10'], length.out = 20),
                               rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                               rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20)) %>%
          mutate(ll = NA)
      } else {
        slices <- slice_design(center = mle,
                               Ri1 = seq(from = 0.9 * mle['Ri1'], to = 1.1 * mle['Ri1'], length.out = 20),
                               Ri2 = seq(from = 0.9 * mle['Ri2'], to = 1.1 * mle['Ri2'], length.out = 20),
                               I10 = seq(from = 0.9 * mle['I10'], to = 1.1 * mle['I10'], length.out = 20),
                               I20 = seq(from = 0.9 * mle['I20'], to = 1.1 * mle['I20'], length.out = 20),
                               R10 = seq(from = 0.9 * mle['R10'], to = 1.1 * mle['R10'], length.out = 20),
                               R20 = seq(from = 0.9 * mle['R20'], to = 1.1 * mle['R20'], length.out = 20),
                               R120 = seq(from = 0.9 * mle['R120'], to = 1.1 * mle['R120'], length.out = 20),
                               rho1 = seq(from = 0.9 * mle['rho1'], to = 1.1 * mle['rho1'], length.out = 20),
                               rho2 = seq(from = 0.9 * mle['rho2'], to = 1.1 * mle['rho2'], length.out = 20)) %>%
          mutate(ll = NA)
      }
      
    }
    
    # Calculate log likelihoods:
    for (i in 1:nrow(slices)) {
      coef(resp_mod, estpars) <- unname(slices[i, estpars])
      x0_trans <- coef(resp_mod, estpars, transform = TRUE)
      slices$ll[i] <- -obj_fun(par = x0_trans)
    }
    rm(i, x0_trans)
    
    # Check that any NAs are due to initial conditions or other unrealistic parameter values:
    if (sens == 'no_rsv_immune') {
      nas_in_ll <- slices %>%
        filter(is.na(ll)) %>%
        mutate(init_sum = I10 + I20 + R10)
    } else {
      nas_in_ll <- slices %>%
        filter(is.na(ll)) %>%
        mutate(init_sum = I10 + I20 + R10 + R20 + R120)
    }
    
    if (fit_shared) {
      expect_true(all(nas_in_ll$init_sum > 1.0 | nas_in_ll$rho1 > 1.0 | nas_in_ll$rho2 > 1.0 |
                        nas_in_ll$theta_lambda1 > 1.0 | nas_in_ll$theta_lambda2 > 1.0))
    } else {
      expect_true(all(nas_in_ll$init_sum > 1.0 | nas_in_ll$rho1 > 1.0 | nas_in_ll$rho2 > 1.0))
    }
    
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

# Add names to lists:
names(res_list_full) = names(res_list) = names(mle_list) = names(slice_list) = virus_seasons

# Write results to file:
write_rds(res_list_full, 'traj_match_round1_byvirseas_FULL.rds')
write_rds(res_list, 'traj_match_round1_byvirseas_TOP.rds')
write_rds(mle_list, 'traj_match_round1_byvirseas_MLE.rds')
write_rds(slice_list, 'traj_match_round1_byvirseas_SLICE.rds')

print('Done.')

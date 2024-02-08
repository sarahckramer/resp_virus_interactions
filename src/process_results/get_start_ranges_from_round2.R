# ---------------------------------------------------------------------------------------------------------------------
# Get start ranges for from round2 trajectory matching results (for profile likelihoods, bootstraps, or to rerun round2)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set directory where results from round2 fits are stored:
res_dir <- 'results/round2_fit/sens/hk_plus_canada/round2_1_cont2/'

# Which round of fits?:
which_round <- str_split(res_dir, '_')[[1]][which(!is.na(as.numeric(str_split(res_dir, '_')[[1]])))]

# Are results from a sensitivity analysis?:
if (str_detect(res_dir, 'sinusoidal')) {
  sens <- 'sinusoidal_forcing'
} else if (str_detect(res_dir, 'h3_covar')) {
  sens <- 'h3_covar'
} else if (str_detect(res_dir, 'less_circ_h3')) {
  sens <- 'less_circ_h3'
} else if (str_detect(res_dir, 'no_rsv_immune')) {
  sens <- 'no_rsv_immune'
} else if (str_detect(res_dir, 'no_ah')) {
  sens <- 'no_ah'
} else if (str_detect(res_dir, 'rhino_covar')) {
  sens <- 'rhino_covar'
} else if (str_detect(res_dir, 'no_int')) {
  sens <- 'no_int'
} else {
  sens <- 'main'
}

# Check that directory for storing results exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/round2_CIs/')) {
  dir.create('results/round2_CIs/')
}

if (sens == 'main') {
  
  if (str_detect(res_dir, 'age_structured')) {
    
    if(!dir.exists('results/round2_CIs/sens/')) {
      dir.create('results/round2_CIs/sens/')
    }
    if(!dir.exists('results/round2_CIs/sens/age_structured/')) {
      dir.create('results/round2_CIs/sens/age_structured/')
    }
    
    new_dir <- paste0('results/round2_CIs/sens/age_structured/from_2_', which_round, '/')
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
  } else {
    
    if(!dir.exists('results/round2_CIs/sens/')) {
      dir.create('results/round2_CIs/sens/')
    }
    if(!dir.exists('results/round2_CIs/sens/hk_plus_canada/')) {
      dir.create('results/round2_CIs/sens/hk_plus_canada/')
    }
    
    new_dir <- paste0('results/round2_CIs/sens/hk_plus_canada/from_2_', which_round, '/')
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
    
  }
  
} else {
  
  if(!dir.exists('results/round2_CIs/sens/')) {
    dir.create('results/round2_CIs/sens/')
  }
  if(!dir.exists(paste0('results/round2_CIs/sens/', sens, '/'))) {
    dir.create(paste0('results/round2_CIs/sens/', sens, '/'))
  }
  
  new_dir <- paste0('results/round2_CIs/sens/', sens, '/from_2_', which_round, '/')
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }
  
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
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  # expect_equal(df_use, 90)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(table(pars_df$message))
  
  # If only one parameter set is in this range, MLE has not yet been reached; take top 5% of fits instead:
  if (no_best == 1) {
    is_mle <- FALSE
    print('MLE not reached!')
    no_best <- 25
  } else {
    is_mle <- TRUE
  }
  print(no_best)
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  print(summary(pars_top$loglik))
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # If none remaining, print warning:
  if(nrow(pars_top) == 0) {
    print('No convergence among best-fit runs!')
  }
  print(no_best)
  print(summary(pars_top$loglik))
  
  # Store original results before removing unrealistic parameter values:
  pars_top_orig <- pars_top
  
  # # Remove where d2 > 10 and theta_lambda2 != 1.0:
  # pars_top <- pars_top %>%
  #   filter(!(d2 > 10.0 & theta_lambda2 < 0.99))
  
  # Set unrealistic values to NA:
  if (!str_detect(filepath, 'no_int')) {
    pars_top$delta1[pars_top$delta1 > 7.0] <- NA
    pars_top$d2[pars_top$d2 > 10.0] <- NA
  }
  pars_top$hk_rho1[pars_top$hk_rho1 == 1.0] <- NA
  pars_top$hk_rho2[pars_top$hk_rho2 == 1.0] <- NA
  pars_top$can_rho2[pars_top$can_rho2 == 1.0] <- NA
  
  # Since phi=0 is equivalent to phi=52.25, don't use full range; transform so that we can select from only best-supported range:
  par(mfrow = c(2, 2))
  hist(pars_top$hk_phi, breaks = 50)
  hist(pars_top$can_phi, breaks = 50)
  pars_top <- pars_top %>%
    mutate(hk_phi = if_else(hk_phi < 5, hk_phi + 52.25, hk_phi),
           can_phi = if_else(can_phi < 5, can_phi + 52.25, can_phi))
  hist(pars_top$hk_phi, breaks = 50)
  hist(pars_top$can_phi, breaks = 50)
  
  # If using sinusoidal forcing, do the same for phi1 and phi2:
  if ('can_phi1' %in% names(pars_top)) {
    
    par(mfrow = c(4, 2))
    hist(pars_top$hk_phi1, breaks = 50)
    hist(pars_top$hk_phi2, breaks = 50)
    hist(pars_top$can_phi1, breaks = 50)
    hist(pars_top$can_phi2, breaks = 50)
    
    pars_top <- pars_top %>%
      mutate(hk_phi1 = if_else(hk_phi1 < 5, hk_phi1 + 52.25, hk_phi1),
             hk_phi2 = if_else(hk_phi2 < 5, hk_phi2 + 52.25, hk_phi2),
             can_phi1 = if_else(can_phi1 < 5, can_phi1 + 52.25, can_phi1),
             can_phi2 = if_else(can_phi2 < 5, can_phi2 + 52.25, can_phi2))
    
    hist(pars_top$hk_phi1, breaks = 50)
    hist(pars_top$hk_phi2, breaks = 50)
    hist(pars_top$can_phi1, breaks = 50)
    hist(pars_top$can_phi2, breaks = 50)
    
  }
  
  # Return formatted results:
  return(list(pars_top, pars_top_orig, is_mle))
  
}

# Read in results:
res <- load_and_format_mega_results(res_dir)

# Check that best-fit parameter values do not lead trajectories to drop below 0:
res_orig <- res[[2]]

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
if (sens == 'no_rsv_immune') {
  unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')
}

global_estpars <- c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')
shared_estpars <- c('rho1', 'rho2', 'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
true_estpars <- c(global_estpars, shared_estpars, unit_estpars)

sens <- 'sinusoidal_forcing'
source('src/functions/setup_global_likelilhood.R')

traj_list <- lapply(1:length(c(seasons_hk, seasons_can)), function(ix) {
  
  if (ix <= length(seasons_hk)) {
    
    pars_temp <- res_orig %>%
      select(all_of(global_estpars),
             paste0('hk_', shared_estpars),
             paste0(seasons_hk[ix], '_hk_', unit_estpars))
    
  } else {
    
    pars_temp <- res_orig %>%
      select(all_of(global_estpars),
             paste0('can_', shared_estpars),
             paste0(seasons_can[ix - length(seasons_hk)], '_can_', unit_estpars))
    
  }
  
  names(pars_temp) <- true_estpars
  
  p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(pars_temp))
  for (param in names(pars_temp)) {
    p_mat[param, ] <- pars_temp %>% pull(param)
  }

  trajectory(object = po_list[[ix]],
             params = p_mat,
             format = 'data.frame') %>%
    select(!(H1_tot:H2_tot)) %>%
    pivot_longer(X_SS:H2,
                 names_to = 'state')
  
})

expect_false(any(lapply(traj_list, function(ix) {
  any(ix[, 'value'] < 0)
}) %>%
  unlist()))

# Are results the MLE?
is_mle <- res[[3]]

res_dir_comp <- str_split(res_dir, '_')[[1]]
res_dir_comp[which(!is.na(as.numeric(res_dir_comp)))] <- as.character(as.numeric(which_round) - 1)
res_dir_prev <- paste(res_dir_comp, collapse = '_')
rm(res_dir_comp)

is_mle_prev <- try(
  load_and_format_mega_results(res_dir_prev)[[3]]
)

if (inherits(is_mle_prev, 'try-error')) {
  is_mle_prev <- FALSE
}

# Save MLEs:
if (is_mle & is_mle_prev) {
  
  res_orig <- res_orig %>%
    select(-loglik)
  
  if (str_detect(res_dir, 'sens')) {
    
    if (str_detect(res_dir, 'hk_plus_canada')) {
      write_rds(res_orig, file = 'results/round2_fit/sens/hk_plus_canada/MLEs.rds')
    } else {
      write_rds(res_orig, file = paste0(paste(str_split(res_dir, '/')[[1]][1:(length(str_split(res_dir, '/')[[1]]) - 2)], collapse = '/'), '/MLEs.rds'))
    }
    
  } else if (str_detect(res_dir, 'age_structured')) {
    write_rds(res_orig, file = 'results/age_structured_SA/MLEs.rds')
  } else {
    write_rds(res_orig, file = 'results/MLEs.rds')
  }
  
}

# Get results for determining start ranges:
res <- res[[1]] %>%
  select(-loglik)

# Get minimum and maximum start values:
ci_start <- as.data.frame(rbind(summarise(res, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                summarise(res, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))

# Any parameters where minimum and maximum are equal?:
no_range <- c()
for (i in 1:ncol(ci_start)) {
  if (identical(ci_start[1, i], ci_start[2, i])) {
    no_range <- c(no_range, i)
  }
}
if (length(no_range) > 0) {
  print('Warning: Some parameters have same min and max values!')
  print(names(ci_start)[no_range])
}

# Possible that d2 ranges are missing if all top fits were > 10; if so, replace:
if (any(ci_start == Inf)) {
  if (any(ci_start$d2 == Inf)) {
    ci_start$d2 <- c(0, 10)
  } else if (any(ci_start$delta1 == Inf)) {
    ci_start$delta1 <- c(7 / 60, 7)
  } else {
    print('Range is NA for some other parameter.')
  }
}

# Check that sums of initial conditions can't sum to >1:
init_cond_estpars <- c('I10', 'I20', 'R10', 'R20', 'R120')

if (str_detect(res_dir, 'hk_plus_canada')) {
  
  sums <- ci_start %>%
    mutate(minmax = c('min', 'max')) %>%
    select(contains(init_cond_estpars), minmax) %>%
    pivot_longer(-minmax) %>%
    mutate(location = str_remove(str_sub(name, 8, 10), '_'),
           season = str_sub(name, 1, 6)) %>%
    group_by(location, season, minmax) %>%
    summarise(sum = sum(value))
  
} else {
  
  sums <- ci_start %>%
    mutate(minmax = c('min', 'max')) %>%
    select(contains(init_cond_estpars), minmax) %>%
    pivot_longer(-minmax) %>%
    mutate(season = str_sub(name, 1, 6)) %>%
    group_by(season, minmax) %>%
    summarise(sum = sum(value))
  
}

if (any(sums %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}

if (any(sums %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  ids <- sums %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    select(location, season) %>%
    mutate(id = paste(season, location, sep = '_')) %>%
    pull(id)
  
  for (id_temp in ids) {
    
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start[2, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(id_temp))
    red_needed <- sums %>%
      mutate(id = paste(season, location, sep = '_')) %>%
      filter(id == id_temp,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start[1, ] %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(id_temp))
    
    if (!all(new_upper_bounds > orig_lower_bounds)) {
      new_upper_bounds_try <- orig_upper_bounds
      new_upper_bounds_try[-which(orig_lower_bounds >= new_upper_bounds)] <- orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] - (red_needed * (orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] / sum(orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)])))
      new_upper_bounds <- new_upper_bounds_try
      rm(new_upper_bounds_try)
    }
    
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start[2, which(str_detect(names(ci_start), id_temp) &
                        (str_detect(names(ci_start), 'R10') |
                           str_detect(names(ci_start), 'R20') |
                           str_detect(names(ci_start), 'R120')))] <- new_upper_bounds
    expect_lt(ci_start[2, ] %>% select(contains(id_temp)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(id_temp)
}

# Write start ranges to file:
write_rds(ci_start, file = paste0(new_dir, 'round2CI_startvals.rds'))

# Clean up:
rm(list = ls())

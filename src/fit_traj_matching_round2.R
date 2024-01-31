# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# Round 2: Fit all seasons simultaneously, constraining shared parameters to be the same for all seasons
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(nloptr)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
which_round <- as.integer(Sys.getenv("WHICHROUND")); print(which_round)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
sens <- as.character(Sys.getenv("SENS")); print(sens)

# # Set parameters for local run:
# jobid <- 1
# no_jobs <- 1
# sobol_size <- 10
# which_round <- 1
# search_type <- 'round1_CIs'
# sens <- 'sinusoidal_forcing' # 'main', 'less_circ_h3', 'sinusoidal_forcing', 'no_ah', 'no_int', 'no_rsv_immune', 'h3_covar', 'rhino_covar'

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'

seasons_hk <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
seasons_can <- c('s10-11', 's11-12', 's12-13', 's13-14')

time_max <- 23.75 # Maximal execution time (in hours)

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

# ---------------------------------------------------------------------------------------------------------------------

# Functions

create_obj_fxn <- function(po, estpars) {
  # Creates the objective function for a given season
  # param po: A pomp model object for a specific season
  # param estpars: A vector listing the parameters to be fit
  # returns: The negative log likelihood
  
  ofun <- traj_objfun(data = po, 
                      est = estpars, 
                      partrans = po@partrans,
                      verbose=TRUE)
  
  return(ofun)
}

calculate_global_loglik <- function(trans_vals) {
  # Calculates the log-likelihood for each season, and combines to yield a global log-likelihood
  # param trans_vals: Unnamed vector of transformed parameters; fxn only works if this is the only input?
  # returns: The global, negative log-likelihood
  
  # Add names to vector:
  if(is.null(names(trans_vals))) names(trans_vals) <- x0_trans_names
  
  # Split params into global, shared, and unit params:
  shared_hk_in <- as.data.frame(sapply(paste0('^', 'hk_'), grep, x = names(trans_vals)))
  unit_hk_in <- as.data.frame(sapply(paste0('^', seasons_hk, '_hk_'), grep, x = names(trans_vals)))
  names(unit_hk_in) <- seasons_hk
  
  shared_can_in <- as.data.frame(sapply(paste0('^', 'can_'), grep, x = names(trans_vals)))
  unit_can_in <- as.data.frame(sapply(paste0('^', seasons_can, '_can_'), grep, x = names(trans_vals)))
  names(unit_can_in) <- seasons_can
  
  global_params <- trans_vals[-c(unlist(shared_hk_in), unlist(shared_can_in), unique(unlist(unit_hk_in)), unique(unlist(unit_can_in)))]
  
  shared_hk_params <- trans_vals[unlist(shared_hk_in)]
  names(shared_hk_params) <- gsub('hk_', '', names(shared_hk_params))
  shared_hk_params <- c(global_params, shared_hk_params)
  
  shared_can_params <- trans_vals[unlist(shared_can_in)]
  names(shared_can_params) <- gsub('can_', '', names(shared_can_params))
  shared_can_params <- c(global_params, shared_can_params)
  
  unit_hk_params <- list()
  for (i in 1:length(seasons_hk)) {
    unit <- trans_vals[unit_hk_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_hk_params[[i]] <- unit
    }
  }
  
  unit_can_params <- list()
  for (i in 1:length(seasons_can)) {
    unit <- trans_vals[unit_can_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_can_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season (HK):
  units_hk_ll <- rep(NA, length(seasons_hk))
  
  for (i in 1:length(seasons_hk)) {
    if(length(unit_hk_params) > 0) {
      params_temp <- c(shared_hk_params, unit_hk_params[[i]])
    } else {
      params_temp <- c(shared_hk_params)
    }
    
    units_hk_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Get -ll for each season (Canada):
  units_can_ll <- rep(NA, length(seasons_can))
  
  for (i in 1:length(seasons_can)) {
    if (length(unit_can_params) > 0) {
      params_temp <- c(shared_can_params, unit_can_params[[i]])
    } else {
      params_temp <- c(shared_can_params)
    }
    
    units_can_ll[i] <- obj_fun_list[[i + length(seasons_hk)]](params_temp)
  }
  
  # Calculate -ll for each location:
  hk_ll <- sum(units_hk_ll)
  can_ll <- sum(units_can_ll)
  
  # Calculate global -ll:
  glob_ll <- hk_ll + can_ll
  return(glob_ll)
}

transform_params <- function(orig_vals, po1, po2, seas1, seas2, params_all, params_global, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po1: pomp object with correct parameter transformations (Hong Kong)
  # param po2: pomp object with correct parameter transformations (Canada)
  # param seas1: Vector containing all seasons of interest (Hong Kong)
  # params seas2: Vector containing all seasons of interest (Canada)
  # param params_all: Names of all parameters to be estimated
  # param params_global: Names of all global parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_global <- orig_vals[which(params_global %in% params_all)]
  
  coef(po1, params_global) <- orig_vals_global
  trans_vals_global <- coef(po1, params_global, transform = TRUE)
  trans_vals <- trans_vals_global
  
  coef(po2, params_global) <- orig_vals_global
  
  for (i in c('hk', 'can')) {
    
    orig_vals_shared <- orig_vals[grep(paste0('^', i, '_'), params_all, value = TRUE)]
    names(orig_vals_shared) <- gsub(paste0(i, '_'), '', names(orig_vals_shared))
    
    if (i == 'hk') {
      
      coef(po1, names(orig_vals_shared)) <- orig_vals_shared
      trans_vals_shared <- coef(po1, names(orig_vals_shared), transform = TRUE)
      
    } else if (i == 'can') {
      
      coef(po2, names(orig_vals_shared)) <- orig_vals_shared
      trans_vals_shared <- coef(po2, names(orig_vals_shared), transform = TRUE)
      
    }
    
    names(trans_vals_shared) <- paste0(i, '_', names(orig_vals_shared))
    trans_vals <- c(trans_vals, trans_vals_shared)
    
  }
  
  po_save1 <- po1
  po_save2 <- po2
  
  for (i in 1:length(c(seas1, seas2))) {
    
    po1 <- po_save1
    po2 <- po_save2
    
    if (i <= length(seas1)) {
      
      orig_vals_unit <- orig_vals[grep(paste0('^', seas1[i], '_hk'), params_all, value = TRUE)]
      names(orig_vals_unit) <- gsub(paste0(seas1[i], '_hk_'), '', names(orig_vals_unit))
      
      coef(po1, names(orig_vals_unit)) <- orig_vals_unit
      trans_vals_unit <- coef(po1, names(orig_vals_unit), transform = TRUE)
      names(trans_vals_unit) <- paste0(seas1[i], '_hk_', names(orig_vals_unit))
      
    } else {
      
      orig_vals_unit <- orig_vals[grep(paste0('^', seas2[i - length(seas1)], '_can'), params_all, value = TRUE)]
      names(orig_vals_unit) <- gsub(paste0(seas2[i - length(seas1)], '_can_'), '', names(orig_vals_unit))
      
      coef(po2, names(orig_vals_unit)) <- orig_vals_unit
      trans_vals_unit <- coef(po2, names(orig_vals_unit), transform = TRUE)
      names(trans_vals_unit) <- paste0(seas2[i - length(seas1)], '_can_', names(orig_vals_unit))
      
    }
    
    trans_vals <- c(trans_vals, trans_vals_unit)
    
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po1, po2, seas1, seas2, params_all, params_global, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po1: pomp object with correct parameter transformations (Hong Kong)
  # param po2: pomp object with correct parameter transformations (Canada)
  # param seas1: Vector containing all seasons of interest (Hong Kong)
  # params seas2: Vector containing all seasons of interest (Canada)
  # param params_all: Names of all parameters to be estimated
  # param params_global: Names of all global parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_global <- trans_vals[which(params_global %in% params_all)]
  
  coef(po1, params_global, transform = TRUE) <- trans_vals_global
  orig_vals_global <- coef(po1, params_global)
  orig_vals <- orig_vals_global
  
  coef(po2, params_global, transform = TRUE) <- trans_vals_global
  
  for (i in c('hk', 'can')) {
    
    trans_vals_shared <- trans_vals[grep(paste0('^', i, '_'), params_all, value = TRUE)]
    names(trans_vals_shared) <- gsub(paste0(i, '_'), '', names(trans_vals_shared))
    
    if (i == 'hk') {
      
      coef(po1, names(trans_vals_shared), transform = TRUE) <- trans_vals_shared
      orig_vals_shared <- coef(po1, params_shared)
      
    } else if (i == 'can') {
      
      coef(po2, names(trans_vals_shared), transform = TRUE) <- trans_vals_shared
      orig_vals_shared <- coef(po2, params_shared)
      
    }
    
    names(orig_vals_shared) <- paste0(i, '_', names(trans_vals_shared))
    orig_vals <- c(orig_vals, orig_vals_shared)
    
  }
  
  po_save1 <- po1
  po_save2 <- po2
  
  for (i in 1:length(c(seas1, seas2))) {
    
    po1 <- po_save1
    po2 <- po_save2
    
    if (i <= length(seas1)) {
      
      trans_vals_unit <- trans_vals[grep(paste0('^', seas1[i], '_hk'), params_all, value = TRUE)]
      names(trans_vals_unit) <- gsub(paste0(seas1[i], '_hk_'), '', names(trans_vals_unit))
      
      coef(po1, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
      orig_vals_unit <- coef(po1, names(trans_vals_unit))
      
      names(orig_vals_unit) <- paste0(seas1[i], '_hk_', names(trans_vals_unit))
      
    } else {
      
      trans_vals_unit <- trans_vals[grep(paste0('^', seas2[i - length(seas1)], '_can'), params_all, value = TRUE)]
      names(trans_vals_unit) <- gsub(paste0(seas2[i - length(seas1)], '_can_'), '', names(trans_vals_unit))
      
      coef(po2, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
      orig_vals_unit <- coef(po2, names(trans_vals_unit))
      
      names(orig_vals_unit) <- paste0(seas2[i - length(seas1)], '_can_', names(trans_vals_unit))
      
    }
    
    orig_vals <- c(orig_vals, orig_vals_unit)
    
  }
  
  return(orig_vals)
}

# ---------------------------------------------------------------------------------------------------------------------

# Fit using trajectory matching

# Loop through years and construct pomp models:
po_list <- vector('list', length(c(seasons_hk, seasons_can)))
for (yr_index in 1:length(po_list)) {
  
  if (yr_index <= length(seasons_hk)) {
    
    yr <- seasons_hk[yr_index]
    fit_canada <- FALSE
    vir1 <- 'flu_h1_plus_b'
    
  } else {
    
    yr <- seasons_can[yr_index - length(seasons_hk)]
    fit_canada <- TRUE
    vir1 <- 'flu'
    
  }
  
  print(fit_canada)
  print(yr)
  
  # Load data and create pomp object:
  source('src/resp_interaction_model.R')
  
  # Check whether any data for given season:
  if (exists('resp_mod')) {
    
    # Add pomp object to list:
    po_list[[yr_index]] <- resp_mod
    
  }
  
  # Remove pomp object before repeating loop:
  rm(resp_mod)
  
}

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# ---------------------------------------------------------------------------------------------------------------------

# Choose parameters to estimate:
global_estpars <- c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')
shared_estpars <- c('rho1', 'rho2', 'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

shared_sp_estpars <- c()
for (i in c('hk', 'can')) {
  shared_sp_estpars <- c(shared_sp_estpars, paste(i, shared_estpars, sep = '_'))
}

unit_sp_estpars <- c()
for (i in 1:length(seasons_hk)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons_hk[i], 'hk', unit_estpars, sep = '_'))
}
for (i in 1:length(seasons_can)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons_can[i], 'can', unit_estpars, sep = '_'))
}
rm(i)

true_estpars <- c(global_estpars, shared_estpars, unit_estpars)
estpars <- c(global_estpars, shared_sp_estpars, unit_sp_estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Set upper/lower values for global params:
start_range <- data.frame(rho1 = c(0, 1.0),
                          rho2 = c(0, 1.0),
                          theta_lambda1 = c(0, 1.0),
                          theta_lambda2 = c(0, 1.0),
                          theta_rho1 = c(0, 5.0),
                          theta_rho2 = c(0, 5.0),
                          delta1 = c(7 / 60, 7),
                          d2 = c(0, 10),
                          alpha = c(0, 0.5),
                          phi = c(0, 52.25),
                          eta_temp1 = c(-0.5, 0.5),
                          eta_temp2 = c(-0.5, 0.5),
                          eta_ah1 = c(-0.5, 0.5),
                          eta_ah2 = c(-0.5, 0.5),
                          b1 = c(0.05, 0.2),
                          b2 = c(0.05, 0.2),
                          phi1 = c(0, 52.25),
                          phi2 = c(0, 52.25),
                          beta_h3 = c(0, 5.0),
                          beta_rhino = c(0, 5.0))
start_range <- start_range[c(global_estpars, shared_estpars)]
start_range <- start_range %>%
  rename_with(~ paste('hk', ., sep = '_'), .cols = all_of(shared_estpars)) %>%
  bind_cols(start_range[shared_estpars]) %>%
  rename_with(~ paste('can', ., sep = '_'), .cols = all_of(shared_estpars))

# ---------------------------------------------------------------------------------------------------------------------

# Get 95% CI from round 1 for unit params:
tj_res_list_hk <- read_rds('results/round1_fitsharedFALSE/traj_match_round1_byvirseas_TOP.rds')
tj_res_list_can <- read_rds('results/round2_fit/sens/canada/round1_fitsharedFALSE/traj_match_round1_byvirseas_TOP.rds')

ci_list <- vector('list', length(c(seasons_hk, seasons_can)))

for (i in 1:length(ci_list)) {
  
  if (i <= length(seasons_hk)) {
    
    yr <- seasons_hk[i]
    tj_res_temp <- tj_res_list_hk[[which(str_detect(names(tj_res_list_hk), yr))]] %>%
      select(all_of(unit_estpars)) %>% na.omit()
    
  } else {
    
    yr <- seasons_can[i - length(seasons_hk)]
    tj_res_temp <- tj_res_list_can[[which(str_detect(names(tj_res_list_can), yr))]] %>%
      select(all_of(unit_estpars)) %>% na.omit()
    
  }
  
  ci_temp <- as.data.frame(rbind(summarise(tj_res_temp, across(.cols = everything(), min)),
                                 summarise(tj_res_temp, across(.cols = everything(), max))))
  
  # Check that initial conditions can't sum to >1:
  sums <- ci_temp %>% select(I10:R120) %>% rowSums()
  
  if (sums[1] > 1.0) {
    print('Lower bounds sum to more than 1!')
    stop()
  }
  
  if (sums[2] > 1.0) {
    
    # Reduce the upper bounds of R10/R20/R120 proportionally:
    orig_upper_bounds <- ci_temp[2, ] %>% select(R10:R120)
    red_needed <- sums[2] - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_temp[1, ] %>% select(R10:R120)
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Ensure upper bounds now sum to 1 or less:
    ci_temp[2, c('R10', 'R20', 'R120')] <- new_upper_bounds
    expect_lt(ci_temp[2, ] %>% select(I10:R120) %>% sum(), 1.0)
    
  }
  
  # Store in list:
  ci_list[[i]] <- ci_temp
}
rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Get data frame of all ranges:
if (search_type == 'round2_CIs') {
  
  start_range <- read_rds(paste0('results/round2_CIs/from_2_', which_round - 1, '/round2CI_startvals_H1_plus_B.rds'))
  
} else {
  
  for (i in 1:length(ci_list)) {
    
    start_range_temp <- ci_list[[i]]
    
    if (i <= length(seasons_hk)) {
      names(start_range_temp) <- paste(seasons_hk[i], 'hk', unit_estpars, sep = '_')
    } else {
      names(start_range_temp) <- paste(seasons_can[i - length(seasons_hk)], 'can', unit_estpars, sep = '_')
    }
    
    start_range <- start_range %>%
      bind_cols(start_range_temp)
    
    rm(start_range_temp)
    
  }
  rm(i)
  
  start_range <- start_range[, estpars]
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Get starting values for each parameter:
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

if (search_type == 'round2_CIs') {
  start_values <- start_values %>%
    mutate(phi = if_else(phi > 52.25, phi - 52.25, phi))
  
  if ('phi1' %in% names(start_values)) {
    start_values <- start_values %>%
      mutate(phi1 = if_else(phi1 > 52.25, phi1 - 52.25, phi1),
             phi2 = if_else(phi2 > 52.25, phi2 - 52.25, phi2))
  }
}

# Check that starting values and estpars are correct:
print(start_range)
print(summary(start_values))
print(estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# Set maximal execution time for each estimation:
nmins_exec <- time_max * 60 / (sobol_size / no_jobs)
print(sprintf("Max estimation time=%.1f min", nmins_exec))

# Get unique identifiers:
sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)

# Fit:
for (i in seq_along(sub_start)) {
  
  print(paste0('Estimation: ', sub_start[i]))
  
  # Get start values:
  x0 <- as.numeric(start_values[sub_start[i], ])
  x0_trans <- transform_params(x0, po_list[[1]], po_list[[7]], seasons_hk, seasons_can, estpars, global_estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)
  
  # Check that parameter transformations correct:
  x0_orig <- back_transform_params(x0_trans, po_list[[1]], po_list[[7]], seasons_hk, seasons_can, estpars, global_estpars, shared_estpars)
  expect_equal(x0, unname(x0_orig))
  rm(x0_orig)
  
  # Calculate initial log-likelihood:
  print(-1 * calculate_global_loglik(x0_trans))
  
  # Fit models:
  tic <- Sys.time()
  m <- try(
    nloptr(x0 = x0_trans, 
           eval_f = calculate_global_loglik,
           opts = list(algorithm = "NLOPT_LN_SBPLX",
                       maxtime = 60 * nmins_exec,
                       maxeval = -1, # Negative value: criterion is disabled
                       xtol_rel = -1, # Default value: 1e-4
                       print_level = 0))
  )
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'hours'
  print(etime)
  
  # If estimation is successful, save results:
  if (!inherits(m, 'try-error')) {
    x0_fit <- m$solution
    names(x0_fit) <- x0_trans_names
    x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], po_list[[7]], seasons_hk, seasons_can, estpars, global_estpars, shared_estpars)
    
    out <- list(estpars = x0_fit_untrans,
                ll = -m$objective,
                conv = m$status,
                message = m$message,
                niter = m$iterations,
                etime = as.numeric(etime))
    
    # Write to file:
    saveRDS(out, file = sprintf('results/res_%d_%d.rds',
                                jobid,
                                sub_start[i])
    )
    
    # Print results:
    print(out$ll)
    print(out$estpars, digits = 2)
    print(out$conv)
    print(out$message)
  }
  
}
rm(i)

# Clean up:
rm(list = ls())

print('Done!')

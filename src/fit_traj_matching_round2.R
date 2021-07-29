# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(nloptr)

# Set parameters for run:
jobid <- 1
no_jobs <- 20
time_max <- 11.5 # Maximal execution time (in hours)

debug_bool <- FALSE
# yr <- 2006 # 2004:2014 # 2006 = 2005-06 season
vir1 <- 'flu_B' # 'flu_A', 'flu_B'
vir2 <- 'rsv'

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

seasons <- 2006:2014

sobol_size <- 100
search_type <- 'round1_CIs'
int_eff <- 'susc' # 'susc' or 'sev' - fit impact of interaction on susceptibility or severity?

# CHECK: Fit 'panel' but with no interaction?

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
  
  # Split params into shared and unit params:
  unit_in <- as.data.frame(sapply(paste0('^', seasons, '_'), grep, x = names(trans_vals)))
  if (ncol(unit_in) == 1) {
    unit_in <- as.data.frame(matrix(data = c(unlist(unit_in)), nrow = 1))
  }
  names(unit_in) <- seasons
  shared_params <- trans_vals[-unique(unlist(unit_in))]
  
  unit_params <- list()
  for (i in 1:length(seasons)) {
    unit <- trans_vals[unit_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 2]
      unit_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season:
  units_ll <- rep(NA, length(seasons))
  
  for (i in 1:length(seasons)) {
    if (length(unit_params) > 0) {
      params_temp <- c(shared_params, unit_params[[i]])
    } else {
      params_temp <- c(shared_params)
    }
    
    units_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Calculate global -ll:
  glob_ll <- sum(units_ll)
  return(glob_ll)
}

transform_params <- function(orig_vals, po, seas, params_all, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_shared <- orig_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared) <- orig_vals_shared
  trans_vals_shared <- coef(po, params_shared, transform = TRUE)
  trans_vals <- trans_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    orig_vals_unit <- orig_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(orig_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(orig_vals_unit))
    
    coef(po, names(orig_vals_unit)) <- orig_vals_unit
    trans_vals_unit <- coef(po, names(orig_vals_unit), transform = TRUE)
    names(trans_vals_unit) <- paste0(seas[i], '_', names(orig_vals_unit))
    trans_vals <- c(trans_vals, trans_vals_unit)
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po, seas, params_all, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_shared <- trans_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared, transform = TRUE) <- trans_vals_shared
  orig_vals_shared <- coef(po, params_shared)
  orig_vals <- orig_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    trans_vals_unit <- trans_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(trans_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(trans_vals_unit))
    
    coef(po, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
    orig_vals_unit <- coef(po, names(trans_vals_unit))
    
    names(orig_vals_unit) <- paste0(seas[i], '_', names(trans_vals_unit))
    orig_vals <- c(orig_vals, orig_vals_unit)
  }
  
  return(orig_vals)
}

# ---------------------------------------------------------------------------------------------------------------------

# Fit using trajectory matching

# Loop through years and construct pomp models:
po_list <- vector('list', length(seasons))
for (yr in seasons) {
  print(yr)
  
  # Load data and create pomp object:
  source('src/resp_interaction_model.R')
  
  # Check whether appreciable activity for both viruses:
  if (sum(resp_mod@data[1, ]) <= 100) {
    print('Insufficient virus 1')
  }
  if (sum(resp_mod@data[2, ]) <= 100) {
    print('Insufficient virus 2')
  }
  
  if (sum(resp_mod@data[1, ]) > 100 & sum(resp_mod@data[2, ]) > 100) {
    po_list[[yr - (seasons[1] - 1)]] <- resp_mod
  }
  
}

# Remove empty elements:
seasons <- seasons[lapply(po_list, length) > 0]
po_list <- po_list[lapply(po_list, length) > 0]

# Choose parameters to estimate:
if (int_eff == 'susc') {
  shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_lambda1', 'theta_lambda2')
} else if (int_eff == 'sev') {
  shared_estpars <- c('rho1', 'rho2', 'delta', 'theta_rho1', 'theta_rho2')
} else {
  stop('Unrecognized int_eff value.')
}

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

unit_sp_estpars <- c()
for (i in 1:length(seasons)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons[i], unit_estpars, sep = '_'))
}
rm(i)
# unit_sp_estpars <- paste0(seasons, '_', rep(unit_estpars, each = length(seasons)))

true_estpars <- c(shared_estpars, unit_estpars)
estpars <- c(shared_estpars, unit_sp_estpars)

# Set upper/lower values for global params:
start_range <- data.frame(rho1 = c(0, 1),
                          rho2 = c(0, 1),
                          delta = c(delta_min, 7 / 2),
                          theta_lambda1 = c(0, 1),
                          theta_lambda2 = c(0, 1),
                          theta_rho1 = c(0, 1),
                          theta_rho2 = c(0, 1))

# Set upper/lower values for unit params (broad):
unit_start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                               Ri2 = c(1.0, Ri_max2),
                               I10 = c(0, 1e-3),
                               I20 = c(0, 1e-3),
                               R10 = c(0, 0.3),
                               R20 = c(0, 0.3),
                               R120 = c(0, 0.3))

# Get 99% CI from round 1 for unit params:
tj_res_list <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

tj_res_list <- tj_res_list[str_detect(names(tj_res_list), vir1)]

ci_list <- vector('list', length(seasons))

for (i in 1:length(ci_list)) {
  tj_res_temp <- tj_res_list[[i]] %>%
    select(all_of(unit_estpars))
  
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

# Get data frame of all ranges:
for (i in 1:length(seasons)) {
  
  if (search_type == 'broad') {
    start_range_temp <- unit_start_range
  } else if (search_type == 'round1_CIs') {
    start_range_temp <- ci_list[[i]]
  }
  
  names(start_range_temp) <- paste0(seasons[i], '_', names(unit_start_range))
  start_range <- start_range %>%
    bind_cols(start_range_temp)
  rm(start_range_temp)
  
}
rm(i)

start_range <- start_range[, estpars]

# Get starting values for each parameter:
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

print(start_range)
print(summary(start_values))
print(estpars)

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
  x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)
  
  # Check that parameter transformations correct:
  x0_orig <- back_transform_params(x0_trans, resp_mod, seasons, estpars, shared_estpars)
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
  
  
  
  
  
  x0_fit <- m$solution
  names(x0_fit) <- x0_trans_names
  print(-1 * calculate_global_loglik(x0_fit))
  back_transform_params(x0_fit, resp_mod, seasons, estpars, shared_estpars) %>% print()
  
  
  
}
rm(i)















# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

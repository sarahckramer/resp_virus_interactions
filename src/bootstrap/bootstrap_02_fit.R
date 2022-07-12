# ---------------------------------------------------------------------------------------------------------------------
# Run trajectory matching on synthetic data to obtain CIs for parameter value estimates
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
int_eff <- as.character(Sys.getenv("INTERACTIONEFFECT")); print(int_eff)
vir1 <- as.character(Sys.getenv("VIRUS1")); print(vir1)

# # Set parameters for local run:
# jobid <- 1
# no_jobs <- 1
# vir1 <- 'flu_h1' # 'flu_h1', 'flu_b'
# sobol_size <- 10
# int_eff <- 'susc' # 'susc' or 'sev' - fit impact of interaction on susceptibility or severity?

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
time_max <- 14.75 # Maximal execution time (in hours)

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

lag_val <- 0
prof_lik <- FALSE

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

# Load pomp models

# Loop through years and construct pomp models:
po_list <- vector('list', length(seasons))
for (yr_index in 1:length(seasons)) {
  yr <- seasons[yr_index]
  print(yr)
  
  # Load data and create pomp object:
  source('src/resp_interaction_model.R')
  
  # Check whether any data for given season:
  if (exists('resp_mod')) {
    
    # If doing profile likelihood, set interaction parameter:
    if (prof_lik) {
      coef(resp_mod, prof_param) <- prof_val
    }
    
    # Add pomp object to list:
    po_list[[yr_index]] <- resp_mod
    
  }
  
  # Remove pomp object before repeating loop:
  rm(resp_mod)
  
}

# Remove empty elements:
seasons <- seasons[lapply(po_list, length) > 0]
po_list <- po_list[lapply(po_list, length) > 0]

# ---------------------------------------------------------------------------------------------------------------------

# Assign synthetic data to pomp objects

# Read in synthetic flu/RSV data:
synth_LIST <- read_rds(paste0('results/synth_data_for_bootstrapping_', vir1, '.rds'))

# Set data in pomp objects to synthetic values:
for (i in 1:length(po_list)) {
  
  po_list[[i]]@data <- synth_LIST[[i]][[jobid]]
  
}

# Clean up synthetic data:
rm(synth_LIST)

# ---------------------------------------------------------------------------------------------------------------------

# Get parameter start values

# Choose parameters to estimate:
if (int_eff == 'susc') {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
} else if (int_eff == 'sev') {
  shared_estpars <- c('rho1', 'rho2', 'theta_rho1', 'theta_rho2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
} else {
  stop('Unrecognized int_eff value.')
}

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

unit_sp_estpars <- c()
for (i in 1:length(seasons)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons[i], unit_estpars, sep = '_'))
}
rm(i)

true_estpars <- c(shared_estpars, unit_estpars)
estpars <- c(shared_estpars, unit_sp_estpars)

# Get start ranges for all parameters:
if (vir1 == 'flu_h1') {
  start_range <- read_rds('results/round2_cis/round2CI_startvals_PROF_H1.rds')
} else if (vir1 == 'flu_b') {
  start_range <- read_rds('results/round2_cis/round2CI_startvals_PROF_B.rds')
} else {
  stop('Unknown vir1!')
}
start_range <- start_range[, estpars]

# Get starting values for each parameter:
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)
start_values <- start_values %>%
  mutate(phi = if_else(phi > 52.25, phi - 52.25, phi))

print(start_range)
print(summary(start_values))
print(estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Fit using trajectory matching

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# Set maximal execution time for each estimation:
nmins_exec <- time_max * 60 / (sobol_size / no_jobs)
print(sprintf("Max estimation time=%.1f min", nmins_exec))

# Get unique identifiers:
sub_start <- (1:(sobol_size / no_jobs))

# Fit:
for (i in seq_along(sub_start)) {
  
  print(paste0('Estimation: ', sub_start[i]))
  
  # Get start values:
  x0 <- as.numeric(start_values[sub_start[i], ])
  x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)
  
  # Check that parameter transformations correct:
  x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
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
    x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], seasons, estpars, shared_estpars)
    
    out <- list(estpars = x0_fit_untrans,
                ll = -m$objective,
                conv = m$status,
                message = m$message,
                niter = m$iterations,
                etime = as.numeric(etime))
    
    # Write to file:
    saveRDS(out, file = sprintf('results/res_%s_%s_%d_%d.rds',
                                vir1,
                                int_eff,
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
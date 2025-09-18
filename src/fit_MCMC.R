# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model using MCMC (through package BayesianTools)
# Round 2: Fit all seasons simultaneously, constraining shared parameters to be the same for all seasons
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(parallel)
library(doMC)
library(BayesianTools)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
which_round <- as.integer(Sys.getenv("WHICHROUND")); print(which_round)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
fit_canada <- as.logical(Sys.getenv("FITCANADA")); print(fit_canada)

# # Set parameters for local run:
# jobid <- 1
# no_jobs <- 10
# 
# sobol_size <- 500
# which_round <- 1
# search_type <- 'round1_CIs'
# fit_canada <- FALSE

# Set MCMC parameters:
n_chains <- 4
n_iter <- 5e4

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'
sens <- 'main'

if (fit_canada) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
if (fit_canada) {
  seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
}

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
  return(-glob_ll)
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

set_prior <- function(pars_nm) {
  # Sets up priors for all parameters being fit
  # param pars_nm: Names of all fit parameters
  # returns: Prior object for use with BayesianTools
  
  dens_fun <- function(par) {
    names(par) <- pars_nm
    
    get_seas <- unique(str_sub(pars_nm[str_detect(pars_nm, 's1')], 1, 6))
    
    d <- numeric(length(pars_nm))
    
    d['rho1'] <- dnorm(par['rho1'], mean = 0, sd = 10, log = TRUE)
    # d['rho1'] <- dnorm(par['rho1'], mean = -2, sd = 5, log = TRUE)
    d['rho2'] <- dnorm(par['rho2'], mean = 0, sd = 10, log = TRUE)
    d['theta_lambda1'] <- dnorm(par['theta_lambda1'], mean = 0, sd = 10, log = TRUE)
    d['theta_lambda2'] <- dnorm(par['theta_lambda2'], mean = 0, sd = 10, log = TRUE)
    d['delta1'] <- dnorm(par['delta1'], mean = 0, sd = 3, log = TRUE)
    d['d2'] <- dnorm(par['d2'], mean = -2.3, sd = 7, log = TRUE)
    d['alpha'] <- dnorm(par['alpha'], mean = 0, sd = 10, log = TRUE)
    d['phi'] <- dnorm(par['phi'], mean = 0, sd = 10, log = TRUE)
    d['eta_temp1'] <- dunif(par['eta_temp1'], min = -0.5, max = 0.5, log = TRUE)
    d['eta_temp2'] <- dunif(par['eta_temp2'], min = -0.5, max = 0.5, log = TRUE)
    d['eta_ah1'] <- dunif(par['eta_ah1'], min = -0.5, max = 0.5, log = TRUE)
    d['eta_ah2'] <- dunif(par['eta_ah2'], min = -0.5, max = 0.5, log = TRUE)
    
    for (seas in get_seas) {
      d[paste0(seas, '_Ri1')] <- dnorm(par[paste0(seas, '_Ri1')], mean = 0, sd = 10, log = TRUE)
      d[paste0(seas, '_Ri2')] <- dnorm(par[paste0(seas, '_Ri1')], mean = 0, sd = 10, log = TRUE)
      d[paste0(seas, '_I10')] <- dlnorm(-1 * par[paste0(seas, '_I10')], meanlog = 1.95, sdlog = 0.18, log = TRUE)
      d[paste0(seas, '_I20')] <- dlnorm(-1 * par[paste0(seas, '_I10')], meanlog = 1.95, sdlog = 0.18, log = TRUE)
      d[paste0(seas, '_R10')] <- dnorm(par[paste0(seas, '_R10')], mean = 0, sd = 10, log = TRUE)
      d[paste0(seas, '_R20')] <- dnorm(par[paste0(seas, '_R20')], mean = 0, sd = 10, log = TRUE)
      d[paste0(seas, '_R120')] <- dnorm(par[paste0(seas, '_R120')], mean = 0, sd = 10, log = TRUE)
    }
    
    return(sum(d))  # Return joint log-prior (sum of independent components)
  }
  
  samp_fun <- function(n = 1) {
    get_seas <- unique(str_sub(pars_nm[str_detect(pars_nm, 's1')], 1, 6))
    
    samp_mat <- matrix(NA, nrow = n, ncol = length(pars_nm))
    colnames(samp_mat) <- pars_nm
    
    samp_mat[, 'rho1'] <- rnorm(n, mean = 0, sd = 10)
    # samp_mat[, 'rho1'] <- rnorm(n, mean = -2, sd = 5)
    samp_mat[, 'rho2'] <- rnorm(n, mean = 0, sd = 10)
    samp_mat[, 'theta_lambda1'] <- rnorm(n, mean = 0, sd = 10)
    samp_mat[, 'theta_lambda2'] <- rnorm(n, mean = 0, sd = 10)
    samp_mat[, 'delta1'] <- rnorm(n, mean = 0, sd = 3)
    samp_mat[, 'd2'] <- rnorm(n, mean = -2.3, sd = 7)
    samp_mat[, 'alpha'] <- rnorm(n, mean = 0, sd = 10)
    samp_mat[, 'phi'] <- rnorm(n, mean = 0, sd = 10)
    samp_mat[, 'eta_temp1'] <- runif(n, min = -0.5, max = 0.5)
    samp_mat[, 'eta_temp2'] <- runif(n, min = -0.5, max = 0.5)
    samp_mat[, 'eta_ah1'] <- runif(n, min = -0.5, max = 0.5)
    samp_mat[, 'eta_ah2'] <- runif(n, min = -0.5, max = 0.5)
    
    for (seas in get_seas) {
      samp_mat[, paste0(seas, '_Ri1')] <- rnorm(n, mean = 0, sd = 10)
      samp_mat[, paste0(seas, '_Ri2')] <- rnorm(n, mean = 0, sd = 10)
      samp_mat[, paste0(seas, '_I10')] <- -1 * rlnorm(n, meanlog = 1.95, sdlog = 0.18)
      samp_mat[, paste0(seas, '_I20')] <- -1 * rlnorm(n, meanlog = 1.95, sdlog = 0.18)
      samp_mat[, paste0(seas, '_R10')] <- rnorm(n, mean = 0, sd = 10)
      samp_mat[, paste0(seas, '_R20')] <- rnorm(n, mean = 0, sd = 10)
      samp_mat[, paste0(seas, '_R120')] <- rnorm(n, mean = 0, sd = 10)
    }
    
    return(samp_mat)
  }
  
  # prior <- createPrior(
  #   density = dens_fun,   # Log-density function
  #   sampler = samp_fun,   # Sampling function
  #   lower   = bounds$inf, # Lower bounds for each parameter
  #   upper   = bounds$sup  # Upper bounds for each parameter
  # )
  
  prior <- createPrior(
    density = dens_fun,
    sampler = samp_fun,
    lower = c(rep(-Inf, 8), rep(-0.5, 4), rep(-Inf, 42)),
    upper = c(rep(Inf, 8), rep(0.5, 4), rep(Inf, 42))
  )
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Fit using trajectory matching

# Loop through years and construct pomp models:
po_list <- vector('list', length(seasons))
for (yr_index in 1:length(seasons)) {
  yr <- seasons[yr_index]
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

# Clean up:
rm(hk_dat, can_dat, us_dat, dat_pomp, age_structured, sens, yr_index, yr, nrow_check)

# Choose parameters to estimate:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

unit_sp_estpars <- c()
for (i in 1:length(seasons)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons[i], unit_estpars, sep = '_'))
}
rm(i)

true_estpars <- c(shared_estpars, unit_estpars)
estpars <- c(shared_estpars, unit_sp_estpars)

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

# Set upper/lower values for unit params (broad):
unit_start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                               Ri2 = c(1.0, Ri_max2),
                               I10 = c(0, 1e-3),
                               I20 = c(0, 1e-3),
                               R10 = c(0, 0.3),
                               R20 = c(0, 0.3),
                               R120 = c(0, 0.3))

# Get 95% CI from round 1 for unit params:
tj_res_list <- read_rds('results/round1_fitsharedFALSE/traj_match_round1_byvirseas_TOP.rds')
ci_list <- vector('list', length(seasons))

for (i in 1:length(ci_list)) {
  yr <- seasons[i]
  
  tj_res_temp <- tj_res_list[[which(str_detect(names(tj_res_list), yr))]] %>%
    select(all_of(unit_estpars)) %>% na.omit()
  
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
  
  ci_list[[i]] <- ci_temp
}
rm(i)

# Get data frame of all ranges:
if (search_type == 'round2_CIs') {
  
  start_range <- read_rds(paste0('results/round2_CIs/from_2_', which_round - 1, '/round2CI_startvals.rds'))
  
} else {
  
  for (i in 1:length(seasons)) {
    
    if (search_type == 'broad') {
      start_range_temp <- unit_start_range
    } else if (search_type == 'round1_CIs') {
      start_range_temp <- ci_list[[i]]
    } else if (search_type == 'round2_CIs') {
      print('ERROR: Round2 CIs used.')
    } else {
      stop('Unrecognized search type!')
    }
    
    names(start_range_temp) <- paste0(seasons[i], '_', names(unit_start_range))
    start_range <- start_range %>%
      bind_cols(start_range_temp)
    rm(start_range_temp)
    
  }
  rm(i)
  
  start_range <- start_range[, estpars]
  
}

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

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# Get unique identifiers:
sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)
print(sub_start)

# Fit:
for (i in seq_along(sub_start)) {
  
  print(paste0('Estimation: ', sub_start[i]))
  
  
  x0 <- as.numeric(start_values[sub_start[i], ])
  x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)
  
  # Check that parameter transformations correct:
  x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
  expect_equal(x0, unname(x0_orig))
  rm(x0_orig)
  
  # Calculate initial log-likelihood:
  print(calculate_global_loglik(x0_trans))
  
  # Get initial chain values:
  start_chains <- replicate(n = n_chains, expr = x0_trans) %>% t()
  
  # Get prior:
  pr <- set_prior(pars_nm = estpars)
  
  # Set up models:
  bayesianSetup <- createBayesianSetup(
    likelihood = calculate_global_loglik,
    prior = pr,
    names = estpars
  )
  
  # Fit:
  tic <- Sys.time()
  m <- runMCMC(bayesianSetup = bayesianSetup,
               sampler = 'DEzs',
               settings = list(
                 iterations = n_chains * n_iter,
                 startValue = start_chains
               ))
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'hours'
  print(etime)
  
  # Extract results:
  est_mcmc_det <- getSample(sampler = m, coda = TRUE, start = n_iter / 2)
  
  dic <- DIC(sampler = m, start = n_iter / 2)
  mle <- MAP(m)
  
  print(mle)
  print(MCMCvis::MCMCsummary(object = est_mcmc_det, round = 4))
  
  # MCMCvis::MCMCtrace(object = est_mcmc_det, pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE)
  # samples_df <- as.data.frame(getSample(m, start = n_iter / 2, coda = FALSE))
  
  # For now, save whole fit object:
  saveRDS(m, file = sprintf('results/mod_%s_%d_%d.rds',
                            vir1,
                            jobid,
                            sub_start[i])
  )
  
}
rm(i)

# Clean up:
rm(list = ls())

print('Done!')

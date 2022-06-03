# ---------------------------------------------------------------------------------------------------------------------
# Code to fit flu/RSV interaction model using MCMC
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(MCMCpack)
library(fmcmc)
library(matrixcalc)
library(distr)

# Get cluster environmental variables:
# jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
# no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
# sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
int_eff <- as.character(Sys.getenv("INTERACTIONEFFECT")); print(int_eff)
vir1 <- as.character(Sys.getenv("VIRUS1")); print(vir1)
num_chains <- as.integer(Sys.getenv("CHAINS")); print(num_chains)
burnin_val <- as.integer(Sys.getenv("BURNIN")); print(burnin_val)
mcmc_val <- as.integer(Sys.getenv("MCMC")); print(mcmc_val)
thin_val <- as.integer(Sys.getenv("THIN")); print(thin_val)
tune_val <- as.numeric(Sys.getenv("TUNE")); print(tune_val)

# # Set parameters for local run:
# # jobid <- 1
# # no_jobs <- 1
# vir1 <- 'flu_h1'
# 
# # sobol_size <- 5
# int_eff <- 'susc' # 'susc' or 'sev' - fit impact of interaction on susceptibility or severity?
# 
# num_chains <- 4
# burnin_val <- 1000
# mcmc_val <- 10000
# thin_val <- 10
# tune_val <- 1.0

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

lag_val <- 0

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

calculate_global_loglik <- function(trans_vals, x0_trans_names = x0_trans_names, seasons = seasons, obj_fun_list = obj_fun_list) {
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
  return(-1 * glob_ll)
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

# Prep models and parameters for fitting

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

# Remove empty elements:
seasons <- seasons[lapply(po_list, length) > 0]
po_list <- po_list[lapply(po_list, length) > 0]

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

# Get starting parameter sets (top 5 fits from trajectory matching):
if (vir1 == 'flu_h1') {
  start_values <- read_rds('results/round2_cis/round2_topfits_flu_h1.rds')
} else if (vir1 == 'flu_b') {
  start_values <- read_rds('results/round2_cis/round2_topfits_flu_b.rds')
} else {
  stop('Unknown vir1!')
}

start_values <- start_values[, estpars]

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# Where start values are very low, set to higher value to allow better mixing:
for (i in 1:nrow(start_values)) {

  x0 <- as.numeric(start_values[i, ])
  x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)

  x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
  expect_equal(x0, unname(x0_orig))
  rm(x0_orig)

  # start_values[i, ][abs(x0) < 0.005 &
  #                     (str_detect(x0_trans_names, 'theta') |
  #                        (str_detect(x0_trans_names, 'R') & str_detect(x0_trans_names, '0')))] <- 0.005
  start_values[i, 3:4][abs(x0)[3:4] < 0.005] <- 0.005

  # # Check that this does not result in total initial conditions > 1:
  # for (season in seasons) {
  #
  #   vals_temp <- start_values[i, str_detect(names(start_values), season) &
  #                               str_detect(names(start_values), '0')]
  #
  #   if (sum(vals_temp) > 0.99) {
  #     print(paste0('Initial conditions in season ', season, ' greater than 1!'))
  #   }
  #
  # }

  # # Determine a good cutoff value where likelihoods begin to change:
  # par(mfrow = c(3, 2))
  #
  # low_params <- names(start_values[i, ][abs(x0) < 0.005 &
  #                                         (str_detect(x0_trans_names, 'theta') |
  #                                            (str_detect(x0_trans_names, 'R') & str_detect(x0_trans_names, '0')))])
  # test_vals <- c(1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, seq(0.01, 0.05, by = 0.01))
  #
  # for (param in low_params) {
  #
  #   ll <- c()
  #
  #   for (val in test_vals) {
  #     x0_temp <- x0
  #     x0_temp[which(x0_trans_names == param)] <- val
  #
  #     x0_trans_temp <- transform_params(x0_temp, po_list[[1]], seasons, estpars, shared_estpars)
  #
  #     ll <- c(ll, calculate_global_loglik(x0_trans_temp, x0_trans_names, seasons, obj_fun_list))
  #   }
  #
  #   plot(test_vals, ll, type = 'b', pch = 20, main = param)
  #
  # }

}

print(summary(start_values))
print(estpars)

# ---------------------------------------------------------------------------------------------------------------------

# # Set priors on interaction parameters
# 
# # Function for PDF of log(exp) distribution:
# log_exp <- function(x, lambda = 5) {
#   lambda * exp(x - lambda * exp(x))
# }
# 
# dist_log_exp <- AbscontDistribution(d = log_exp, withStand = TRUE)
# # rdist_log_exp <- r(dist_log_exp)
# ddist_log_exp <- d(dist_log_exp)
# 
# # Get prior likelihood:
# prior_like <- function(trans_vals, x0_trans_names = x0_trans_names) {
#   names(trans_vals) <- x0_trans_names
#   
#   str1 <- trans_vals['theta_lambda1']
#   str2 <- trans_vals['theta_lambda2']
#   dur1 <- trans_vals['delta1']
#   dur2 <- trans_vals['d2']
#   
#   prior_thetalambda1 <- dnorm(str1, mean = 0, sd = 1.6, log = TRUE)
#   prior_thetalambda2 <- dnorm(str2, mean = 0, sd = 1.6, log = TRUE)
#   prior_delta1 <- ddist_log_exp(dur1) %>% log()
#   prior_d2 <- dnorm(dur2, mean = 0, sd = 1.5, log = TRUE)
#   
#   ll <- prior_thetalambda1 + prior_thetalambda2 + prior_delta1 + prior_d2
#   return(unname(ll))
# }
# 
# # Get posterior likelihood:
# posterior_like <- function(trans_vals, x0_trans_names = x0_trans_names, seasons = seasons, obj_fun_list = obj_fun_list) {
#   
#   return(calculate_global_loglik(trans_vals, x0_trans_names, seasons, obj_fun_list) + prior_like(trans_vals, x0_trans_names))
#   
# }

# ---------------------------------------------------------------------------------------------------------------------

# # Calculate variance-covariance matrix for MCMC proposals
# 
# # Get maximum likelihood estimate:
# mle <- start_values[1, ]
# expect_true(all(names(mle) == names(start_values)))
# 
# # Get estimate of Hessian matrix:
# mle_trans <- transform_params(as.numeric(mle), po_list[[1]], seasons, estpars, shared_estpars)
# mle_orig <- back_transform_params(mle_trans, po_list[[1]], seasons, estpars, shared_estpars)
# expect_equal(as.numeric(mle), unname(mle_orig))
# rm(mle_orig)
# 
# mle_trans_names <- names(mle_trans)
# print(calculate_global_loglik(mle_trans, mle_trans_names, seasons, obj_fun_list))
# 
# fit_w_hessian <- optim(par = mle_trans,
#                        fn = calculate_global_loglik,
#                        x0_trans_names = mle_trans_names,
#                        seasons = seasons,
#                        obj_fun_list = obj_fun_list,
#                        method = 'BFGS',
#                        hessian = TRUE,
#                        control = list(maxit = 0))
# H_mat <- fit_w_hessian$hessian
# is.negative.definite(H_mat)
# 
# # # Try different methods for calculating Hessian:
# # library(numDeriv)
# # H_mat2 <- numDeriv::hessian(func = calculate_global_loglik,
# #                             x0_trans_names = mle_trans_names,
# #                             seasons = seasons,
# #                             obj_fun_list = obj_fun_list,
# #                             x = mle_trans,
# #                             method = 'Richardson')
# # is.negative.definite(H_mat2)
# 
# # Get variance-covariance matrix:
# hess.new <- H_mat
# CC <- NULL
# if (max(diag(hess.new) == 0)) {
#   for (i in 1:nrow(hess.new)) {
#     if (hess.new[i, i] == 0) {
#       hess.new[i, i] <- -1e-06
#     }
#   }
# }
# while (is.null(CC)) {
#   hess.flag <- 1
#   hess.new <- hess.new - diag(diag(0.01 * abs(H_mat)))
#   try(CC <- chol(-1 * hess.new), silent = TRUE)
# }
# rm(i)
# is.negative.definite(hess.new)
# 
# T_mat <- diag(1.0, nrow = nrow(H_mat))
# V_mat <- T_mat %*% solve(-1 * hess.new) %*% T_mat
# 
# # V = T (-1 x H)^(-1) T
# # T is the diagonal positive definite matrix formed from "tune"
# # H is the approximate Hessian of "fun" evaluated at its mode

# ---------------------------------------------------------------------------------------------------------------------

# # Fit using MCMC
# 
# # Loop through start sets and fit:
# for (i in 1:num_chains) {
# 
#   print(paste0('Estimation: ', i))
# 
#   # Get start values:
#   x0 <- as.numeric(start_values[i, ])
#   x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
#   x0_trans_names <- names(x0_trans)
# 
#   # Check that parameter transformations correct:
#   x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
#   expect_equal(x0, unname(x0_orig))
#   rm(x0_orig)
# 
#   # Calculate initial log-likelihood:
#   print(calculate_global_loglik(x0_trans, x0_trans_names, seasons, obj_fun_list))
#   
#   # Fit models:
#   tic <- Sys.time()
#   m <- MCMCmetrop1R(fun = calculate_global_loglik,
#                     x0_trans_names = x0_trans_names,
#                     seasons = seasons,
#                     obj_fun_list = obj_fun_list,
#                     theta.init = x0_trans,
#                     burnin = burnin_val,
#                     mcmc = mcmc_val,
#                     thin = thin_val,
#                     tune = tune_val,
#                     verbose = round(mcmc_val / 4),
#                     logfun = TRUE,
#                     # force.samp = TRUE,
#                     # optim.method = 'BFGS',
#                     # optim.control = list(maxit = 0),
#                     V = V_mat)
#   toc <- Sys.time()
#   etime <- toc - tic
#   units(etime) <- 'hours'
#   print(etime)
# 
#   # Save results to file:
#   saveRDS(m, file = sprintf('results/res_mcmc_%s_%s_%d_%d_%d_%d_%f.rds',
#                             vir1,
#                             int_eff,
#                             i,
#                             mcmc_val,
#                             burnin_val,
#                             thin_val,
#                             tune_val))
# 
#   # Print results:
#   print(calculate_global_loglik(m[nrow(m), ], x0_trans_names, seasons, obj_fun_list))
#   print(back_transform_params(colMeans(m), po_list[[1]], seasons, estpars, shared_estpars))
# 
# }
# rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Fit using fmcmc package

# Set up matrix of initial values:
init_trans <- matrix(NA, nrow = num_chains, ncol = length(estpars))

for (i in 1:num_chains) {
  x0 <- as.numeric(start_values[i, ])
  x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
  x0_trans_names <- names(x0_trans)

  x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
  expect_equal(x0, unname(x0_orig))
  rm(x0_orig)

  print(calculate_global_loglik(x0_trans, x0_trans_names, seasons, obj_fun_list))
  # print(posterior_like(x0_trans, x0_trans_names, seasons, obj_fun_list))

  init_trans[i, ] <- x0_trans
}
rm(i)

# Add column names:
colnames(init_trans) <- x0_trans_names

# Fit model:
tic <- Sys.time()
m <- MCMC(initial = init_trans[1, ],
          fun = calculate_global_loglik,
          nsteps = mcmc_val + burnin_val,
          x0_trans_names = x0_trans_names,
          seasons = seasons,
          obj_fun_list = obj_fun_list,
          nchains = num_chains,
          burnin = burnin_val,
          thin = thin_val,
          # kernel = kernel_normal(scale = 0.01),
          kernel = kernel_ram(),
          multicore = FALSE,
          conv_checker = convergence_gelman(10000))
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'hours'
print(etime)

# Save results to file:
saveRDS(m, file = sprintf('results/res_mcmc_%s_%s_%d_%d_%d.rds',
                          vir1,
                          int_eff,
                          mcmc_val,
                          burnin_val,
                          thin_val))

# Save resulting kernel:
kernel_store <- get_kernel()
saveRDS(kernel_store, file = sprintf('results/kernel_ram_%s_%s_%d_%d_%d.rds',
                                     vir1,
                                     int_eff,
                                     mcmc_val,
                                     burnin_val,
                                     thin_val))

# Print results:
for (i in 1:num_chains) {
  print(calculate_global_loglik(m[[i]][nrow(m[[i]]), ], x0_trans_names, seasons, obj_fun_list))
}
for (i in 1:num_chains) {
  print(back_transform_params(colMeans(m[[i]]), po_list[[1]], seasons, estpars, shared_estpars))
}
rm(i)

print('Done!')

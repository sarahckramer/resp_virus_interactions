# ---------------------------------------------------------------------------------------------------------------------
# Code to assess fit and convergence of MCMC runs
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(coda)
library(distr)
library(RColorBrewer)

# Set parameters for pomp:
vir1 <- 'flu_h1'
vir2 <- 'rsv'

debug_bool <- FALSE
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

# Function for PDF of log(exp) distribution:
log_exp <- function(x, lambda = 5) {
  lambda * exp(x - lambda * exp(x))
}

dist_log_exp <- AbscontDistribution(d = log_exp, withStand = TRUE)
ddist_log_exp <- d(dist_log_exp)

# Get prior likelihood:
prior_like <- function(trans_vals, x0_trans_names = x0_trans_names) {
  names(trans_vals) <- x0_trans_names
  
  str1 <- trans_vals['theta_lambda1']
  str2 <- trans_vals['theta_lambda2']
  dur1 <- trans_vals['delta1']
  dur2 <- trans_vals['d2']
  
  prior_thetalambda1 <- dnorm(str1, mean = 0, sd = 1.6, log = TRUE)
  prior_thetalambda2 <- dnorm(str2, mean = 0, sd = 1.6, log = TRUE)
  prior_delta1 <- ddist_log_exp(dur1) %>% log()
  prior_d2 <- dnorm(dur2, mean = 0, sd = 1.5, log = TRUE)
  
  ll <- prior_thetalambda1 + prior_thetalambda2 + prior_delta1 + prior_d2
  return(unname(ll))
}

# Get posterior likelihood:
posterior_like <- function(trans_vals, x0_trans_names = x0_trans_names, seasons = seasons, obj_fun_list = obj_fun_list) {
  
  return(calculate_global_loglik(trans_vals, x0_trans_names, seasons, obj_fun_list) + prior_like(trans_vals, x0_trans_names))
  
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

# Load results

# Read in MCMCpack results:
res_mcmc1 = res_mcmc2 = res_mcmc3 = vector('list', length = 4)
for (i in 1:length(res_mcmc1)) {
  res_mcmc1[[i]] <- read_rds(paste0('results/mcmc/res_mcmc_flu_h1_susc_350000_50000_10_0.050000/res_mcmc_flu_h1_susc_', i, '_350000_50000_10_0.050000.rds'))
  res_mcmc2[[i]] <- read_rds(paste0('results/mcmc/res_mcmc_flu_h1_susc_500000_100000_10_0.02/res_mcmc_flu_h1_susc_', i, '_500000_100000_10_0.020000.rds'))
  res_mcmc3[[i]] <- read_rds(paste0('results/mcmc/res_mcmc_flu_h1_susc_500000_100000_10_0.5/res_mcmc_flu_h1_susc_', i, '_500000_100000_10_0.500000.rds'))
}

res_mcmc1 <- as.mcmc.list(res_mcmc1)
res_mcmc2 <- as.mcmc.list(res_mcmc2)
res_mcmc3 <- as.mcmc.list(res_mcmc3)

# Read in fmcmc results:
res_fmcmc <- read_rds('results/mcmc/res_mcmc_flu_h1_susc_350000_100000_10.rds')
res_fmcmc_old <- read_rds('results/mcmc/res_mcmc_flu_h1_susc_300000_10000_10.rds')

# Get names of estimated parameters:
estpars <- colnames(res_fmcmc[[1]])
shared_estpars <- estpars[which(!str_detect(estpars, 's1'))]
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# ---------------------------------------------------------------------------------------------------------------------

# Assess fit

# Get list of all results:
res_list <- list(res_mcmc1, res_mcmc2, res_mcmc3, res_fmcmc, res_fmcmc_old)

# Visualize log-likelihood over time for each chain:
logliks_fxn_list = logliks_post_list = vector('list', length = length(res_list))
thin_for_plotting <- 100
for (i in 1:length(res_list)) {
  
  logliks_fxn_temp = logliks_post_temp = matrix(NA, nrow = nrow(res_list[[i]][[1]]) / thin_for_plotting, ncol = length(res_list[[i]]))
  
  for (j in 1:length(res_list[[i]])) {
    
    logliks_fxn_temp[, j] <- apply(res_list[[i]][[j]][seq(thin_for_plotting, nrow(res_list[[i]][[j]]), by = thin_for_plotting), ], 1, calculate_global_loglik,
                                   x0_trans_names = estpars, seasons = seasons, obj_fun_list = obj_fun_list)
    logliks_post_temp[, j] <- apply(res_list[[i]][[j]][seq(thin_for_plotting, nrow(res_list[[i]][[j]]), by = thin_for_plotting), ], 1, posterior_like,
                                    x0_trans_names = estpars, seasons = seasons, obj_fun_list = obj_fun_list)
    
    # logliks_fxn_temp[, j] <- apply(res_list[[i]][[j]], 1, calculate_global_loglik,
    #                                x0_trans_names = estpars, seasons = seasons, obj_fun_list = obj_fun_list)
    # logliks_post_temp[, j] <- apply(res_list[[i]][[j]], 1, posterior_like,
    #                                 x0_trans_names = estpars, seasons = seasons, obj_fun_list = obj_fun_list)
    
  }
  
  logliks_fxn_list[[i]] <- logliks_fxn_temp
  logliks_post_list[[i]] <- logliks_post_temp
  rm(logliks_fxn_temp, logliks_post_temp)
  
}
rm(i, j)

par(mfrow = c(3, 2))
for (i in 1:length(logliks_fxn_list)) {
  matplot(logliks_fxn_list[[i]], type = 'b', pch = 20)
}
par(mfrow = c(3, 2))
for (i in 1:length(logliks_post_list)) {
  matplot(logliks_post_list[[i]], type = 'b', pch = 20)
}
rm(i)

# Un-transform fit parameter values:
res_list_untransformed <- vector('list', length = length(res_list))
for (i in 1:length(res_list)) {
  
  res_list_untransformed_temp <- res_list[[i]]
  
  for (j in 1:length(res_list[[i]])) {
    
    res_list_untransformed_temp[[j]] <- apply(res_list[[i]][[j]], 1, back_transform_params,
                                              po = po_list[[1]], seas = seasons, params_all = estpars, params_shared = shared_estpars) %>%
      t() %>%
      mcmc()
    
  }
  
  res_list_untransformed[[i]] <- res_list_untransformed_temp
  rm(res_list_untransformed_temp)
  
}
rm(i, j)

# Calculate effective sample size:
for (i in 1:length(res_list)) {
  print(effectiveSize(res_list[[i]]))
}
rm(i)

# Check correlations between parameters:
pdf('results/mcmc/pairs_plots.pdf', width = 60, height = 45)
for (i in 1:length(res_list)) {
  
  pairs(res_list[[i]][[1]][seq(thin_for_plotting, nrow(res_list[[i]][[1]]), by = thin_for_plotting), ])
  # for (j in 1:length(res_list[[i]])) {
  #   pairs(res_list[[i]][[j]][seq(thin_for_plotting, nrow(res_list[[i]][[j]]), by = thin_for_plotting), ])
  # }
  # Produces a very large file - plot only first chain for each run for now
  
}
dev.off()
rm(i)

pdf('results/mcmc/pairs_plots_alt.pdf', width = 30, height = 25)
for (i in 1:length(res_list)) {
  crosscorr.plot(res_list[[i]], col = brewer.pal(10, 'RdBu'))
}
dev.off()
rm(i)

# for (i in 1:length(res_list)) {
#   crosscorr(res_list[[i]]) %>% print()
# }
# rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Assess convergence

# Plot trace plots of transformed parameter values:
pdf('results/mcmc/traceplots_flu_h1_susc_mcmc1.pdf', width = 10, height = 10)
plot(res_list[[1]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_mcmc2.pdf', width = 10, height = 10)
plot(res_list[[2]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_mcmc3.pdf', width = 10, height = 10)
plot(res_list[[3]])
dev.off()

# pdf('results/mcmc/traceplots_flu_h1_susc_fmcmc.pdf', width = 10, height = 10)
# plot(res_list[[4]])
# dev.off()
# 
# pdf('results/mcmc/traceplots_flu_h1_susc_fmcmc_old.pdf', width = 10, height = 10)
# plot(res_list[[5]])
# dev.off()

# Plot trace plots of un-transformed parameter values:

pdf('results/mcmc/traceplots_flu_h1_susc_mcmc1_UNTRANSFORMED.pdf', width = 10, height = 10)
plot(res_list_untransformed[[1]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_mcmc2_UNTRANSFORMED.pdf', width = 10, height = 10)
plot(res_list_untransformed[[2]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_mcmc3_UNTRANSFORMED.pdf', width = 10, height = 10)
plot(res_list_untransformed[[3]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_fmcmc_UNTRANSFORMED.pdf', width = 10, height = 10)
plot(res_list_untransformed[[4]])
dev.off()

pdf('results/mcmc/traceplots_flu_h1_susc_fmcmc_old_UNTRANSFORMED.pdf', width = 10, height = 10)
plot(res_list_untransformed[[5]])
dev.off()

# Calculate Gelman-Rubin diagnostic to check for convergence:
pdf('results/mcmc/gelman_plots.pdf', width = 10, height = 10)
for (i in 1:length(res_list)) {
  gelman.plot(res_list[[i]])
}
dev.off()

for (i in 1:length(res_list)) {
  gelman.diag(res_list[[i]])$mpsrf %>% print()
  # Same value output by fmcmc convergence checker
}
rm(i)

# Check autocorrelation of chains:
for (i in 1:length(res_list)) {
  autocorr.plot(res_list[[i]])
}

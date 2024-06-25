# ---------------------------------------------------------------------------------------------------------------------
# Demo of model fitting procedure

# Tested on Windows 10, in R version 4.4.0
# Run time: ~___
# Required packages: tidyverse, testthat, pomp, nloptr, 
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Get run start time:
tic <- Sys.time()

# Load libraries:
library(tidyverse)
library(testthat)
library(pomp)
library(nloptr)

# ---------------------------------------------------------------------------------------------------------------------

# Build model

# Set relevant parameters:
sens <- 'main' # main analysis, not a sensitivity analysis
sobol_size <- 500 # how many starting parameter sets to use for fitting?

# List seasons for which data are available:
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

# # Get synthetic data set:
# synth_LIST <- read_rds('results/data_for_synthetic_testing.rds')
# 
# seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
# 
# dat_list <- vector('list', length = 6)
# 
# for (i in 1:length(seasons)) {
#   
#   yr <- seasons[i]
#   hk_dat_temp <- hk_dat %>% filter(season == yr)
#   synth_dat_temp <- synth_LIST[[1]][[i]][[1]] %>% t()
#   
#   hk_dat_temp[, c('n_P1', 'n_P2')] <- synth_dat_temp
#   dat_list[[i]] <- hk_dat_temp
#   
# }
# 
# synth_dat <- bind_rows(dat_list)
# 
# dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
# 
# synth_dat <- synth_dat %>%
#   rename('i_ILI' = 'GOPC') %>%
#   mutate(i_ILI = i_ILI / 1000) %>%
#   inner_join(dat_clim,
#              by = c('Year' = 'year',
#                     'Week' = 'week')) %>%
#   select(time:pop, temp, ah, rh) %>%
#   mutate(h3_inc = 0,
#          rhino_inc = 0)
# 
# write_csv(synth_dat, file = 'src/demo/demo_data.csv')

# Read in data:
df <- read_csv('src/demo/demo_data.csv')
# Note: these are synthetic data, generated at the MLEs obtained from the Hong Kong data

# n_T: number of tests performed for flu and RSV
# n_P1: number of observed influenza cases
# n_P2: number of observed RSV cases
# i_ILI: proportion of consultations at public out-patient clinics due to influenza-like illness
# pop: population of Hong Kong
# temp: mean temperature (standardized)
# ah: absolute humidity (standardized)

# Data frames used for analysis also require columns for the standardized proportion of
# tests positive for H3N2 (h3_inc) and rhinovirus (rhino_inc); these are set to zero since
# these columns are not used in the demo.

# Load required functions:
source('src/functions/functions_flu_RSV.R')
# source('src/functions/test_code.R')

# Create models for each season:
po_list <- vector('list', length(seasons))

for (yr_index in 1:length(seasons)) {
  
  yr <- seasons[yr_index]
  print(yr)
  
  df_temp <- df %>% filter(season == yr)
  
  resp_mod <- create_SITRxSITR_mod(dat = df_temp,
                                   Ri1_max = 2.0,
                                   Ri2_max = 3.0,
                                   d2_max = 10.0,
                                   debug_bool = FALSE,
                                   sens = sens,
                                   loc = 'hk')
  
  po_list[[yr_index]] <- resp_mod
  
  rm(resp_mod)
  
}
rm(yr_index, yr, df_temp, resp_mod)

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# ---------------------------------------------------------------------------------------------------------------------

# Fit models to each season separately ("round 1")

# For the purposes of the demo, we will fit to only one season:
resp_mod <- po_list[[3]] # model and data from 2015-16 season

# Set seed:
set.seed(749501349)

# Get starting values of all parameters:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')

start_range <- data.frame(Ri1 = c(1.0, 2.0),
                          Ri2 = c(1.0, 3.0),
                          I10 = c(0, 1e-3),
                          I20 = c(0, 1e-3),
                          R10 = c(0, 0.3),
                          R20 = c(0, 0.3),
                          R120 = c(0, 0.3),
                          rho1 = c(0, 1.0),
                          rho2 = c(0, 1.0),
                          theta_lambda1 = c(0, 1.0),
                          theta_lambda2 = c(0, 1.0),
                          theta_rho1 = c(0, 1.0),
                          theta_rho2 = c(0, 1.0),
                          delta1 = c(7 / 60, 7),
                          d2 = c(7 / 60, 7),
                          alpha = c(0, 0.5),
                          phi = c(0, 52.25),
                          eta_temp1 = c(-0.5, 0.5),
                          eta_temp2 = c(-0.5, 0.5),
                          eta_ah1 = c(-0.5, 0.5),
                          eta_ah2 = c(-0.5, 0.5))
start_range <- start_range[, estpars]

start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

# Create objective function for call to nloptr:
obj_fun <- traj_objfun(data = resp_mod,
                       est = estpars,
                       partrans = resp_mod@partrans,
                       verbose = TRUE)

# Fit for first 10 starting parameter sets:
sub_start <- 1:10
res_list <- vector('list', length = 10)

for (i in seq_along(sub_start)) {
  
  print(paste0('Estimation: ', sub_start[i]))
  
  # Get param start values:
  x0 <- as.numeric(start_values[sub_start[i], ])
  coef(resp_mod, estpars) <- x0
  x0_trans <- coef(resp_mod, estpars, transform = TRUE)
  
  # Run trajectory matching using subplex algorithm:
  # http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
  m <- try(
    nloptr(x0 = unname(x0_trans),
           eval_f = obj_fun,
           opts = list(algorithm = 'NLOPT_LN_SBPLX',
                       maxtime = 120.0,
                       maxeval = -1, # disabled
                       xtol_rel = -1, # disabled; default: 1e-4
                       print_level = 0))
  )
  
  # If estimation is successful, save results:
  if (!inherits(m, 'try-error')) {
    coef(resp_mod, estpars, transform = TRUE) <- m$solution
    
    # Collect all results:
    out <- list(estpars = coef(resp_mod, estpars),
                ll = -m$objective,
                message = m$message)
    
    # Write to file:
    res_list[[i]] <- out
    
    # Print results:
    print(out$ll)
    print(out$estpars, digits = 2)
    print(out$message)
  }
  
}

# Clean up:
rm(start_range, start_values, sub_start, i, x0, x0_trans, m, out)

# ---------------------------------------------------------------------------------------------------------------------

# Process "round 1" results

# Compile:
res_r1 <- lapply(res_list, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_list, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_list, getElement, 'message') %>%
              unlist()) %>%
  mutate(year = 's15-16') %>%
  select(year, Ri1:message)
rm(res_list)

# Remove where no convergence occurs:
res_r1 <- res_r1 %>%
  filter(!str_detect(message, 'maxtime')) %>%
  select(-message)

# Check that fit parameters do not produce simulations where state variables go negative:
p_mat <- parmat(coef(resp_mod), nrep = 10)

for (param in estpars) {
  p_mat[param, ] <- res_r1 %>% pull(param)
}

expect_equal(trajectory(resp_mod, params = p_mat, format = 'data.frame') %>%
               pivot_longer(X_SS:H2, names_to = 'state') %>%
               filter(value < 0) %>%
               nrow(), 0)

rm(p_mat, param)

# Arrange results by log-likelihood and keep only best estimates:
res_r1 <- res_r1 %>%
  arrange(desc(loglik))

no_best <- nrow(subset(res_r1, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = length(estpars))))
res_r1 <- res_r1[1:no_best, ]
# here all 10 fall within an acceptable range of the max log likelihood; in the real analysis,
# this is done with fits from all 500 starting parameter sets

# Clean up:
rm(resp_mod, estpars, no_best)

# ---------------------------------------------------------------------------------------------------------------------

# Fit models to all seasons simultaneously ("round 2")

# Load functions used to calculate log-likelihood of all seasons:
source('src/demo/fxns_demo.R')

# List all parameters to estimate:
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

# Get starting values of all parameters:

# Typically, we would get the start ranges of all parameters either from the "round 1" results
# above, or from a previous attempt to fit all seasons simultaneously. We would then draw 500
# starting parameter sets using the "sobol_design" function, as we did for round 1 above.
# Here, since we do not have time to run the full fitting procedure, we will use starting values
# that we know yield good results. Run time is about 20-30 minutes.

x0 <- read_rds('src/demo/starting_params.rds')

# set.seed(749501349)
# start_range <- read_rds('results/synthetic_testing/set_1/round2_CIs/from_2_4/round2CI_startvals.rds')
# start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
#                              upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
#                              nseq = sobol_size)
# 
# x0 <- as.numeric(start_values[485, ])
# # x0 <- as.numeric(out$estpars)
# write_rds(x0, file = 'src/demo/starting_params.rds')

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
})

# Fit:
x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
x0_trans_names <- names(x0_trans)

x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
expect_equal(x0, unname(x0_orig))
rm(x0_orig)

m <- try(
  nloptr(x0 = x0_trans, 
         eval_f = calculate_global_loglik,
         opts = list(algorithm = "NLOPT_LN_SBPLX",
                     maxtime = 60 * 60, # Max run time: 60 minutes
                     maxeval = -1, # Negative value: criterion is disabled
                     xtol_rel = -1, # Default value: 1e-4
                     print_level = 0))
)

if (!inherits(m, 'try-error')) {
  x0_fit <- m$solution
  names(x0_fit) <- x0_trans_names
  x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], seasons, estpars, shared_estpars)
  
  out <- list(estpars = x0_fit_untrans,
              ll = -m$objective,
              message = m$message)
  
  print(out$ll)
  print(out$estpars, digits = 2)
  print(out$message)
}

# Clean up:
rm(unit_estpars, unit_sp_estpars, x0_trans, x0_trans_names, m, x0_fit, x0_fit_untrans)

# ---------------------------------------------------------------------------------------------------------------------

# Process "round 2" results

# Reformat:
res_r2 <- out$estpars %>%
  t() %>%
  bind_cols('loglik' = out$ll) %>%
  bind_cols('message' = out$message)

# Check for convergence:
expect_false(str_detect(res_r2$message, 'maxtime'))

# Check that fit parameters do not produce simulations where state variables go negative:
traj_list <- lapply(1:length(seasons), function(ix) {
  
  pars_temp <- res_r2 %>% select(all_of(shared_estpars), contains(seasons[ix]))
  names(pars_temp) <- true_estpars
  
  coef(po_list[[ix]], true_estpars) <- pars_temp
  
  trajectory(po_list[[ix]], format = 'data.frame') %>%
    pivot_longer(X_SS:H2, names_to = 'state')
  
})

expect_false(any(lapply(traj_list, function(ix) {
  any(ix[, 'value'] < 0)
}) %>%
  unlist()))

# With only one result, there is not much to do here. Typically, the best fits from the 500
# starting parameter sets would be used as the starting parameter sets for further rounds of
# fitting, or, if the MLE has been reached, as the parameters used to generate synthetic data
# for parametric bootstrapping.

# ---------------------------------------------------------------------------------------------------------------------

# Bootstrapping

# Set seed:
set.seed(9489703)

# Set key parameters:
n <- 500 # How many synthetic datasets to create?

# Get MLE from "round 2" fitting:
# Here, we will use the parameter values fit in the previous step.
print(res_r2)

# For each season, set parameter values in pomp models to MLE:
for (i in 1:length(seasons)) {
  
  season <- seasons[i]
  
  res_temp <- res_r2 %>%
    select(all_of(shared_estpars), contains(season)) %>%
    rename_with(~ str_remove(.x, paste0(season, '_')))
  
  coef(po_list[[i]], names(res_temp)) <- res_temp
  
  rm(res_temp)
  
}
rm(i)

# Generate synthetic data:
synth_LIST <- vector('list', length = length(po_list))

for (i in 1:length(synth_LIST)) {
  
  par_mat <- parmat(params = coef(po_list[[i]]))
  sim_determ <- trajectory(po_list[[i]], params = par_mat, format = 'array')
  
  synth_list_TEMP <- vector('list', length = n)
  
  for (j in 1:n) {
    
    out_temp <- rmeasure(object = po_list[[i]], x = sim_determ,
                         time = time(po_list[[i]]), params = par_mat) %>%
      matrix(nrow = 2, byrow = FALSE)
    
    rownames(out_temp) <- c('n_P1', 'n_P2')
    out_temp[is.nan(out_temp)] <- NA
    
    synth_list_TEMP[[j]] <- out_temp
    rm(out_temp)
    
  }
  
  synth_LIST[[i]] <- synth_list_TEMP
  rm(synth_list_TEMP, sim_determ, par_mat)
  
}
rm(i, j)

# Fit model to synthetic data:

# During the true analysis, starting parameter sets are chosen from the range of parameters
# fit during "round 2." Specifically, ten starting parameter sets are chosen, and the model
# is fit to each of the 500 synthetic parameter sets. The code below shows how to set this
# up using the first of 500 synthetic parameter sets and the same starting parameter set
# used above. As the fitting procedure itself is the same as performed above for "round 2,"
# the actual fitting step below is commented out, in the interest of time. In the true
# analysis, the best of the 10 fits for each synthetic dataset is chosen, and confidence
# intervals are obtained by calculating the highest posterior density for each parameter.

for (i in 1:length(po_list)) {
  po_list[[i]]@data <- synth_LIST[[i]][[1]]
}

x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
x0_trans_names <- names(x0_trans)

# m <- try(
#   nloptr(x0 = x0_trans, 
#          eval_f = calculate_global_loglik,
#          opts = list(algorithm = "NLOPT_LN_SBPLX",
#                      maxtime = 60 * 60,
#                      maxeval = -1, # Negative value: criterion is disabled
#                      xtol_rel = -1, # Default value: 1e-4
#                      print_level = 0))
# )
# 
# if (!inherits(m, 'try-error')) {
#   x0_fit <- m$solution
#   names(x0_fit) <- x0_trans_names
#   x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], seasons, estpars, shared_estpars)
#   
#   out <- list(estpars = x0_fit_untrans,
#               ll = -m$objective,
#               message = m$message)
# }

# ---------------------------------------------------------------------------------------------------------------------

# Wrapping up

# Get total demo run time:
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'mins'
print(etime)

# Clean up:
rm(list = ls())

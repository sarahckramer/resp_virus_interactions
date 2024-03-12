# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)

# Set vaccination coverage levels and time points:
vacc_cov_vec <- round(seq(0.05, 1.0, by = 0.05), digits = 2) # seq(0.1, 1.0, by = 0.1)
vacc_time_vec <- round(seq(0, 52, by = 1)) # seq(0, 52, by = 2)

# Set parameters for run:
vir2 <- 'rsv'

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

debug_bool <- FALSE

# Choose season and vaccine coverage level:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
# jobid <- 1

p_vacc <- vacc_cov_vec[(jobid - 1) %% 20 + 1]; print(p_vacc)

# Which assumptions about vaccine efficacy/duration are made?
sens_sim <- as.character(Sys.getenv("SENS")); print(sens_sim)

# Set vaccination efficacy against flu:
if (sens_sim == 'vacceff60') {
  vacc_eff <- 0.6
} else if (sens_sim == 'vacceff95') {
  vacc_eff <- 0.95
} else {
  vacc_eff <- 0.8
}
print(vacc_eff)

# Set interaction parameter names:
int_params <- c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')

# Read in MLEs:
mle_hk <- read_rds('results/MLEs_hk.rds')[1, ]
mle_can <- read_rds('results/MLEs_canada.rds')[1, ]

# ---------------------------------------------------------------------------------------------------------------------

# Get values for impact of NATURAL influenza infection:
int_param_vals <- mle_hk %>%
  select(all_of(int_params))

if (sens_sim == 'fit_can') {
  int_param_vals <- mle_can %>%
    select(all_of(int_params))
}

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code (Hong Kong / "subtropical" scenario)

# Set parameters for run:
vir1 <- 'flu_h1_plus_b'
sens <- 'main'
fit_canada <- FALSE
fit_us <- FALSE

# Get year:
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
yr <- seasons[ceiling(jobid / 20)]; print(yr)

# Specify shared parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')

# Get MLEs for each season:
mle <- mle_hk

# Perform model checks:
source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')

# Set desired model parameter values:
model_params <- mle %>%
  dplyr::select(all_of(shared_estpars), contains(yr)) %>%
  rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
  unlist()

if (sens_sim == 'vacc_can') {
  temp_eff <- mle_can %>% select('theta_lambda1') %>% unlist()
  model_params <- c(model_params, unname(temp_eff), unname(model_params['delta1']), vacc_eff)
  rm(temp_eff)
} else if (sens_sim == 'deltavacc1month') {
  model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 30, vacc_eff)
} else if (sens_sim == 'deltavacc6months') {
  model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
} else {
  model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
}

names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')

resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                      Ri1_max = Ri_max1,
                                      Ri2_max = Ri_max2,
                                      d2_max = d2_max,
                                      t0_eff = 0,
                                      debug_bool = debug_bool,
                                      sens = sens,
                                      test_diff = FALSE)

resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)

model_params <- parmat(params = coef(resp_mod), nrep = 2)

# Also run where 1) RSV has no impact on flu, and 2) no flu is circulating:
# model_params['theta_lambda2', ] <- 1.0
# model_params['I10', ] <- 0

# Set vaccination coverage:
model_params['p_vacc', ] <- c(0, p_vacc)

# Initiate results data frame:
res <- NULL

# Loop through vaccination time points:
for (t_vacc in vacc_time_vec) {
  
  # Run deterministic model:
  sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = FALSE) %>%
    dplyr::select(time, H1:H2, .id) %>%
    mutate(vacc_cov = p_vacc,
           vacc_time = t_vacc,
           season = yr)
  
  res <- res %>% bind_rows(sim_temp)
  
}

# Write simulation to file:
write_rds(res, paste0('results/vaccination_simulation_study/simulations/sim_determ_', yr, '_', p_vacc * 100, 'perc_SUBTROPICAL.rds'))

# # Check that, if p_vacc = 0 (no vaccination), all vaccine timepoints yield same results:
# res_comp1 <- res %>% filter(.id == 1, vacc_time == min(vacc_time))
# for (t in unique(res$vacc_time)[-which.min(res$vacc_time)]) {
#   res_comp2 <- res %>% filter(.id == 1 & vacc_time == t)
#   print(all.equal(res_comp1$H1, res_comp2$H1))
#   print(all.equal(res_comp1$H2, res_comp2$H2))
# }
# rm(t, res_comp1, res_comp2)

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code (Canada / "temperate" scenario)

# Set parameters for run:
vir1 <- 'flu'
sens <- 'sinusoidal_forcing'
fit_canada <- TRUE

# Get year:
seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
yr <- seasons[ceiling(jobid / 20)]; print(yr)

# Specify shared parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')

# Check whether all seasons already run:
if (!is.na(yr)) {
  
  # Get MLEs for each season:
  mle <- mle_can
  
  # Perform model checks:
  source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')
  
  # Set desired model parameter values:
  model_params <- mle %>%
    dplyr::select(all_of(shared_estpars), contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  model_params[int_params] <- int_param_vals %>% unlist()
  
  if (sens_sim == 'vacc_can') {
    temp_eff <- mle_can %>% select('theta_lambda1') %>% unlist()
    model_params <- c(model_params, unname(temp_eff), unname(model_params['delta1']), vacc_eff)
    rm(temp_eff)
  } else if (sens_sim == 'deltavacc1month') {
    model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 30, vacc_eff)
  } else if (sens_sim == 'deltavacc6months') {
    model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
  } else {
    model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
  }
  
  names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')
  
  resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                        Ri1_max = Ri_max1,
                                        Ri2_max = Ri_max2,
                                        d2_max = d2_max,
                                        t0_eff = 0,
                                        debug_bool = debug_bool,
                                        sens = sens,
                                        test_diff = TRUE)
  
  resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)
  
  model_params <- parmat(params = coef(resp_mod), nrep = 2)
  
  # Also run where 1) RSV has no impact on flu, and 2) no flu is circulating:
  # model_params['theta_lambda2', ] <- 1.0
  # model_params['I10', ] <- 0
  
  # Set vaccination coverage:
  model_params['p_vacc', ] <- c(0, p_vacc)
  
  # Initiate results data frame:
  res <- NULL
  
  # Loop through vaccination time points:
  for (t_vacc in vacc_time_vec) {
    
    # Run deterministic model:
    sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = TRUE) %>%
      dplyr::select(time, H1:H2, .id) %>%
      mutate(vacc_cov = p_vacc,
             vacc_time = t_vacc,
             season = yr)
    
    res <- res %>% bind_rows(sim_temp)
    
  }
  
  # Write simulation to file:
  write_rds(res, paste0('results/vaccination_simulation_study/simulations/sim_determ_', yr, '_', p_vacc * 100, 'perc_TEMPERATE.rds'))
  
}

# Print message when completed:
print('Done!')

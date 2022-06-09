# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)

# Set vaccination efficacy against flu:
vacc_eff <- 0.6

# Set parameters for run:
vir1 <- 'flu_h1'
vir2 <- 'rsv'
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

debug_bool <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code

# Get MLEs for each season:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

# Create list to store results of simulations:
res_list <- vector('list', length(seasons))

# Loop through seasons and run:
for (yr_index in 1:length(seasons)) {
  
  # Set current season:
  yr <- seasons[yr_index]
  
  # Perform model checks:
  source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')
  
  # Set desired model parameter values:
  model_params <- mle %>%
    dplyr::select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff, 0)
  names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff', 'p_vacc')
  
  # Initiate results data frame:
  res <- NULL
  
  # Loop through vaccine coverage levels:
  for (p_vacc in seq(0, 1.0, by = 0.05)) {
    
    # Update vaccination coverage:
    model_params['p_vacc'] <- p_vacc
    
    # Loop through vaccination time points:
    for (t_vacc in seq(1, 52, by = 1)) {
      
      # Run deterministic model:
      sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool) %>%
        dplyr::select(time, H1:H2) %>%
        mutate(vacc_cov = p_vacc,
               vacc_time = t_vacc)
      
      res <- res %>% bind_rows(sim_temp)
      
    }
    
  }
  
  # Store results for current season:
  res_list[[yr_index]] <- res
  
}

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

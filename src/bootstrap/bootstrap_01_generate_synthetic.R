# ---------------------------------------------------------------------------------------------------------------------
# Generate synthetic data for use in bootstrapping
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(9489703)

# Load libraries:
library(tidyverse)

# Set key parameters:
n <- 500 # How many synthetic datasets to create?
sens <- 'main'
fit_canada <- FALSE
fit_us <- FALSE

if (fit_us) {
  region_num <- 8
  region <- paste0('Region ', region_num)
}

if (fit_canada | fit_us) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

# ---------------------------------------------------------------------------------------------------------------------

# Set up pomp models with MLE of model parameters

# Get MLE:
if (sens != 'main') {
  
  if (fit_canada) {
    mle <- read_rds(paste0('results/round2_fit/sens/canada/MLEs_flu_h1_plus_b.rds'))[1, ]
  } else if (fit_us) {
    mle <- read_rds(paste0('results/round2_fit/sens/us/region_', region_num, '/MLEs_flu.rds'))[1, ]
  } else {
    mle <- read_rds(paste0('results/round2_fit/sens/', sens, '/MLEs_flu_h1_plus_b.rds'))[1, ]
  }
  
} else {
  mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')[1, ]
}

# Get names of shared parameters:
if (sens == 'sinusoidal_forcing') {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
  if (fit_us) {
    shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                        'b1', 'b2', 'phi1', 'phi2')
  }
} else if (sens == 'h3_covar') {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2',
                      'beta_h3')
} else if (sens == 'no_ah') {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2')
} else if (sens == 'no_int') {
  shared_estpars <- c('rho1', 'rho2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
} else {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
}

# Set model parameters:
debug_bool <- FALSE
vir2 <- 'rsv'
prof_lik <- FALSE

seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
if (sens == 'less_circ_h3') {
  seasons <- c('s17-18', 's18-19')
}
if (fit_canada) {
  seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
}
if (fit_us) {
  seasons <- c('s10-11', 's11-12', 's12-13', 's13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
}

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

# Load models for each season:
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

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# For each season, set parameter values to MLE:
for (i in 1:length(seasons)) {
  
  season <- seasons[i]
  
  mle_temp <- mle %>%
    select(all_of(shared_estpars), contains(season)) %>%
    rename_with(~ str_remove(.x, paste0(season, '_')))
  
  coef(po_list[[i]], names(mle_temp)) <- mle_temp
  
  rm(mle_temp)
}
rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Generate synthetic data

# Draw number of observed, lab-confirmed cases (stochastic) from deterministic simulations:
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

# Plot synthetic data vs. observed data:
to_plot <- sample(1:n, size = 6)
for (i in 1:length(seasons)) {
  
  obs_temp <- po_list[[i]]@data %>%
    t() %>%
    as_tibble() %>%
    mutate(time = 1:ncol(po_list[[i]]@data)) %>%
    pivot_longer(cols = -time, names_to = 'virus', values_to = 'val') %>%
    mutate(virus = if_else(virus == 'n_P1', 'Flu', 'RSV'))
  
  synth_temp <- NULL
  
  for (j in to_plot) {
    
    synth_temp <- rbind(synth_temp,
                        synth_LIST[[i]][[j]] %>%
                          t() %>%
                          as_tibble() %>%
                          mutate(time = 1:ncol(po_list[[i]]@data),
                                 sim = j) %>%
                          pivot_longer(cols = -c(time, sim), names_to = 'virus', values_to = 'val') %>%
                          mutate(virus = if_else(virus == 'n_P1', 'Flu', 'RSV'))
    )
    
  }
  
  p_temp <- ggplot() + geom_point(data = obs_temp, aes(x = time, y = val, col = virus)) +
    geom_line(data = synth_temp, aes(x = time, y = val, col = virus)) +
    facet_wrap(~ sim) +
    theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time', y = '# of Observed Cases', col = 'Virus', title = seasons[i])
  print(p_temp)
  
}

# Save synthetic data:
if (sens != 'main') {
  
  if (fit_canada) {
    
    write_rds(synth_LIST, 'results/round2_fit/sens/canada/synth_data_for_bootstrapping_flu.rds')
    
  } else if (fit_us) {
    
    write_rds(synth_LIST, paste0('results/round2_fit/sens/us/region_', region_num, '/synth_data_for_bootstrapping_flu.rds'))
    
  } else {
    
    write_rds(synth_LIST, paste0('results/round2_fit/sens/', sens, '/synth_data_for_bootstrapping_', vir1, '.rds'))
    
  }
  
} else {
  
  write_rds(synth_LIST, paste0('results/synth_data_for_bootstrapping_', vir1, '.rds'))
  
}

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Code to fit flu/RSV interaction model using panelPomp
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(panelPomp)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
int_eff <- as.character(Sys.getenv("INTERACTIONEFFECT")); print(int_eff)
vir1 <- as.character(Sys.getenv("VIRUS1")); print(vir1)
prof_lik <- as.logical(Sys.getenv("PROFLIK")); print(prof_lik)

mif_particles_number <- as.integer(Sys.getenv("MIFPARTICLES")); print(mif_particles_number)
pf_particles_number <- as.integer(Sys.getenv("PFPARTICLES")); print(pf_particles_number)
filtering_number <- as.integer(Sys.getenv("FILTERING")); print(filtering_number)

# # Set parameters for local run:
# jobid <- 1
# no_jobs <- 20
# vir1 <- 'flu_h1' # 'flu_A', 'flu_B'
# 
# sobol_size <- 100
# search_type <- 'round2_CIs'
# int_eff <- 'susc' # 'susc' or 'sev' - fit impact of interaction on susceptibility or severity?
# prof_lik <- FALSE
# 
# mif_particles_number <- 400
# pf_particles_number <- 1000
# filtering_number <- 20

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

lag_val <- 0

# if (prof_lik) {
#   prof_param <- 'theta_lambda1'
#   # prof_param <- 'theta_lambda2'
#   # prof_param <- 'delta1'
#   # prof_param <- 'd2'
#   
#   if (prof_param == 'delta1') {
#     prof_val <- (7 / seq(5, 255, by = 5))[jobid]
#   } else if (prof_param == 'd2') {
#     prof_val <- c(0.01, seq(0.1, 0.9, by = 0.1), seq(1, 5, by = 0.1))[jobid]
#   }
#   else {
#     prof_val <- seq(0.0, 1.0, by = 0.02)[jobid]
#   }
#   print(prof_val)
#   
#   jobid_orig <- jobid
#   jobid <- 1
#   
# } else {
#   jobid_orig <- jobid
# }

# ---------------------------------------------------------------------------------------------------------------------

# Create panelPomp object

# Load pomp objects for each season:
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

# Name list of pomps:
names(po_list) <- seasons

# # Explore impact of sigmaSE on simulations:
# resp_mod <- po_list[[1]]
# sigma_val <- 0.1
# coef(resp_mod)[c('Ri1', 'Ri2', 'sigmaSE', 'R120')] <- c(1.25, 1.5, sigma_val, 0.5)
# simulate(resp_mod, nsim = 20, format = 'data.frame') %>%
#   select(time:.id, n_P1:n_P2) %>%
#   pivot_longer(n_P1:n_P2, names_to = 'vir', values_to = 'cases') %>%
#   ggplot(aes(x = time, y = cases, col = vir, group = paste(vir, .id))) + geom_line() + theme_classic() + labs(title = sigma_val)

# Construct panelPomp object:
shared_param_vals <- c(7 / 5, 1.0, 1.0, 1.0, 0.5, 0.15, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7 / 5, 7 / 10)
names(shared_param_vals) <- c('delta1', 'd2', 'theta_lambda1', 'theta_lambda2', 'rho1', 'rho2', 'alpha', 'phi',
                              'theta_rho1', 'theta_rho2', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2',
                              'beta_sd1', 'beta_sd2', 'gamma1', 'gamma2')

resp_mod <- panelPomp(po_list,
                      shared = shared_param_vals)

# ---------------------------------------------------------------------------------------------------------------------

# Get estpars and start values

# Choose parameters to estimate:
if (int_eff == 'susc') {
  if (prof_lik) {
    estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
    estpars <- estpars[estpars != prof_param]
  } else {
    estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
    # estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
    #                     'alpha', 'phi', 'eta_temp1', 'eta_temp2')
  }
} else if (int_eff == 'sev') {
  if (prof_lik) {
    estpars <- c('rho1', 'rho2', 'theta_rho1', 'theta_rho2', 'delta1', 'd2',
                 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
    estpars <- estpars[estpars != prof_param]
  } else {
    estpars <- c('rho1', 'rho2', 'theta_rho1', 'theta_rho2', 'delta1', 'd2',
                 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
  }
} else {
  stop('Unrecognized int_eff value.')
}

shared_estpars <- estpars
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
for (i in 1:length(seasons)) {
  estpars <- c(estpars, paste0(unit_estpars, '[', seasons[i], ']'))
}
rm(i)

# Set upper/lower values for global params:
start_range <- data.frame(rho1 = c(0, 1.0),
                          rho2 = c(0, 1.0),
                          theta_lambda1 = c(0, 1.0),
                          theta_lambda2 = c(0, 1.0),
                          theta_rho1 = c(0, 1.0),
                          theta_rho2 = c(0, 1.0),
                          delta1 = c(7 / 60, 7),
                          d2 = c(0, 10),
                          alpha = c(0, 0.5),
                          phi = c(0, 52.25),
                          eta_temp1 = c(-0.5, 0.5),
                          eta_temp2 = c(-0.5, 0.5),
                          eta_ah1 = c(-0.5, 0.5),
                          eta_ah2 = c(-0.5, 0.5),
                          beta_sd1 = c(0, 0.5),
                          beta_sd2 = c(0., 0.5))

# Set upper/lower values for unit params (broad):
unit_start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                               Ri2 = c(1.0, Ri_max2),
                               I10 = c(0, 1e-3),
                               I20 = c(0, 1e-3),
                               R10 = c(0, 0.3),
                               R20 = c(0, 0.3),
                               R120 = c(0, 0.3))

# Get data frame of all ranges:
if (search_type == 'round2_CIs') {
  
  if (vir1 == 'flu_h1') {
    
    if (int_eff == 'susc') {
      start_range <- read_rds('results/round2_cis/round2CI_startvals_PROF_H1.rds')
    } else if (int_eff == 'sev') {
      stop('SEV not yet implemented!')
    } else {
      stop('Unrecognized int_eff!')
    }
    
    for (season in seasons) {
      names(start_range)[str_detect(names(start_range), season)] <- paste0(str_remove(names(start_range)[str_detect(names(start_range), season)],
                                                                                      paste0(season, '_')), '[', season, ']')
    }
    start_range <- cbind(start_range, beta_sd1 = c(0.01, 0.5), beta_sd2 = c(0.01, 0.5))
    
  } else if (vir1 == 'flu_b') {
    
    if (int_eff == 'susc') {
      start_range <- read_rds('results/round2_cis/round2CI_startvals_PROF_B.rds')
    } else if (int_eff == 'sev') {
      stop('SEV not yet implemented!')
    } else {
      stop('Unrecognized int_eff!')
    }
    
    for (season in seasons) {
      names(start_range)[str_detect(names(start_range), season)] <- paste0(str_remove(names(start_range)[str_detect(names(start_range), season)],
                                                                                      paste0(season, '_')), '[', season, ']')
    }
    start_range <- cbind(start_range, beta_sd1 = c(0.01, 0.5), beta_sd2 = c(0.01, 0.5))
    
  } else {
    stop('Unknown vir1!')
  }
  
} else {
  
  for (i in 1:length(seasons)) {
    
    if (search_type == 'broad') {
      start_range_temp <- unit_start_range
    } else if (search_type == 'round2_CIs') {
      print('ERROR: Round2 CIs used.')
    } else {
      stop('Unrecognized search type!')
    }
    
    names(start_range_temp) <- paste0(names(unit_start_range), '[', seasons[i], ']')
    start_range <- start_range %>%
      bind_cols(start_range_temp)
    rm(start_range_temp)
    
  }
  rm(i)
  
}

start_range <- start_range[, estpars]

# Get starting values for each parameter:
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

print(estpars)
print(start_range)
print(summary(start_values))

# ---------------------------------------------------------------------------------------------------------------------

# Run model fitting

# Get unique identifiers:
sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)

# Loop through start values and perform MIF:
for (i in seq_along(sub_start)) {
  
  print(paste0('Estimation: ', sub_start[i]))
  
  # Get param start values:
  x0 <- as.numeric(start_values[sub_start[i], ])
  coef(resp_mod)[estpars] <- x0
  params <- coef(resp_mod)
  print(params)
  
  # Run MIF:
  tic <- Sys.time()
  
  if (int_eff == 'susc') {

    mf <- try (

      mif2(resp_mod,
           Np = mif_particles_number,
           Nmif = filtering_number,
           cooling.type = 'geometric',
           cooling.fraction.50 = 0.5,
           rw.sd = rw.sd(rho1 = 0.02,
                         rho2 = 0.02,
                         delta1 = 0.02,
                         d2 = 0.02,
                         theta_lambda1 = 0.02,
                         theta_lambda2 = 0.02,
                         alpha = 0.02,
                         phi = 0.02,
                         eta_temp1 = 0.02,
                         eta_temp2 = 0.02,
                         eta_ah1 = 0.02,
                         eta_ah2 = 0.02,
                         # beta_sd1 = 0.02,
                         # beta_sd2 = 0.02,
                         Ri1 = 0.02,
                         Ri2 = 0.02,
                         I10 = ivp(0.1),
                         I20 = ivp(0.1),
                         R10 = ivp(0.1),
                         R20 = ivp(0.1),
                         R120 = ivp(0.1))) %>%
        mif2(Nmif = filtering_number,
             cooling.fraction.50 = 0.25)

    )

  } else if (int_eff == 'sev') {

    mf <- try (

      mif2(resp_mod,
           Np = mif_particles_number,
           Nmif = filtering_number,
           cooling.type = 'geometric',
           cooling.fraction.50 = 0.5,
           rw.sd = rw.sd(rho1 = 0.02,
                         rho2 = 0.02,
                         delta1 = 0.02,
                         d2 = 0.02,
                         theta_rho1 = 0.02,
                         theta_rho2 = 0.02,
                         alpha = 0.02,
                         phi = 0.02,
                         eta_temp1 = 0.02,
                         eta_temp2 = 0.02,
                         eta_ah1 = 0.02,
                         eta_ah2 = 0.02,
                         # beta_sd1 = 0.02,
                         # beta_sd1 = 0.02,
                         Ri1 = 0.02,
                         Ri2 = 0.02,
                         I10 = ivp(0.1),
                         I20 = ivp(0.1),
                         R10 = ivp(0.1),
                         R20 = ivp(0.1),
                         R120 = ivp(0.1))) %>%
        mif2(Nmif = filtering_number,
             cooling.fraction.50 = 0.25)

    )

  } else {
    stop('Unrecognized int_eff!')
  }

  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)

  # If estimation is successful, save results:
  if (!inherits(mf, 'try-error')) {
    tic <- Sys.time()

    # ll <- replicate(10, mf %>% pfilter(Np = pf_particles_number) %>% logLik()) %>%
    #   logmeanexp(se = TRUE)
    ll <- replicate(10, mf %>% pfilter(Np = pf_particles_number) %>% unitlogLik()) %>%
      panel_logmeanexp(MARGIN = 1, se = TRUE)
    
    print(ll)

    toc <- Sys.time()
    etime <- toc - tic
    units(etime) <- 'mins'
    print(etime)

    mf.ll <- mf %>%
      coef() %>%
      bind_rows() %>%
      bind_cols(loglik = ll[1],
                loglik.se = ll[2]) %>%
      select(all_of(estpars),
             loglik, loglik.se)

    out <- list(mf = mf,
                # res = mf.ll,
                Nf = filtering_number,
                Np_mif = mif_particles_number,
                Np_pf = pf_particles_number)

    # Write to file:
    saveRDS(out,
            file = sprintf('results/mif_res_%s_%s_%s.rds',
                           vir1,
                           int_eff,
                           sub_start[i])
    )

    # Print results:
    print(out$res %>%
            select(shared_estpars,
                   loglik,
                   loglik.se))

  }
  
}

print('Done!')

# # Plot traces:
# out_df <- NULL
# for (i in 1:length(out)) {
#   df_temp <- out[[i]]$mf %>%
#     traces() %>%
#     as_tibble() %>%
#     rownames_to_column(var = 'iteration') %>%
#     mutate(run = as.character(i),
#            iteration = as.numeric(iteration)) %>%
#     select(run, iteration, loglik, all_of(estpars)) %>%
#     pivot_longer(-c(run, iteration), names_to = 'param', values_to = 'value')
#   out_df <- rbind(out_df, df_temp)
# }
# 
# p1 <- ggplot(data = out_df, aes(x = iteration, y = value, group = run, col = run)) +
#   geom_line() + facet_wrap(~ param, scales = 'free') +
#   theme_classic() + scale_color_brewer(palette = 'Set1')
# print(p1)

# ---------------------------------------------------------------------------------------------------------------------

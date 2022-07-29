# ---------------------------------------------------------------------------------------------------------------------
# Code to explore possible reasons for differences in vaccine impact on RSV by season
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(patchwork)

# Set vaccination coverage levels and time points:
i0_vec <- c(1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 2.5e-5, 5e-5, 7.5e-5, 1e-4)#, 2.5e-4, 5e-4)
props <- c(1.0, 0.8, 0.6, 0.4, 0.2, 0)
vacc_cov_vec <- seq(0.1, 0.6, by = 0.1)
vacc_time_vec <- seq(0, 32, by = 4)

# Set vaccination efficacy against flu:
vacc_eff <- 0.6

# Set parameters for run:
vir1 <- 'flu_h1'
vir2 <- 'rsv'
seasons <- c('s13-14', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

debug_bool <- FALSE

# ---------------------------------------------------------------------------------------------------------------------

# Run small simulation study for each season, varying I10 to alter phase difference/early flu AR

# Get MLEs for each season:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

# Loop through seasons and run:
pdf('results/vaccine_simulation_study/plots/sims_change_I10.pdf',
    width = 20, height = 10)
for (yr_index in 1:length(seasons)) {
  
  # Set current season:
  yr <- seasons[yr_index]
  print(yr)
  
  # Perform model checks:
  source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')
  
  # Set desired model parameter values:
  model_params <- mle %>%
    dplyr::select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
  names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')
  
  resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                        Ri1_max = Ri_max1,
                                        Ri2_max = Ri_max2,
                                        d2_max = d2_max,
                                        t0_eff = 0,
                                        debug_bool = debug_bool)
  
  resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)
  
  # Loop through I10 values:
  p_list <- vector('list', length(i0_vec))
  for (i0_val in i0_vec) {
    
    # Print progress:
    print(i0_val)
    
    # Update I10:
    coef(resp_mod, 'I10') <- i0_val
    
    # Run deterministic simulation:
    sim_temp <- trajectory(resp_mod, format = 'data.frame') %>%
      select(H1:time) %>%
      pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val')
    
    # Plot to check effect on phase/early AR
    p_temp <- ggplot(data = sim_temp, aes(x = time, y = val, group = Virus, col = Virus)) +
      geom_line() + theme_classic() + labs(title = i0_val)
    p_list[[which(i0_vec == i0_val)]] <- p_temp
    
    # Calculate difference in peak timing, AR in first few weeks, in absence of vaccination:
    flu_metrics_temp <- sim_temp %>%
      pivot_wider(names_from = Virus, values_from = val) %>%
      summarise(pt1 = which.max(H1), pt2 = which.max(H2),
                ar5 = sum(H1[time <= 5]), ar10 = sum(H1[time <= 10])) %>%
      mutate(pt_diff = pt1 - pt2, i10 = i0_val) %>%
      select(i10, pt_diff, ar5:ar10)
    
    # Get matrix of model parameters for simulation study:
    model_params <- parmat(params = coef(resp_mod), nrep = 2)
    
    # Loop through vaccine coverage levels:
    for (p_vacc in vacc_cov_vec) {
      
      # Print progress:
      print(p_vacc)
      
      # Initiate results data frame:
      res <- NULL
      
      # Update vaccination coverage:
      model_params['p_vacc', ] <- c(0, p_vacc)
      
      # Loop through vaccination time points:
      for (t_vacc in vacc_time_vec) {
        
        # Run deterministic model:
        sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool) %>%
          dplyr::select(time, H1:H2, .id) %>%
          mutate(vacc_cov = p_vacc,
                 vacc_time = t_vacc,
                 season = yr)
        
        res <- res %>% bind_rows(sim_temp)
        
      }
      
      # Join flu metrics to res:
      res <- res %>%
        bind_cols(flu_metrics_temp)
      
      # Write simulation to file:
      write_rds(res, paste0('results/vaccine_simulation_study/simulations/sim_determ_', yr, '_', i0_val, '_', p_vacc * 100, 'perc.rds'))
      
    }
    
    # Clean up:
    rm(res)
    
  }
  
  # Output plots of outbreaks by I10, without vaccination:
  do.call('grid.arrange', p_list)
  
}
dev.off()

# Read in all simulation results:
file_list <- list.files('results/vaccine_simulation_study/simulations/SA_vary_phase_and_early_AR_2/',
                        pattern = 'sim', full.names = TRUE)

res_list <- vector('list', length(file_list))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_sub(file_list[i], 86, 91))
}
rm(i)

# Combine into single tibble:
res <- bind_rows(res_list) %>%
  as_tibble()

# Calculate measures of vaccine impact from simulation results:
res_metrics <- res %>%
  group_by(season, vacc_cov, vacc_time, .id, i10) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2),
            pt_diff = unique(pt_diff), ar5 = unique(ar5), ar10 = unique(ar10)) %>%
  ungroup() %>%
  group_by(season, vacc_cov, vacc_time, i10) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1],
            pt_diff = unique(pt_diff), ar5 = unique(ar5), ar10 = unique(ar10)) %>%
  ungroup()

# Loop through seasons and initial conditions and plot results:
res_metrics <- res_metrics %>%
  pivot_longer(ar1_impact:ar2_impact, names_to = 'metric', values_to = 'val')

facet_labs <- c('Influenza', 'RSV')
names(facet_labs) <- c('ar1_impact', 'ar2_impact')

pdf('results/vaccine_simulation_study/plots/vacc_impact_AR_SA1.pdf',
    width = 20, height = 10)
for (seas in seasons) {
  
  res_temp <- res_metrics %>% filter(season == seas)
  upper_bound_ar <- ceiling(max(res_temp$val) * 10) / 10
  
  p_list <- vector('list', length(i0_vec))
  
  for (i0_val in i0_vec) {
    
    res_temp_by_i0 <- res_temp %>% filter(i10 == i0_val)
    pt_diff_val <- unique(res_temp_by_i0$pt_diff)
    ar5_val <- signif(unique(res_temp_by_i0$ar5), 3)
    ar10_val <- signif(unique(res_temp_by_i0$ar10), 3)
    
    p_temp <- ggplot(data = res_temp_by_i0, aes(x = vacc_time, y = vacc_cov, fill = val)) +
      geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
      theme_classic() +
      scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                           breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
           title = paste0(seas, ' / ', i0_val, ' (pt_diff=', pt_diff_val, ', ar5=', ar5_val, ')'), fill = 'Impact')
    p_list[[which(i0_vec == i0_val)]] <- p_temp
    
  }
  
  do.call('grid.arrange', p_list)
  
}
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Run small simulation study for each season, varying the proportion of recovered allocated to R10/R20 vs. R120

# Loop through seasons and run:
pdf('results/vaccine_simulation_study/plots/sims_change_R_allocation.pdf',
    width = 20, height = 10)
for (yr_index in 1:length(seasons)) {
  
  # Set current season:
  yr <- seasons[yr_index]
  print(yr)
  
  # Perform model checks:
  source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')
  
  # Set desired model parameter values:
  model_params <- mle %>%
    dplyr::select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
  names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')
  
  resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                        Ri1_max = Ri_max1,
                                        Ri2_max = Ri_max2,
                                        d2_max = d2_max,
                                        t0_eff = 0,
                                        debug_bool = debug_bool)
  
  resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)
  
  # Calculate maximum possible value of R120, if we want to hold total proportion immune to flu/RSV constant:
  max_r120 <- min(c(model_params['R10'] + model_params['R120'], model_params['R20'] + model_params['R120']))
  
  # Calculate desired total immune to flu and RSV:
  total_R_flu <- model_params['R10'] + model_params['R120']
  total_R_rsv <- model_params['R20'] + model_params['R120']
  
  # Loop and allocate some proportion of max_r120 to R120 vs. R10/R20:
  p_list <- vector('list', length(props))
  
  for (prop_allocate in props) {
    
    # Print progress:
    print(prop_allocate)
    
    # Update initial conditions:
    coef(resp_mod, c('R10', 'R20', 'R120')) <- c(total_R_flu - prop_allocate * max_r120, total_R_rsv - prop_allocate * max_r120, prop_allocate * max_r120)
    
    # Check that total immune to flu and RSV are correct:
    expect_true(round(sum(coef(resp_mod, c('R10', 'R120'))), 7) == round(total_R_flu, 7))
    expect_true(sum(coef(resp_mod, c('R20', 'R120'))) == total_R_rsv)
    
    # Check that initial conditions still sum to less than 1.0:
    if (sum(coef(resp_mod, c('I10', 'I20', 'R10', 'R20', 'R120'))) < 1.0) {
      
      # Run deterministic simulation:
      sim_temp <- trajectory(resp_mod, format = 'data.frame') %>%
        select(H1:time) %>%
        pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val')
      
      # Plot to check impact on outbreak trajectory:
      p_temp <- ggplot(data = sim_temp, aes(x = time, y = val, group = Virus, col = Virus)) +
        geom_line() + theme_classic() + labs(title = prop_allocate)
      p_list[[which(props == prop_allocate)]] <- p_temp
      
      # Get matrix of model parameters for simulation study:
      model_params <- parmat(params = coef(resp_mod), nrep = 2)

      # Loop through vaccine coverage levels:
      for (p_vacc in vacc_cov_vec) {

        # Print progress:
        print(p_vacc)

        # Initiate results data frame:
        res <- NULL

        # Update vaccination coverage:
        model_params['p_vacc', ] <- c(0, p_vacc)

        # Loop through vaccination time points:
        for (t_vacc in vacc_time_vec) {

          # Run deterministic model:
          sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool) %>%
            dplyr::select(time, H1:H2, .id) %>%
            mutate(vacc_cov = p_vacc,
                   vacc_time = t_vacc,
                   season = yr)

          res <- res %>% bind_rows(sim_temp)

        }

        # Join flu metrics to res:
        res <- res %>%
          mutate(prop = prop_allocate)

        # Write simulation to file:
        write_rds(res, paste0('results/vaccine_simulation_study/simulations/sim_determ_', yr, '_', prop_allocate, '_', p_vacc * 100, 'perc.rds'))

      }

      # Clean up:
      rm(res)

    }

  }
  
  # Remove empty elements of plot list:
  p_list <- p_list[lapply(p_list, length) > 0]
  
  # Output plots of outbreaks by prop_allocate, without vaccination:
  do.call('grid.arrange', p_list)
  
}
dev.off()

# Read in all simulation results:
file_list <- list.files('results/vaccine_simulation_study/simulations/SA_vary_prop_allocated_R120/',
                        pattern = 'sim', full.names = TRUE)

res_list <- vector('list', length(file_list))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_sub(file_list[i], 85, 90))
}
rm(i)

# Combine into single tibble:
res <- bind_rows(res_list) %>%
  as_tibble()

# Calculate measures of vaccine impact from simulation results:
res_metrics <- res %>%
  group_by(season, vacc_cov, vacc_time, .id, prop) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(season, vacc_cov, vacc_time, prop) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

# Loop through seasons and initial conditions and plot results:
res_metrics <- res_metrics %>%
  pivot_longer(ar1_impact:ar2_impact, names_to = 'metric', values_to = 'val')

facet_labs <- c('Influenza', 'RSV')
names(facet_labs) <- c('ar1_impact', 'ar2_impact')

pdf('results/vaccine_simulation_study/plots/vacc_impact_AR_SA2.pdf',
    width = 20, height = 10)
for (seas in seasons) {
  
  res_temp <- res_metrics %>% filter(season == seas)
  upper_bound_ar <- ceiling(max(res_temp$val) * 10) / 10
  
  p_list <- vector('list', length(props))
  
  for (prop_allocate in props) {
    
    res_temp_by_prop <- res_temp %>% filter(prop == prop_allocate)
    
    if (nrow(res_temp_by_prop) > 0) {
      p_temp <- ggplot(data = res_temp_by_prop, aes(x = vacc_time, y = vacc_cov, fill = val)) +
        geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
        theme_classic() +
        scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                             breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
        scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
        labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
             title = paste0(seas, ' / ', prop_allocate), fill = 'Impact')
      p_list[[which(props == prop_allocate)]] <- p_temp
    }
    
  }
  
  p_list <- p_list[lapply(p_list, length) > 0]
  do.call('grid.arrange', p_list)
  
}
dev.off()

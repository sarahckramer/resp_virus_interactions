# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(patchwork)

# Set vaccination coverage levels and time points:
vacc_cov_vec <- seq(0.1, 1.0, by = 0.1) # seq(0.05, 1.0, by = 0.05)
vacc_time_vec <- seq(0, 52, by = 2) # seq(0, 52, by = 1)

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

# Run main simulation study code (Hong Kong / "subtropical" scenario)

# Get MLEs for each season:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

# Loop through seasons and run:
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
  
  model_params <- parmat(params = coef(resp_mod), nrep = 2)
  
  # Also run where 1) RSV has no impact on flu, and 2) no flu is circulating:
  # model_params['theta_lambda2', ] <- 1.0
  # model_params['I10', ] <- 0
  
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
    
    # Write simulation to file:
    write_rds(res, paste0('results/vaccine_simulation_study/simulations/sim_determ_', yr, '_', p_vacc * 100, 'perc.rds'))
    
    # # Check that, if p_vacc = 0 (no vaccination), all vaccine timepoints yield same results:
    # res_comp1 <- res %>% filter(.id == 1, vacc_time == min(vacc_time))
    # for (t in unique(res$vacc_time)[-which.min(res$vacc_time)]) {
    #   res_comp2 <- res %>% filter(.id == 1 & vacc_time == t)
    #   print(all.equal(res_comp1$H1, res_comp2$H1))
    #   print(all.equal(res_comp1$H2, res_comp2$H2))
    # }
    # rm(t, res_comp1, res_comp2)
    
  }
  
  rm(res)
  
}

# Read in all simulation results:
file_list <- list.files('results/vaccine_simulation_study/simulations/run_initial_hk/', pattern = 'sim', full.names = TRUE)

res_list <- vector('list', length(file_list))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_sub(file_list[i], 72, 77))
}
rm(i)

# Combine into single tibble:
res <- bind_rows(res_list) %>%
  as_tibble()

# Calculate measures of vaccine impact from simulation results:
res_metrics <- res %>%
  group_by(season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

# Want vacc_time relative to timing of flu/RSV peak?
# Get peak timing of flu and RSV without vaccination:
pt_novacc <- res %>%
  filter(.id == 1) %>%
  group_by(season, vacc_cov, vacc_time) %>%
  summarise(pt1 = which.max(H1), pt2 = which.max(H2)) %>%
  ungroup() %>%
  select(!vacc_cov) %>% unique()
expect_true(nrow(pt_novacc) == length(seasons) * length(vacc_time_vec))

# And join these to metrics results:
res_metrics <- res_metrics %>%
  left_join(pt_novacc, by = c('season', 'vacc_time')) %>%
  mutate(rel_vacc_time1 = vacc_time - pt1,
         rel_vacc_time2 = vacc_time - pt2) %>%
  select(!c(pt1:pt2))

# # Write results to file:
# write_rds(res_metrics, 'results/res_MET_simulation_study_TEMP.rds')

# Plot results:
res_metrics <- res_metrics %>%
  pivot_longer(ar1_impact:ar2_impact, names_to = 'metric', values_to = 'val')

res_metrics_AVG <- res_metrics %>%
  group_by(vacc_cov, vacc_time, metric) %>%
  summarise(val = median(val))

upper_bound_ar <- ceiling(max(res_metrics_AVG$val[res_metrics_AVG$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10

facet_labs <- c('Influenza', 'RSV')
names(facet_labs) <- c('ar1_impact', 'ar2_impact')

p_lines_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
  geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', col = 'Vacc. Cov.')

p_heat_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = vacc_cov, fill = val)) +
  geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + #theme(axis.ticks = element_blank()) +
  scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
       fill = 'Impact')

print(p_lines_ar / p_heat_ar)

# Plot results by individual season:
pdf('results/vaccine_simulation_study/plots/vacc_impact_AR.pdf',
    width = 20, height = 10)

for (seas in seasons) {
  
  sim_temp <- res %>%
    filter(season == seas, .id == 1, vacc_cov == 1, vacc_time == 0) %>%
    select(time:H2) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val,
           Virus = if_else(Virus == 'H1', 'Flu', 'RSV'))
  res_temp <- res_metrics %>% filter(season == seas)
  
  p_outbreak <- ggplot(data = sim_temp, aes(x = time, y = val, color = Virus)) +
    geom_line() + geom_point() +
    theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  upper_bound_ar <- ceiling(max(res_temp$val[res_temp$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10
  
  # facet_labs <- c('Influenza', 'RSV')
  # names(facet_labs) <- c('ar1_impact', 'ar2_impact')
  
  p_lines_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
    geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', title = seas, col = 'Vacc. Cov.')
  
  p_heat_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = vacc_cov, fill = val)) +
    geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + #theme(axis.ticks = element_blank()) +
    scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                         breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
         title = seas, fill = 'Impact')
  
  print(p_outbreak | (p_lines_ar / p_heat_ar))
  
}

dev.off()

# Plot all outbreaks by vaccination timing and coverage:
pdf('results/vaccine_simulation_study/plots/all_sims.pdf',
    width = 30, height = 11)
for (seas in seasons) {
  
  sim_temp_all <- res %>%
    filter(season == seas & .id == 2) %>%
    select(-c(.id, season)) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val)
  
  p_outbreak_all <- ggplot(data = sim_temp_all, aes(x = time, y = val, col = Virus)) + geom_line() +
    facet_grid(vacc_cov ~ vacc_time) + theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  print(p_outbreak_all)
  
}
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

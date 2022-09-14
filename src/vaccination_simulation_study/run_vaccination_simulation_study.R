# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)

# Set vaccination coverage levels and time points:
vacc_cov_vec <- round(seq(0.05, 1.0, by = 0.05), digits = 2) # seq(0.1, 1.0, by = 0.1)
vacc_time_vec <- round(seq(0, 52, by = 1)) # seq(0, 52, by = 2)

# Set vaccination efficacy against flu:
vacc_eff <- 0.8

# Set parameters for run:
vir1 <- 'flu_h1'
vir2 <- 'rsv'
seasons <- c('s13-14', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

debug_bool <- FALSE

# Choose season and vaccine coverage level:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
# jobid <- 1

yr <- seasons[ceiling(jobid / 20)]
p_vacc <- vacc_cov_vec[(jobid - 1) %% 20 + 1]

print(yr)
print(p_vacc)

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code (Hong Kong / "subtropical" scenario)

# Get MLEs for each season:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

# Perform model checks:
source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')

# Set desired model parameter values:
model_params <- mle %>%
  dplyr::select(rho1:eta_ah2, contains(yr)) %>%
  rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
  unlist()

model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
# model_params <- c(model_params, 0.5, unname(model_params['delta1']), vacc_eff)
# model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
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

# Set vaccination coverage:
model_params['p_vacc', ] <- c(0, p_vacc)

# Initiate results data frame:
res <- NULL

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

# Run main simulation study code ("temperate" scenario)

# Get MLEs for each season:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

# Get "temperate" parameter values:
temp_params <- read_csv('results/vaccination_simulation_study/temperate_params_to_use.csv')

# Get infection data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_rsv
start_week <- 40

hk_dat <- hk_dat %>%
  mutate(season = NA) %>%
  mutate(season = ifelse((Year == 2014 & Week >= start_week) | (Year == 2015 & Week < start_week), 's14-15', season),
         season = ifelse((Year == 2015 & Week >= start_week) | (Year == 2016 & Week < start_week), 's15-16', season),
         season = ifelse((Year == 2016 & Week >= start_week) | (Year == 2017 & Week < start_week), 's16-17', season),
         season = ifelse((Year == 2017 & Week >= start_week) | (Year == 2018 & Week < start_week), 's17-18', season),
         season = ifelse((Year == 2018 & Week >= start_week) | (Year == 2019 & Week < start_week), 's18-19', season),
         season = ifelse((Year == 2019 & Week >= start_week) | (Year == 2020 & Week < start_week), 's19-20', season),
         season = ifelse(Year == 2014 & Week < start_week, 's13-14', season)) %>%
  filter(!is.na(season)) %>%
  filter(season != 's14-15', season != 's19-20') %>%
  group_by(season) %>%
  mutate(pop = max(pop)) %>%
  mutate(n_T = NA, n_P1 = NA, n_P2 = NA, GOPC = NA) %>%
  # mutate(time = time - min(time) + 1) %>%
  ungroup()

hk_dat <- hk_dat %>%
  complete(Week, season) %>%
  filter(!(Week == 53 & Year != 2016)) %>%
  select(time:GOPC, Week:season, pop) %>%
  mutate(Year = if_else(is.na(Year) & season == 's13-14', 2013, Year),
         Year = if_else(is.na(Year) & season == 's15-16', 2015, Year)) %>%
  # time = if_else(is.na(time) & season == 's13-14', Week - 52, time),
  # time = if_else(is.na(time) & season == 's15-16', Week - 45, time)) %>%
  arrange(Year, Week) %>%
  group_by(season) %>%
  mutate(time = if_else(Week >= start_week, Week - start_week + 1, Week + length(start_week:max(Week))),
         pop = min(pop, na.rm = TRUE)) %>%
  ungroup()

hk_dat <- hk_dat %>%
  rename('i_ILI' = 'GOPC')

# Get climate data:
# dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
dat_clim <- read_csv('data/formatted/clim_dat_fr_NORM.csv')
hk_dat_USE <- hk_dat %>%
  inner_join(dat_clim,
             by = c('Year' = 'year',
                    'Week' = 'week')) %>%
  select(time:pop, temp, ah, rh)
rm(dat_clim)

# Run:

# Perform model checks:
source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')

# Set desired model parameter values:
model_params <- mle %>%
  dplyr::select(rho1:eta_ah2, contains(yr)) %>%
  rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
  unlist()
model_params[c('Ri1', 'Ri2', 'I10', 'I20')] <- temp_params %>% filter(season == yr) %>% select(Ri1:I20) %>% unlist()

model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
# model_params <- c(model_params, 1.0, unname(model_params['delta1']), vacc_eff)
# model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')

# Get data frame for current season and create pomp model:
dat_pomp <- hk_dat_USE %>% filter(season == yr)

resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                      Ri1_max = Ri_max1,
                                      Ri2_max = Ri_max2,
                                      d2_max = d2_max,
                                      t0_eff = 0,
                                      debug_bool = debug_bool)

resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)

model_params <- parmat(params = coef(resp_mod), nrep = 2)

# Set vaccination coverage:
model_params['p_vacc', ] <- c(0, p_vacc)

# Initiate results data frame:
res <- NULL

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
write_rds(res, paste0('results/vaccination_simulation_study/simulations/sim_determ_', yr, '_', p_vacc * 100, 'perc_TEMPERATE.rds'))
print('Done!')

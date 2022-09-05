# ---------------------------------------------------------------------------------------------------------------------
# Select parameter sets to yield outbreaks that resemble those in temperate regions
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(pomp)
library(tgp)
library(testthat)
library(patchwork)

# How many parameter sets to try?
n_lhs <- 1000

# Set parameters for run:
vir1 <- 'flu_h1'
vir2 <- 'rsv'
seasons <- c('s13-14', 's15-16', 's16-17', 's17-18', 's18-19')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

debug_bool <- FALSE
lag_val <- 0

# Load relevant functions:
source('src/functions/functions_flu_RSV.R')
source('src/vaccination_simulation_study/functions_simulation_study.R')

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format all necessary data

# Get MLEs for each season in Hong Kong:
mle <- read_rds('results/MLEs_flu_h1.rds')[1, ]

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
  rename('i_ILI' = 'GOPC')# %>%
# mutate(i_ILI = i_ILI / 1000)

# Get climate data:
# dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
dat_clim <- read_csv('data/formatted/clim_dat_fr_NORM.csv')
hk_dat <- hk_dat %>%
  inner_join(dat_clim,
             by = c('Year' = 'year',
                    'Week' = 'week')) %>%
  select(time:pop, temp, ah, rh)
rm(dat_clim)

# ---------------------------------------------------------------------------------------------------------------------

# For each season, try selecting values of I10/I20 that yield temperate patterns
# Want: pt_flu 53-60 (14-21), pt_rsv 46-55 (7-16), pt_diff -2-10

# Loop through seasons:
res_list = p_list = vector('list', length = length(seasons))
for (yr_index in 1:length(seasons)) {
  
  # Get season:
  yr <- seasons[yr_index]
  print(yr)
  
  # Get season-specific data:
  dat_temp <- hk_dat %>% filter(season == yr)
  
  # Get season-specific MLE:
  model_params <- mle %>%
    dplyr::select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  # Load model and update parameter values:
  resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_temp,
                                        Ri1_max = Ri_max1,
                                        Ri2_max = Ri_max2,
                                        d2_max = d2_max,
                                        t0_eff = 0,
                                        debug_bool = debug_bool)
  resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = FALSE)
  
  # Run model at MLE:
  sim_temp <- trajectory(resp_mod, format = 'data.frame') %>%
    select(time, H1:H2)
  
  p1 <- ggplot(data = sim_temp %>%
                 pivot_longer(H1:H2,
                              names_to = 'virus',
                              values_to = 'val'),
               aes(x = time, y = val, color = virus)) +
    geom_line() + theme_classic()
  # print(p1)
  
  # Set seed:
  set.seed(10943765)
  
  # Generate set of model parameters to try (Ri and I0 only):
  param_bound <- cbind(c(1.1, 1.1, log(1e-7), log(1e-7)),
                       c(2.0, 2.0, log(1e-3), log(5e-3)))
  params_test <- lhs(n_lhs, param_bound) %>%
    as.data.frame()
  names(params_test) <- c('Ri1', 'Ri2', 'I10', 'I20')
  params_test$I10 <- exp(params_test$I10)
  params_test$I20 <- exp(params_test$I20)
  
  model_params <- parmat(params = coef(resp_mod), nrep = n_lhs)
  model_params[c('Ri1', 'Ri2', 'I10', 'I20'), ] <- t(params_test)
  
  # # Generate set of model parameters to try (+ R10 and R20):
  # param_bound <- cbind(c(1.1, 1.1, log(1e-7), log(1e-7), 0, 0),
  #                      c(2.0, 2.0, log(1e-4), log(5e-3), 0.6, 0.6))
  # params_test <- lhs(n_lhs, param_bound) %>%
  #   as.data.frame()
  # names(params_test) <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20')
  # params_test$I10 <- exp(params_test$I10)
  # params_test$I20 <- exp(params_test$I20)
  # 
  # params_test$R120 <- coef(resp_mod, 'R120')
  # 
  # check_sum <- params_test %>%
  #   mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
  #   pull(sum) %>%
  #   max()
  # 
  # while (check_sum > 1.0) {
  # 
  #   params_test <- params_test %>%
  #     mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
  #     mutate(R10 = if_else(sum >= 1.0, R10 - ((sum - 0.9995) * (R10 / sum)), R10),
  #            R20 = if_else(sum >= 1.0, R20 - ((sum - 0.9995) * (R20 / sum)), R20),
  #            R120 = if_else(sum >= 1.0, R120 - ((sum - 0.9995) * (R120 / sum)), R120)) %>%
  #   select(-sum)
  # 
  #   check_sum <- params_test %>%
  #     mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
  #     pull(sum) %>%
  #     max()
  # 
  # }
  # 
  # model_params <- parmat(params = coef(resp_mod), nrep = n_lhs)
  # model_params[c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120'), ] <- t(params_test)
  
  # Run simulations:
  sim_temp <- trajectory(resp_mod, params = model_params, format = 'data.frame') %>%
    select(time:.id, H1:H2)
  
  # Calculate outbreak metrics of interest:
  metrics_temp <- sim_temp %>%
    group_by(.id) %>%
    summarise(pt1 = which.max(H1),
              pt2 = which.max(H2),
              ar1 = sum(H1),
              ar2 = sum(H2)) %>%
    mutate(pt_diff = pt1 - pt2)
  
  # Select only "realistic" temperate outbreaks:
  sims_select <- metrics_temp %>%
    filter(ar1 > 1e-4 & ar2 > 1e-4) %>%
    filter(pt_diff >= -2 & pt_diff <= 10) %>%
    # filter(pt1 >= 14 & pt1 <= 21) %>%
    filter(pt2 >= 7 & pt2 <= 16) %>%
    pull(.id)
  print(length(sims_select))
  
  params_select <- params_test[sims_select, ] %>%
    mutate(.id = sims_select) %>%
    inner_join(metrics_temp, by = '.id')
  
  # Plot "realistic" outbreaks:
  p2 <- ggplot(data = sim_temp %>%
                 filter(.id %in% sims_select) %>%
                 pivot_longer(H1:H2, names_to = 'virus', values_to = 'val'),
               aes(x = time, y = val, color = virus)) +
    geom_line() + theme_classic() +
    facet_wrap(~ .id, scales = 'free_y')
  p_list[[yr_index]] <- p2
  
  # Store possible parameter values in list:
  res_list[[yr_index]] <- params_select
  
}

print(p_list)

# ---------------------------------------------------------------------------------------------------------------------

# Choose "temperate" parameter sets

# Select parameter sets to use (FR climate data):
res_to_use <- bind_rows(res_list[[1]] %>% filter(.id == 970) %>% mutate(season = seasons[1]),
                        res_list[[2]] %>% filter(.id == 509) %>% mutate(season = seasons[2]),
                        res_list[[3]] %>% filter(.id == 988) %>% mutate(season = seasons[3]),
                        res_list[[4]] %>% filter(.id == 520) %>% mutate(season = seasons[4]), # or 520 - smaller RSV outbreak, but only difference of 5 vs. 7 weeks
                        res_list[[5]] %>% filter(.id == 846) %>% mutate(season = seasons[5])) %>%
  arrange(season)

write_csv(res_to_use, file = 'results/vaccine_simulation_study/temperate_params_to_use.csv')

# # Select parameter sets to use (HK climate data):
# res_to_use <- res_list[[3]] %>%
#   mutate(season = seasons[3])
# 
# res_to_use <- res_to_use %>%
#   bind_rows(res_list[[1]] %>% filter(.id == 782) %>% mutate(season = seasons[1]),
#             res_list[[2]] %>% filter(.id == 267) %>% mutate(season = seasons[2]),
#             res_list[[4]] %>% filter(.id == 868) %>% mutate(season = seasons[4]),
#             res_list[[5]] %>% filter(.id == 482) %>% mutate(season = seasons[5])) %>%
#   arrange(season)
# 
# write_csv(res_to_use, file = 'results/vaccine_simulation_study/temperate_params_to_use.csv')

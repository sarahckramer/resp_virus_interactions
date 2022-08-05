# ---------------------------------------------------------------------------------------------------------------------
# Set up model of flu/RSV transmission
# Note: Adapted from code of Dr. Matthieu Domenech de Celles
# ---------------------------------------------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)

# Load required functions:
source('src/functions/functions_flu_RSV.R')
source('src/functions/test_code.R')

# Set lag_val if it doesn't exist:
if (!exists('lag_val')) {
  lag_val <- 0
}

# Fit to observed data, not synthetic data from age-structured model:
if (!exists('age_structured')) {
  age_structured <- FALSE
}

# ---------------------------------------------------------------------------------------------------------------------

# Load and format data

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')

# Get data of interest:
dat_pomp <- hk_dat[[paste(str_sub(vir1, 5, str_length(vir1)), vir2, sep = '_')]] %>%
  filter(season == yr)
nrow_check <- nrow(dat_pomp)

if (age_structured) {
  
  dat_pomp <- read_csv('data/age_structured_SA/synthetic_obs_combined.csv') %>%
    filter(season == yr) %>%
    select(time:n_T, obs1_s1b, obs2_s1b, i_ILI:pop) %>%
    rename('n_P1' = 'obs1_s1b',
           'n_P2' = 'obs2_s1b')
  nrow_check <- nrow(dat_pomp)
  print(dat_pomp)
  
}

# Get climate data:
dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
if (lag_val > 0) {
  dat_clim <- dat_clim %>%
    mutate(temp = lag(temp, n = lag_val),
           ah = lag(ah, n = lag_val))
}
dat_pomp <- dat_pomp %>%
  inner_join(dat_clim,
             by = c('Year' = 'year',
                    'Week' = 'week')) %>%
  select(time:pop, temp, ah, rh)
expect_true(nrow(dat_pomp) == nrow_check)
rm(dat_clim)

# If no data for this season, skip:
if (nrow(dat_pomp) > 0) {
  
  # Format data:
  if (!age_structured) {
    dat_pomp <- dat_pomp %>%
      rename('i_ILI' = 'GOPC') %>%
      mutate(i_ILI = i_ILI / 1000)
    # https://www.censtatd.gov.hk/en/
    # https://www.censtatd.gov.hk/en/web_table.html?id=1A#
  }
  
  # Plot data:
  if (debug_bool) {
    # Plot ILI incidence:
    p1 <- ggplot(data = dat_pomp, aes(x = time, y = i_ILI)) + geom_line() +
      labs(x = 'Time (Weeks)', y = 'ILI Incidence Rate (per 1000 Consultations)') +
      theme_classic()
    print(p1)
    
    p2 <- ggplot(data = dat_pomp %>% pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'n_pos'),
                 aes(x = time, y = n_pos / n_T, color = virus)) +
      geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
      theme_classic()
    print(p2)
  }
  
  # ---------------------------------------------------------------------------------------------------------------------
  
  # Create pomp model and run basic model checks
  
  # Create model:
  resp_mod <- create_SITRxSITR_mod(dat = dat_pomp,
                                   Ri1_max = Ri_max1,
                                   Ri2_max = Ri_max2,
                                   d2_max = d2_max,
                                   debug_bool = debug_bool)
  
  # Check transformations:
  check_transformations(resp_mod)
  
  # Check parameters:
  check_params(resp_mod)
  
  # Check initial conditions:
  expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N'))))
  
  # Check constant population size:
  check_correct_N_CONST(resp_mod, unique(dat_pomp$pop))
  
  # Run deterministic simulation:
  sim_determ <- trajectory(object = resp_mod, format = 'data.frame') %>%
    dplyr::select(H1:.id) %>%
    pivot_longer(H1:H2, names_to = 'Vir', values_to = 'Inc')
  p3 <- ggplot(data = sim_determ, aes(x = time, y = Inc, group = Vir, col = Vir)) +
    geom_line() + geom_point() +
    labs(x = 'Time (Weeks)', y = 'Incidence', col = 'Virus') +
    theme_classic()
  if (debug_bool) print(p3)
  
  # Run stochastic simulation and check that obs never more than n_samp:
  p4 <- check_obs_lessthan_samples(resp_mod)
  if (debug_bool) print(p4)
  
  # Check that measurement density model works:
  ll <- logLik(traj_objfun(data = resp_mod))
  if (debug_bool) print(ll)
  
  # Check that dynamics are independent when there is no interaction:
  p5 <- check_independent_dynamics(resp_mod)
  if (debug_bool) print(p5)
  
  # Quick exploration of how interactions impact dynamics:
  p6 <- quick_explore_interaction(resp_mod, c(0, 0.2, 0.4, 0.6, 0.8, 1.0), n_sim = 5)
  if (debug_bool) do.call('grid.arrange', c(p6, ncol = 2))
  
  # Clean up:
  rm(sim_determ, p3, p4, p5, p6, ll)
  
}

# ---------------------------------------------------------------------------------------------------------------------

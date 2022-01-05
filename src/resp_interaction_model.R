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

# Set early start value if it doesn't exist:
if (!exists('early_start_val')) {
  early_start_val <- FALSE
}

# ---------------------------------------------------------------------------------------------------------------------

# Load and format data

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak_ALT_leadNAs.rds')
# fr_dat <- read_rds('data/formatted/GROG_pop_vir_ari_dat_2003-4_2013-14.rds')

# Get data of interest:
dat_pomp <- hk_dat[[paste(str_sub(vir1, 5, str_length(vir1)), vir2, sep = '_')]] %>%
  filter(season == yr)
nrow_check <- nrow(dat_pomp)

# Get climate data:
dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
dat_pomp <- dat_pomp %>%
  inner_join(dat_clim,
             by = c('Year' = 'year',
                    'Week' = 'week')) %>%
  select(time:season, temp, ah, rh)
expect_true(nrow(dat_pomp) == nrow_check)
rm(dat_clim)

# # Format data:
# formatted_dat <- prepare_data('flu_A', vir2, 2006, fr_dat, early_start = early_start_val)
# dat_full <- formatted_dat[[1]]
# dat_pomp <- formatted_dat[[2]]
# rm(fr_dat, formatted_dat)

# If no data for this season, skip:
if (nrow(dat_pomp) > 0) {
  
  # Format data:
  dat_pomp <- dat_pomp %>%
    rename('i_ARI' = 'GOPC') %>%
    mutate(i_ARI = i_ARI / 1000,
           pop = 7071600)
  # https://www.censtatd.gov.hk/en/
  # https://www.censtatd.gov.hk/en/web_table.html?id=1A#
  
  # Plot data:
  if (debug_bool) {
    # Plot ILI incidence:
    p1 <- ggplot(data = dat_pomp, aes(x = time, y = i_ARI)) + geom_line() +
      labs(x = 'Time (Weeks)', y = 'ILI Incidence Rate (per 1000 Consultations)') +
      theme_classic()
    print(p1)
    
    p2 <- ggplot(data = dat_pomp %>% pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'n_pos'),
                 aes(x = time, y = n_pos / n_T, color = virus)) +
      geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
      theme_classic()
    print(p2)
  }
  
  # if(debug_bool) {
  #   # Plot ARI incidence rate (based on total or effective pop size) 
  #   p1 <- ggplot(data = dat_full %>% pivot_longer(cols = c("i_ARI", "i_ARI_wrong"), names_to = "var", values_to = "val"), 
  #                mapping = aes(x = week_date, y = 100 * val, color = var)) + 
  #     geom_line() + 
  #     labs(x = "Time (weeks)", y = "ARI incidence rate (per week per 100)") +
  #     theme_classic()
  #   print(p1)
  #   
  #   # Plot positivity fraction
  #   p2 <- ggplot(data = dat_full, 
  #                mapping = aes(x = week_date, y = n_pos / n_samp, color = virus)) + 
  #     geom_line() + 
  #     labs(x = "Time (weeks)", y = "Positivity fraction") +
  #     theme_classic()
  #   print(p2)
  # }
  
  # ---------------------------------------------------------------------------------------------------------------------
  
  # Create pomp model and run basic model checks
  
  # Create model:
  resp_mod <- create_SITRxSITR_mod(dat = dat_pomp,
                                   Ri1_max = Ri_max1,
                                   Ri2_max = Ri_max2,
                                   delta_min = delta_min,
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

# ---------------------------------------------------------------------------------------------------------------------
# Set up model of flu/RSV transmission
# Note: Adapted from code of Dr. Matthieu Domenech de Celles
# ---------------------------------------------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(pomp)
library(viridis)

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
fr_dat <- read_rds('data/formatted/GROG_pop_vir_ari_dat_2003-4_2013-14.rds')

# # Visualize:
# ggplot(data = fr_dat, aes(x = week_date, y = ira_rate, group = agecat, col = agecat)) + geom_line() + facet_grid(area ~ type, scales = 'free_y') + theme_classic()
# ggplot(data = fr_dat, aes(x = week_date, y = n_pos, group = agecat, col = agecat)) + geom_line() + facet_grid(area ~ type, scales = 'free_y') + theme_classic()

# Format data:
formatted_dat <- prepare_data(vir1, vir2, yr, fr_dat, early_start = early_start_val)
dat_full <- formatted_dat[[1]]
dat_pomp <- formatted_dat[[2]]
rm(fr_dat, formatted_dat)

# Plot data:
if(debug_bool) {
  # Plot ARI incidence rate (based on total or effective pop size) 
  p1 <- ggplot(data = dat_full %>% pivot_longer(cols = c("i_ARI", "i_ARI_wrong"), names_to = "var", values_to = "val"), 
               mapping = aes(x = week_date, y = 100 * val, color = var)) + 
    geom_line() + 
    labs(x = "Time (weeks)", y = "ARI incidence rate (per week per 100)") +
    theme_classic()
  print(p1)
  
  # Plot positivity fraction
  p2 <- ggplot(data = dat_full, 
               mapping = aes(x = week_date, y = n_pos / n_samp, color = virus)) + 
    geom_line() + 
    labs(x = "Time (weeks)", y = "Positivity fraction") +
    theme_classic()
  print(p2)
}

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

# # Run stochastic simulation and check that obs never more than n_samp:
# p4 <- check_obs_lessthan_samples(resp_mod)
# if (debug_bool) print(p4)

# Check that measurement density model works:
ll <- logLik(traj_objfun(data = resp_mod))
if (debug_bool) print(ll)

# Check that dynamics are independent when there is no interaction:
p5 <- check_independent_dynamics(resp_mod)
if (debug_bool) print(p5)

# Clean up:
rm(sim_determ, p3, p5, ll)

# ---------------------------------------------------------------------------------------------------------------------

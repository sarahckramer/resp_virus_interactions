# ---------------------------------------------------------------------------------------------------------------------
# Generate synthetic data for simulation study
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tgp)

# Set seed:
set.seed(1902574195)

# Set global parameters:
n_lhs <- 500
n_sim <- 5

# Load relevant functions:
source('src/functions/functions_sim_dat.R')

# ---------------------------------------------------------------------------------------------------------------------

# Prep pomp object

# Set necessary parameter values:
vir1 <- 'flu_A' # 'flu_A', 'flu_B'
vir2 <- 'rsv'
yr <- 2006

debug_bool <- FALSE

Ri_max1 <- 2.0
Ri_max2 <- 2.0
delta_min <- 7 / 60.0

# Load pomp object:
source('src/resp_interaction_model.R')

# ---------------------------------------------------------------------------------------------------------------------

# Get realistic parameter value sets with no interaction

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 0.75, 0.75, 0.75))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = TRUE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]

# Update reporting rate for RSV:
parms_matrix['rho2', ] <- 0.10

# Store parameter matrix:
parms_matrix_ORIG <- parms_matrix
rm(parm_sets)

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# # Plot results:
# plot_res(sim_res[[2]])

# Choose 1-3 specific outbreaks (higher and lower R20.tot?):
sim_dat <- sim_res[[1]] %>% filter(real)
sim_met <- sim_res[[2]] %>%
  filter(real) %>%
  mutate(R10.tot = R10 + R120,
         R20.tot = R20 + R120)

ids_r2_high <- sim_met %>%
  filter(R20.tot > 0.75 & R20.tot < 0.85) %>%
  pull(id.parm) %>%
  unique() # IDs 48, 243, 333, 375, 468
ids_r2_low <- sim_met %>%
  filter(R20.tot < 0.35) %>%
  pull(id.parm) %>%
  unique()

sim_met %>% filter(id.parm %in% ids_r2_high) %>% select(id.parm, Ri1:R120, R10.tot:R20.tot) %>% unique()
sim_met %>% filter(id.parm %in% ids_r2_low) %>% select(id.parm, Ri1:R120, R10.tot:R20.tot) %>% unique()

sim_dat_temp <- sim_dat %>% filter(id.parm %in% ids_r2_high)
ggplot(data = sim_dat_temp) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()
sim_dat_temp <- sim_dat %>% filter(id.parm %in% ids_r2_low)
ggplot(data = sim_dat_temp) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()

# 375 has low R120, but flu outbreak gets too big
# for high R20, 243 and 468 appear the most realistic
# 169 seems to have most realistic combo of Ri1/Ri2 and realistic R10/R20, but similar timing of flu and RSV
# 240 maybe looks a bit better; 183 even better?

# Reduce data to chosen parameter values:
ids_to_use <- c(183, 240, 243, 468)
parms_to_use <- parms_matrix[, ids_to_use]

# Clean up:
rm(sim_res, sim_dat, sim_met, ids_r2_high, ids_r2_low, sim_dat_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Now with "realistic" values for R10/R20

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Adjust values of R20/R120:
parms_matrix['R10', ] <- 0.50
parms_matrix['R20', ] <- 0.35

# Adjust rho2:
parms_matrix['rho2', ] <- 0.10

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# # Plot results:
# plot_res(sim_res[[2]], plot_title = 'Realistic R', incl_R2 = FALSE)

# Choose 1 specific outbreak (perhaps with similar params to those above?):
sim_dat <- sim_res[[1]] %>% filter(real)
sim_met <- sim_res[[2]] %>% filter(real)

ggplot(data = sim_dat) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()

# potential_ids <- c(37, 71, 87, 95, 128, 214, 303, 326, 354, 375, 380, 442, 460, 490)
potential_ids <- c(22, 24, 47, 67, 70, 82, 87, 100, 113, 154, 189, 202, 214, 220, 231, 256, 314, 368, 374, 416, 447, 448)
potential_ids <- c(24, 47, 100, 154, 256, 447, 448)
ggplot(data = sim_dat %>% filter(id.parm %in% potential_ids)) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()
sim_met %>% filter(id.parm %in% potential_ids) %>% select(id.parm, Ri1:I20) %>% unique()

# Reduce data to chosen parameter values:
ids_to_use <- c(24, 47)
parms_to_use <- cbind(parms_to_use, parms_matrix[, ids_to_use])

# Clean up:
rm(sim_res, sim_dat, sim_met)

# ---------------------------------------------------------------------------------------------------------------------

# Now with R20=R120=0

# Set parameter names:
param_names <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')

# Set parameter ranges:
param_bound <- cbind(c(1.0, 1.0, 0, 0, 0),
                     c(Ri_max1, Ri_max2, 1e-3, 1e-2, 0.75))

# Draw parameter sets using LHS:
parm_sets <- generate_param_lhs(param_names, param_bound, n_lhs, resp_mod, adjust_init = FALSE)
parms_df <- parm_sets[[1]]
parms_matrix <- parm_sets[[2]]
rm(parm_sets)

# Adjust rho2:
parms_matrix['rho2', ] <- 0.10

# Run simulations:
sim_res <- run_and_format_simulations(resp_mod, parms_matrix, parms_df, n_sim)

# Identify realistic simulations:
sim_res <- classify_realistic(sim_res[[1]], sim_res[[2]], n_lhs)

# Plot results:
plot_res(sim_res[[2]], plot_title = 'R20=R120=0', incl_R2 = FALSE)

# Choose 1-3 specific outbreaks (higher and lower R20.tot?):
sim_dat <- sim_res[[1]] %>% filter(real)
sim_met <- sim_res[[2]] %>% filter(real)

sim_met %>% filter(R10 > 0.35 & R10 < 0.55) %>% unique()

# 1.3, 1.6, 0.0002, 0.005, 0.5 works okay with rho2 = 0.1?

ggplot(data = sim_dat) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()

potential_ids <- c(34, 78, 89, 186, 191, 236, 365, 388)
ggplot(data = sim_dat %>% filter(id.parm %in% potential_ids)) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
  theme_classic()
sim_met %>% filter(id.parm %in% potential_ids) %>% select(id.parm, Ri1:R10) %>% unique()

# Reduce data to chosen parameter values:
ids_to_use <- c(78)
parms_to_use <- cbind(parms_to_use, parms_matrix[, ids_to_use])

# Clean up:
rm(sim_res, sim_dat, sim_met)

# ---------------------------------------------------------------------------------------------------------------------

# # Generate "pandemic" season
# 
# # Get results from single season runs:
# r1_res <- read_rds('results/120821_rho_015/traj_match_round1_byvirseas_MLE.rds') %>%
#   bind_rows()
# 
# ri1_pan <- r1_res$Ri1[5]
# i01_pan <- r1_res$I10[5]
# r10_pan <- r1_res$R10[5]
# i02_pan <- r1_res$I20[5]
# 
# # r1_res <- r1_res[1:9, ]
# r1_res <- r1_res %>%
#   filter(I10 < 0.0007) %>%
#   select(Ri2, I20, R20, R120)
# # r1_res <- r1_res[5, ]
# 
# # Add to parameter matrix:
# parms_to_use <- cbind(parms_to_use, coef(resp_mod))
# parms_to_use[c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120'), 8] <- c(ri1_pan, median(r1_res$Ri2), i01_pan, median(r1_res$I20),
#                                                                           r10_pan, median(r1_res$R20), median(r1_res$R120))
# parms_to_use['rho2', 8] <- 0.1
# 
# # Fit with rho2=0.15, not 0.10 - adjust R120 downward:
# parms_to_use['R120', 8] <- 0.4
# # parms_to_use['R20', 8] <- 0.3329981
# # parms_to_use['Ri2', 8] <- 1.618022
# parms_to_use['I20', 8] <- median(r1_res$I20) #i02_pan # 0.0001551751
# parms_to_use['theta_lambda1', 8] <- 0.0
# parms_to_use['delta', 8] <- 7 / 60
# 
# # convert to data frame:
# parms_df <- parms_to_use %>%
#   t() %>%
#   as_tibble() %>%
#   select(Ri1:Ri2, I10:R120) %>%
#   mutate(id.parm = seq(1, ncol(parms_to_use)))
# 
# # Run simulations:
# sim_res <- run_and_format_simulations(resp_mod, parms_to_use, parms_df, n_sim)
# sim_dat <- sim_res[[1]]
# 
# # Plot:
# ggplot(data = sim_dat %>% filter(id.parm == 8)) + geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
#   geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') + facet_wrap(~ id.parm) +
#   theme_classic() + scale_y_sqrt()
# 
# # Make sure to run the pandemic seasons forward using sampling scheme from 2010!
# # Other seasons can maybe each take a random year?
# 
# # Seems like, to explain that level of a delay in the RSV peak, would need an unrealistically long interaction time
# # Do we perhaps hypothesize that an interaction alone couldn't have driven the difference?
# 
# # It is likely that, since the pandemic began before week 40, a decrease in initial I for RSV could also be due to
# # an interaction like this. But would the fitting process be capable of distinguishing between a lower I20 with
# # and without an interaction?
# # This might be something to test in a separate synthetic analysis, and not in the "main" one.

# ---------------------------------------------------------------------------------------------------------------------

# Check whether "realistic" parameters are consistently somewhat realistic regardless of year's sampling pattern








# ---------------------------------------------------------------------------------------------------------------------

# Compile and save all parameter sets








# which of these will actually allow an effect of theta_lambda? - check!; how does using different years impact outcomes?

# ---------------------------------------------------------------------------------------------------------------------

# Generate outbreaks with no interaction, range of theta_lambda1, and range of theta_lambda2










# ---------------------------------------------------------------------------------------------------------------------

# Clean up



# ---------------------------------------------------------------------------------------------------------------------


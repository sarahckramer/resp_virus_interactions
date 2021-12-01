# ---------------------------------------------------------------------------------------------------------------------
# Compile results from s13-14 using leading NAs and results from other seasons
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# ---------------------------------------------------------------------------------------------------------------------

# For runs including interaction

# Read in top parameter fits:
pars_list <- read_rds('results/round1_interaction/traj_match_round1_byvirseas_TOP.rds')
pars_list_1314 <- read_rds('results/round1_interaction_leadNAs/traj_match_round1_byvirseas_TOP.rds')

# Replace old results for s13-14 with results using leading NAs:
pars_list[names(pars_list) %in% names(pars_list_1314)] <- pars_list_1314

# Write new results to file:
write_rds(pars_list, file = 'results/traj_match_round1_byvirseas_TOP.rds')

# ---------------------------------------------------------------------------------------------------------------------

# Repeat for runs without interaction

# Read in top parameter fits:
pars_list <- read_rds('results/round1_rho-logit/traj_match_round1_byvirseas_TOP.rds')
pars_list_1314 <- read_rds('results/round1_rho-logit_leadNAs/traj_match_round1_byvirseas_TOP.rds')

# Replace old results for s13-14 with results using leading NAs:
pars_list[names(pars_list) %in% names(pars_list_1314)] <- pars_list_1314

# Write new results to file:
write_rds(pars_list, file = 'results/traj_match_round1_byvirseas_TOP_noInt.rds')

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

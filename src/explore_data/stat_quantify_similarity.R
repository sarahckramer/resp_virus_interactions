# ---------------------------------------------------------------------------------------------------------------------
# Statistically quantify similarity between time series of various pathogens
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in and format data:
dat_hk <- read_csv('data/formatted/dat_hk.csv')
dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_h1, n_b:n_rsv, n_rhino_est_rnd, GOPC) %>%
  mutate(prop_h1 = n_h1 / n_samp * 100,
         prop_b = n_b / n_samp * 100,
         prop_h1_plus_b = n_h1_b / n_samp * 100,
         prop_rsv = n_rsv / n_samp * 100,
         prop_rhino = n_rhino_est_rnd / n_samp * 100) %>%
  select(-n_samp) %>%
  select(Time:n_rhino_est_rnd, prop_h1:prop_rhino, GOPC)

# Calculate correlations between pathogens:
cor.test(dat_hk$prop_h1_plus_b, dat_hk$prop_rsv, method = 'kendall')
cor.test(dat_hk$prop_h1_plus_b, dat_hk$prop_rhino, method = 'kendall')
cor.test(dat_hk$prop_rsv, dat_hk$prop_rhino, method = 'kendall')

# Calculate peak timing difference between flu and RSV:
# This is done in "calculate_outbreak_metrics.R"

# Clean up:
rm(list = ls())

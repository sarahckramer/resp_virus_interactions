# ---------------------------------------------------------------------------------------------------------------------
# Statistically quantify similarity between time series of various pathogens
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(NSM3)

# Read in and format data:
dat_hk <- read_csv('data/formatted/dat_hk.csv')
dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  dplyr::select(Time:n_samp, n_h1:n_rsv, n_rhino_est_rnd, GOPC) %>%
  mutate(prop_h1 = n_h1 / n_samp * 100,
         prop_h3 = n_h3 / n_samp * 100,
         prop_b = n_b / n_samp * 100,
         prop_h1_plus_b = n_h1_b / n_samp * 100,
         prop_rsv = n_rsv / n_samp * 100,
         prop_rhino = n_rhino_est_rnd / n_samp * 100) %>%
  dplyr::select(-n_samp) %>%
  dplyr::select(Time:n_rhino_est_rnd, prop_h1:prop_rhino, GOPC)

dat_can <- read_csv('data/formatted/dat_canada.csv')
dat_can <- dat_can %>%
  dplyr::select(time, year:week, n_T1:i_ILI) %>%
  mutate(prop_P1 = n_P1 / n_T1 * 100,
         prop_P2 = n_P2 / n_T2 * 100,
         prop_rhino = rhino / n_test_rhino * 100) %>%
  dplyr::select(-c(n_T1, n_T2, n_test_rhino)) %>%
  dplyr::select(time:n_P2, prop_P1:prop_rhino, i_ILI)

# Calculate correlations between pathogens:
cor.test(dat_hk$prop_h1_plus_b, dat_hk$prop_rsv, method = 'kendall')
kendall.ci(dat_hk$prop_h1_plus_b, dat_hk$prop_rsv)

cor.test(dat_hk$prop_h1_plus_b, dat_hk$prop_rhino, method = 'kendall')
cor.test(dat_hk$prop_rsv, dat_hk$prop_rhino, method = 'kendall')

cor.test(dat_can$prop_P1, dat_can$prop_P2, method = 'kendall')
kendall.ci(dat_can$prop_P1, dat_can$prop_P2)

cor.test(dat_can$prop_P1, dat_can$prop_rhino, method = 'kendall')
cor.test(dat_can$prop_P2, dat_can$prop_rhino, method = 'kendall')

# Calculate peak timing difference between flu and RSV:
# This is done in "calculate_outbreak_metrics.R"

# Clean up:
rm(list = ls())

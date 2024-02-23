# ---------------------------------------------------------------------------------------------------------------------
# Calculate lag-one autocorrelation to quantify data noisiness
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in and format data:
dat_hk <- read_csv('data/formatted/dat_hk.csv')
dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_h1, n_b:n_rsv, GOPC) %>%
  mutate(prop_h1 = n_h1 / n_samp * 100,
         prop_b = n_b / n_samp * 100,
         prop_h1_plus_b = n_h1_b / n_samp * 100,
         prop_rsv = n_rsv / n_samp * 100) %>%
  select(-n_samp) %>%
  select(Time:n_rsv, prop_h1:prop_rsv, GOPC)

dat_can <- read_csv('data/formatted/dat_canada.csv')
dat_can <- dat_can %>%
  select(time, year:week, n_T1:n_T2, n_P1:n_P2, i_ILI) %>%
  mutate(prop_P1 = n_P1 / n_T1 * 100,
         prop_P2 = n_P2 / n_T2 * 100) %>%
  select(-c(n_T1, n_T2)) %>%
  select(time:n_P2, prop_P1:prop_P2, i_ILI)

# Calculate lag-one autocorrelation for raw observation numbers:
cor(dat_hk$n_h1[-length(dat_hk$n_h1)], dat_hk$n_h1[-1])
cor(dat_hk$n_b[-length(dat_hk$n_b)], dat_hk$n_b[-1])
cor(dat_hk$n_h1_b[-length(dat_hk$n_h1_b)], dat_hk$n_h1_b[-1])
cor(dat_hk$n_rsv[-length(dat_hk$n_rsv)], dat_hk$n_rsv[-1])

cor(dat_can$n_P1[-length(dat_can$n_P1)], dat_can$n_P1[-1])
cor(dat_can$n_P2[-length(dat_can$n_P2)], dat_can$n_P2[-1])

# Calculate lag-one autocorrelation for positivity rates:
cor(dat_hk$prop_h1[-length(dat_hk$prop_h1)], dat_hk$prop_h1[-1])
cor(dat_hk$prop_b[-length(dat_hk$prop_b)], dat_hk$prop_b[-1])
cor(dat_hk$prop_h1_plus_b[-length(dat_hk$prop_h1_plus_b)], dat_hk$prop_h1_plus_b[-1])
cor(dat_hk$prop_rsv[-length(dat_hk$prop_rsv)], dat_hk$prop_rsv[-1])

cor(dat_can$prop_P1[-length(dat_can$prop_P1)], dat_can$prop_P1[-1])
cor(dat_can$prop_P2[-length(dat_can$prop_P2)], dat_can$prop_P2[-1])

# Calculate lag-one autocorrelation for ILI:
cor(dat_hk$GOPC[-length(dat_hk$GOPC)], dat_hk$GOPC[-1])
cor(dat_can$i_ILI[-length(dat_can$i_ILI)], dat_can$i_ILI[-1])

# Clean up:
rm(list = ls())

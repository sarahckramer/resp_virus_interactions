# ---------------------------------------------------------------------------------------------------------------------
# Generate synthetic covariate data (ILI rates and test numbers by age group) for Hong Kong
# ---------------------------------------------------------------------------------------------------------------------

# Setup

rm(list = ls())
source("src/age_structured_SA/s-base_packages.R")
library(socialmixr)

# Functions

round_preserve_sum <- function(x, digits = 0) {
  # Function to round values in a vector while maintaining the same sum
  # Source: https://www.r-bloggers.com/2016/07/round-values-while-preserve-their-rounded-sum-in-r/
  # param x: Vector of values to be rounded
  # param digits: Round to how many decimal places?
  # returns: Vector with rounded values
  
  if (!all(is.na(x))) {
    
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y / up
    
  } else {
    x
  }
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Load and format French data

# Read in age-structured data from France:
fr_dat <- read_rds('data/formatted/GROG_pop_vir_ari_dat_2003-4_2013-14.rds')

# Format:
fr_dat <- fr_dat %>%
  mutate(virus = fct_collapse(.f = type,
                              'rsv' = 'rsv',
                              'flu_B' = 'flu_b',
                              'flu_A' = c('flu_ah1', 'flu_ah3', 'flu_a_nst'))) %>%
  group_by(epiyear, week_no, week_date, virus, type, agecat) %>%
  summarise(pop_tot = sum(pop),
            pop_eff = sum(pop[!is.na(ira_n)]),
            n_samp = sum(n_samp, na.rm = TRUE),
            n_pos = sum(n_pos, na.rm = TRUE),
            n_ARI = sum(ira_n, na.rm = TRUE)) %>%
  group_by(epiyear, week_no, week_date, virus, agecat) %>%
  summarise(pop_tot = unique(pop_tot),
            pop_eff = unique(pop_eff),
            n_ARI = unique(n_ARI),
            n_samp = unique(n_samp),
            n_pos = sum(n_pos)) %>%
  group_by(epiyear) %>% mutate(week_no = if_else(week_no < 40, week_no + max(week_no), week_no)) %>%
  ungroup() %>%
  mutate(i_ARI = n_ARI / pop_eff,
         prop_pos = n_pos / n_samp) %>%
  # select(-c(pop_eff, n_ARI)) %>%
  filter(epiyear >= 2006)

# Remove pandemic:
fr_dat <- fr_dat %>%
  filter(epiyear != 2010)

# Remove flu_B:
fr_dat <- fr_dat %>%
  filter(virus %in% c('flu_A', 'rsv'))

# ---------------------------------------------------------------------------------------------------------------------

# Calculate rates of ARI and testing by age group

# Calculate average attack rate in each age group, both actual and relative to population-level attack rate:
ari_rate_per_pop <- fr_dat %>%
  filter(virus == 'flu_A') %>%
  group_by(agecat) %>%
  summarise(n_ARI = sum(n_ARI),
            pop_eff = sum(pop_eff),
            i_ARI = sum(n_ARI) / sum(pop_eff)) %>%
  ungroup() %>%
  mutate(i_ARI_rel = (n_ARI / pop_eff) / (sum(n_ARI) / sum(pop_eff))) %>%
  select(agecat, i_ARI:i_ARI_rel)

# Calculate average rate of testing by age group:
test_rate_per_pop <- fr_dat %>%
  filter(virus == 'flu_A') %>%
  group_by(agecat) %>%
  summarise(i_samp = sum(n_samp) / sum(pop_eff)) %>%
  ungroup()

# ---------------------------------------------------------------------------------------------------------------------

# Load and format Hong Kong data

# Load infection data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_rsv %>%
  select(time:n_T, GOPC:season) %>%
  rename('i_ILI' = 'GOPC') %>%
  mutate(i_ILI = i_ILI / 1000)

# Get population size by age group:
age_limits <- c(0, 1, 5, 16, 65) # Age limits for age categories
data(polymod)
CM_all <- contact_matrix(survey = get_survey(survey = "https://doi.org/10.5281/zenodo.3874808"), 
                         age.limits = age_limits,
                         survey.pop = data.frame(lower.age.limit = age_limits, population = c(45.8, 183.2, 578.8, 5153.8, 1451.5) * 1e3),
                         symmetric = T)
N <- CM_all$demography$population 
rm(age_limits, CM_all)

# ---------------------------------------------------------------------------------------------------------------------

# Proportionally allocate cases and tests in Hong Kong

# Inflate or reduce ILI attack rate based on relative attack rates in France:
hk_dat <- hk_dat %>%
  mutate(i_ILI1 = i_ILI * ari_rate_per_pop$i_ARI_rel[1],
         i_ILI2 = i_ILI * ari_rate_per_pop$i_ARI_rel[1],
         i_ILI3 = i_ILI * ari_rate_per_pop$i_ARI_rel[2],
         i_ILI4 = i_ILI * ari_rate_per_pop$i_ARI_rel[3],
         i_ILI5 = i_ILI * ari_rate_per_pop$i_ARI_rel[4])

# Alternative scenario: Assign same ARI rate to all age groups:
# For this, just use the values in column i_ARI

# Allocate total tests proportionally to test rates in France:
hk_dat <- hk_dat %>%
  mutate(n_T1 = test_rate_per_pop$i_samp[1] * N[1],
         n_T2 = test_rate_per_pop$i_samp[1] * N[2],
         n_T3 = test_rate_per_pop$i_samp[2] * N[3],
         n_T4 = test_rate_per_pop$i_samp[3] * N[4],
         n_T5 = test_rate_per_pop$i_samp[4] * N[5],
         n_T_tot = n_T1 + n_T2 + n_T3 + n_T4 + n_T5) %>%
  mutate(n_T1 = n_T * (n_T1 / n_T_tot),
         n_T2 = n_T * (n_T2 / n_T_tot),
         n_T3 = n_T * (n_T3 / n_T_tot),
         n_T4 = n_T * (n_T4 / n_T_tot),
         n_T5 = n_T * (n_T5 / n_T_tot)) %>%
  select(-n_T_tot)

hk_dat[, c('n_T1', 'n_T2', 'n_T3', 'n_T4', 'n_T5')] <- apply(hk_dat[, c('n_T1', 'n_T2', 'n_T3', 'n_T4', 'n_T5')], 1, round_preserve_sum) %>% t()

hk_dat %>%
  mutate(n_T_tot = n_T1 + n_T2 + n_T3 + n_T4 + n_T5,
         diff_0 = (n_T_tot - n_T == 0)) %>%
  pull(diff_0) %>%
  all(na.rm = TRUE) %>%
  expect_true()

# Alternative scenario: Only two youngest age groups report:
hk_dat <- hk_dat %>%
  mutate(n_T1_s2 = test_rate_per_pop$i_samp[1] * N[1],
         n_T2_s2 = test_rate_per_pop$i_samp[1] * N[2],
         n_T_tot = n_T1_s2 + n_T2_s2) %>%
  mutate(n_T1_s2 = n_T * (n_T1_s2 / n_T_tot),
         n_T2_s2 = n_T * (n_T2_s2 / n_T_tot)) %>%
  select(-n_T_tot)

hk_dat[, c('n_T1_s2', 'n_T2_s2')] <- apply(hk_dat[, c('n_T1_s2', 'n_T2_s2')], 1, round_preserve_sum) %>% t()

hk_dat %>%
  mutate(n_T_tot = n_T1_s2 + n_T2_s2,
         diff_0 = (n_T_tot - n_T == 0)) %>%
  pull(diff_0) %>%
  all(na.rm = TRUE) %>%
  expect_true()

# Reorder columns:
hk_dat <- hk_dat %>%
  select(time:Year, Week:season, i_ILI, n_T, i_ILI1:n_T2_s2)

# Save covariate "data":
if (!dir.exists('data/age_structured_SA/')) {
  dir.create('data/age_structured_SA/')
}
write_csv(hk_dat, file = 'data/age_structured_SA/synthetic_covariate_data.csv')

# Clean up:
rm(list = ls())

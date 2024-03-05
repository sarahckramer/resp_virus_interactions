# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates, peak timings, and outbreak concentration for flu and RSV outbreaks
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv
can_dat <- read_csv('data/formatted/dat_canada.csv')
us_dat <- read_rds('data/formatted/dat_us_byRegion.rds')$'Region 8'

# Format and join:
hk_dat <- hk_dat %>%
  select(time, season, pop, n_T:n_P2) %>%
  mutate(n_T1 = n_T,
         n_T2 = n_T) %>%
  select(time:pop, n_T1:n_T2, n_P1:n_P2) %>%
  mutate(loc = 'hk')

can_dat <- can_dat %>%
  select(time, season, pop, n_T1:n_T2, n_P1:n_P2) %>%
  mutate(loc = 'can')

us_dat <- us_dat %>%
  select(time, season, pop, n_T1:n_P2) %>%
  mutate(loc = 'us')

dat <- rbind(hk_dat, can_dat, us_dat)
rm(hk_dat, can_dat, us_dat)

# Calculate peak timing:
metrics_pt <- dat %>%
  mutate(prop1 = n_P1 / n_T1,
         prop2 = n_P2 / n_T2) %>%
  pivot_longer(prop1:prop2, names_to = 'vir', values_to = 'prop_pos') %>%
  group_by(vir, loc, season) %>%
  summarise(pt = which.max(prop_pos)) %>%
  mutate(pt = if_else(loc == 'hk', pt + 45, pt),
         pt = if_else(loc == 'can', pt + 34, pt),
         pt = if_else(loc == 'us', pt + 39, pt)) %>%
  mutate(vir = if_else(vir == 'prop1', 'influenza', 'rsv'))

# Calculate attack rates:
metrics_ar <- dat %>%
  pivot_longer(n_P1:n_P2, names_to = 'vir', values_to = 'number') %>%
  pivot_longer(n_T1:n_T2, names_to = 'vir1', values_to = 'n_T') %>%
  mutate(vir = str_sub(vir, -1),
         vir1 = str_sub(vir1, -1)) %>%
  filter(vir == vir1) %>%
  mutate(prop = number / pop) %>%
  group_by(vir, loc, season) %>%
  summarise(ar = sum(prop, na.rm = TRUE),
            ar_prop = sum(number, na.rm = TRUE) / sum(n_T, na.rm = TRUE)) %>%
  mutate(vir = if_else(vir == '1', 'influenza', 'rsv'))

# Calculate peak timing differences:
metrics_pt_diff <- metrics_pt %>%
  pivot_wider(names_from = vir, values_from = pt) %>%
  mutate(pt_diff = influenza - rsv) %>%
  select(loc:season, pt_diff) %>%
  mutate(vir = 'influenza')

# Calculate the number of weeks containing 75% of reported cases:
metrics_conc <- NULL
for (lc in unique(dat$loc)) {
  
  dat_loc <- dat %>%
    filter(loc == lc)
  
  for (yr in unique(dat_loc$season)) {
    
    dat_temp <- dat_loc %>%
      filter(season == yr)
    
    for (vir in c('n_P1', 'n_P2')) {
      
      case_counts_temp <- dat_temp %>% pull(vir)
      case_counts_temp[is.na(case_counts_temp)] <- 0
      
      target_sum <- sum(case_counts_temp)
      sum_cases <- 0
      which_weeks <- c()
      
      while (sum_cases < 0.75 * target_sum) {
        
        which_max <- which.max(case_counts_temp)
        max_val <- case_counts_temp[which_max]
        
        sum_cases <- sum_cases + max_val
        which_weeks <- c(which_weeks, which_max)
        
        case_counts_temp[which_max] <- 0
        
      }
      
      metrics_conc <- rbind(metrics_conc,
                            c(lc,
                              yr,
                              vir,
                              length(which_weeks),
                              min(which_weeks) + length(which_weeks) - 1 == max(which_weeks)))
      
      
    }
    
  }
  
}
rm(dat_loc, dat_temp, case_counts_temp, lc, max_val, sum_cases, target_sum, vir, which_max, which_weeks, yr)

metrics_conc <- metrics_conc %>%
  as_tibble() %>%
  rename('loc' = 'V1',
         'season' = 'V2',
         'vir' = 'V3',
         'duration' = 'V4',
         'consecutive' = 'V5') %>%
  mutate(vir = if_else(vir == 'n_P1', 'influenza', 'rsv'))

# Join all metrics:
metrics <- metrics_pt %>%
  left_join(metrics_ar, by = c('loc', 'season', 'vir')) %>%
  left_join(metrics_pt_diff, by = c('loc', 'season', 'vir')) %>%
  left_join(metrics_conc, by = c('loc', 'season', 'vir')) %>%
  select(loc:season, vir, ar:ar_prop, pt, pt_diff, duration:consecutive) %>%
  mutate(duration = as.numeric(duration),
         consecutive = as.numeric(as.logical(consecutive))) %>%
  mutate(ar = ar * 100,
         ar_prop = ar_prop * 100)

expect_true(nrow(metrics) == nrow(metrics_pt))
expect_true(nrow(metrics) == nrow(metrics_ar))
expect_true(nrow(metrics) == nrow(metrics_conc))
expect_true(nrow(metrics) == nrow(metrics_pt_diff) * 2)

rm(metrics_pt, metrics_ar, metrics_pt_diff, metrics_conc)

# Print values:
metrics <- metrics %>%
  pivot_longer(ar:consecutive, names_to = 'metric', values_to = 'val')

metrics %>%
  filter(vir == 'influenza') %>%
  group_by(metric, loc) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()

metrics %>%
  filter(vir == 'rsv', metric == 'pt') %>%
  mutate(val = val - 52) %>%
  mutate(val = if_else(season == 's16-17' & loc == 'hk', val - 1, val),
         val = if_else(season == 's14-15' & loc == 'us', val - 1, val)) %>%
  group_by(loc) %>%
  summarise(min = min(val),
            median = median(val),
            mean = mean(val),
            max = max(val)) %>%
  print()

metrics %>%
  filter(vir == 'rsv') %>%
  group_by(metric, loc) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()

metrics %>%
  filter(vir == 'rsv', metric == 'pt') %>%
  mutate(val = val - 52) %>%
  mutate(val = if_else(season == 's16-17' & loc == 'hk', val - 1, val),
         val = if_else(season == 's14-15' & loc == 'us', val - 1, val)) %>%
  group_by(loc) %>%
  summarise(min = min(val),
            median = median(val),
            mean = mean(val),
            max = max(val)) %>%
  print()

# Plot "realistic" values:
p1 <- ggplot(data = metrics %>% filter(metric != 'consecutive'),
             aes(x = loc, y = val)) +
  geom_violin(fill = 'gray90') +
  facet_grid(metric ~ vir, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value')
print(p1)

# Clean up:
rm(list = ls())

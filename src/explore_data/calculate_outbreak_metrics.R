# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates, peak timings, and outbreak concentration for flu and RSV outbreaks
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
can_dat <- read_csv('data/formatted/dat_canada.csv')

# Compile relevant data:
hk_dat <- hk_dat$h1_rsv %>%
  rename(n_h1 = n_P1,
         n_rsv = n_P2) %>%
  select(time, Year, Week, season, n_h1, n_rsv, n_T, GOPC, pop) %>%
  full_join(hk_dat$b_rsv %>%
              rename(n_b = n_P1,
                     n_rsv = n_P2) %>%
              select(time, season, n_b, n_rsv, n_T, GOPC, pop),
            by = c('time', 'season')) %>%
  full_join(hk_dat$h1_plus_b_rsv %>%
              rename(n_h1_plus_b = n_P1,
                     n_rsv = n_P2) %>%
              select(time, season, n_h1_plus_b, n_rsv, n_T, GOPC, pop),
            by = c('time', 'season'))

expect_true(all.equal(hk_dat$n_rsv.x, hk_dat$n_rsv.y))
expect_true(all.equal(hk_dat$n_T.x, hk_dat$n_T.y))
expect_true(all.equal(hk_dat$GOPC.x, hk_dat$GOPC.y))
expect_true(all.equal(hk_dat$pop.x, hk_dat$pop.y))

expect_true(all.equal(hk_dat$n_rsv, hk_dat$n_rsv.y))
expect_true(all.equal(hk_dat$n_T, hk_dat$n_T.y))
expect_true(all.equal(hk_dat$GOPC, hk_dat$GOPC.y))
expect_true(all.equal(hk_dat$pop, hk_dat$pop.y))

hk_dat <- hk_dat %>%
  select(time:n_h1, n_b, n_h1_plus_b:pop) %>%
  arrange(season)

# Calculate peak timing and attack rate metrics:
hk_metrics_pt_ar <- hk_dat %>%
  pivot_longer(n_h1:n_rsv, names_to = 'vir', values_to = 'number') %>%
  mutate(prop = number / pop,
         prop_pos = number / n_T) %>%
  select(-GOPC) %>%
  group_by(vir, season) %>%
  summarise(pt = which.max(prop_pos) + 45,
            ar = sum(prop, na.rm = TRUE),
            ar_prop = sum(number, na.rm = TRUE) / sum(n_T, na.rm = TRUE))

can_metrics_pt_ar <- can_dat %>%
  select(-c(n_test_rhino, rhino)) %>%
  pivot_longer(n_P1:n_P2, names_to = 'vir', values_to = 'number') %>%
  pivot_longer(n_T1:n_T2, names_to = 'vir1', values_to = 'n_T') %>%
  mutate(vir = str_sub(vir, 4),
         vir1 = str_sub(vir1, 4)) %>%
  filter(vir == vir1) %>%
  select(-vir1, i_ILI) %>%
  mutate(vir = if_else(vir == 1, 'influenza', 'rsv')) %>%
  mutate(prop = number / pop,
         prop_pos = number / n_T) %>%
  group_by(vir, season) %>%
  summarise(pt = which.max(prop_pos) + 34,
            ar = sum(prop, na.rm = TRUE),
            ar_prop = sum(number, na.rm = TRUE) / sum(n_T, na.rm = TRUE))

# Calculate peak timing differences:
hk_metrics_pt_diff <- hk_metrics_pt_ar %>%
  select(-c(ar, ar_prop)) %>%
  pivot_wider(names_from = vir, values_from = pt) %>%
  mutate(n_h1 = n_h1 - n_rsv,
         n_b = n_b - n_rsv,
         n_h1_plus_b = n_h1_plus_b - n_rsv) %>%
  select(-n_rsv) %>%
  pivot_longer(n_b:n_h1_plus_b, names_to = 'vir', values_to = 'pt_diff')

can_metrics_pt_diff <- can_metrics_pt_ar %>%
  select(-c(ar, ar_prop)) %>%
  pivot_wider(names_from = vir, values_from = pt) %>%
  mutate(pt_diff = influenza - rsv,
         vir = 'influenza') %>%
  select(season, vir, pt_diff)

# Calculate the number of weeks containing 75% of reported cases:
hk_metrics_conc <- NULL
for (yr in unique(hk_dat$season)) {
  
  hk_dat_temp <- hk_dat %>%
    filter(season == yr)
  
  for (vir in c('n_h1', 'n_b', 'n_h1_plus_b', 'n_rsv')) {
    
    case_counts_temp <- hk_dat_temp %>% pull(vir)
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
    
    hk_metrics_conc <- rbind(hk_metrics_conc,
                             c(vir,
                               yr,
                               length(which_weeks),
                               min(which_weeks) + length(which_weeks) - 1 == max(which_weeks)))
    
  }
  
  rm(case_counts_temp, target_sum, sum_cases, which_weeks, which_max, max_val)
  
}
rm(hk_dat_temp, yr, vir)

hk_metrics_conc <- hk_metrics_conc %>%
  as_tibble() %>%
  rename('vir' = 'V1',
         'season' = 'V2',
         'duration' = 'V3',
         'consecutive' = 'V4')

can_metrics_conc <- NULL
for (yr in unique(can_dat$season)) {
  
  can_dat_temp <- can_dat %>%
    filter(season == yr)
  
  for (vir in c('influenza', 'rsv')) {
    
    if (vir == 'influenza') {
      case_counts_temp <- can_dat_temp %>% pull(n_P1)
    } else if (vir == 'rsv') {
      case_counts_temp <- can_dat_temp %>% pull(n_P2)
    }
    
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
    
    can_metrics_conc <- rbind(can_metrics_conc,
                              c(vir,
                                yr,
                                length(which_weeks),
                                min(which_weeks) + length(which_weeks) - 1 == max(which_weeks)))
    
  }
  
  rm(case_counts_temp, target_sum, sum_cases, which_weeks, which_max, max_val)
  
}
rm(can_dat_temp, yr, vir)

can_metrics_conc <- can_metrics_conc %>%
  as_tibble() %>%
  rename('vir' = 'V1',
         'season' = 'V2',
         'duration' = 'V3',
         'consecutive' = 'V4')

# Join all metrics:
metrics_hk <- hk_metrics_pt_ar %>%
  left_join(hk_metrics_pt_diff, by = c('vir', 'season')) %>%
  left_join(hk_metrics_conc, by = c('vir', 'season')) %>%
  select(vir:season, ar:ar_prop, pt, pt_diff, duration:consecutive) %>%
  mutate(duration = as.numeric(duration),
         consecutive = as.numeric(as.logical(consecutive))) %>%
  mutate(ar = ar * 100,
         ar_prop = ar_prop * 100)
rm(hk_metrics_pt_ar, hk_metrics_pt_diff, hk_metrics_conc)

metrics_can <- can_metrics_pt_ar %>%
  left_join(can_metrics_pt_diff, by = c('vir', 'season')) %>%
  left_join(can_metrics_conc, by = c('vir', 'season')) %>%
  select(vir:season, ar:ar_prop, pt, pt_diff, duration:consecutive) %>%
  mutate(duration = as.numeric(duration),
         consecutive = as.numeric(as.logical(consecutive))) %>%
  mutate(ar = ar * 100,
         ar_prop = ar_prop * 100)
rm(can_metrics_pt_ar, can_metrics_pt_diff, can_metrics_conc)

# Combine
metrics <- metrics_hk %>%
  filter(vir %in% c('n_h1_plus_b', 'n_rsv')) %>%
  mutate(vir = if_else(vir == 'n_h1_plus_b', 'influenza', 'rsv'),
         location = 'hk') %>%
  bind_rows(metrics_can %>%
              mutate(location = 'can'))

# Print values:
metrics <- metrics %>%
  pivot_longer(ar:consecutive, names_to = 'metric', values_to = 'val')
metrics_hk_full <- metrics_hk %>%
  pivot_longer(ar:consecutive, names_to = 'metric', values_to = 'val')
rm(metrics_hk, metrics_can)

metrics %>%
  filter(vir == 'influenza') %>%
  group_by(location, metric) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()

metrics %>%
  filter(vir == 'influenza', metric == 'pt') %>%
  mutate(val = if_else(season == 's16-17', val - 53, val - 52)) %>%
  group_by(location) %>%
  summarise(min = min(val),
            median = median(val),
            mean = mean(val),
            max = max(val)) %>%
  print()

metrics %>%
  filter(vir == 'rsv') %>%
  group_by(location, metric) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()

metrics %>%
  filter(vir == 'rsv', metric == 'pt') %>%
  mutate(val = if_else(season == 's16-17', val - 53, val - 52)) %>%
  group_by(location) %>%
  summarise(min = min(val),
            median = median(val),
            mean = mean(val),
            max = max(val)) %>%
  print()

# Plot "realistic" values:
p1 <- ggplot(data = metrics %>% filter(metric != 'consecutive'),
       aes(x = vir, y = val)) +
  geom_violin(fill = 'gray90') +
  facet_grid(metric ~ location, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value')
print(p1)

p2 <- ggplot(data = metrics_hk_full %>% filter(metric != 'consecutive'),
             aes(x = vir, y = val)) +
  geom_violin(fill = 'gray90') +
  facet_wrap(~ metric, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value')
print(p2)

# Clean up:
rm(list = ls())

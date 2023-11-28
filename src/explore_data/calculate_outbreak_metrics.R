# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates, peak timings, and outbreak concentration for flu and RSV outbreaks
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')

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
metrics_pt_ar <- hk_dat %>%
  pivot_longer(n_h1:n_rsv, names_to = 'vir', values_to = 'number') %>%
  mutate(prop = number / pop,
         prop_pos = number / n_T) %>%
  select(-GOPC) %>%
  group_by(vir, season) %>%
  summarise(pt = which.max(prop_pos) + 45,
            ar = sum(prop, na.rm = TRUE),
            ar_prop = sum(number, na.rm = TRUE) / sum(n_T, na.rm = TRUE))

# Calculate peak timing differences:
metrics_pt_diff <- metrics_pt_ar %>%
  select(-c(ar, ar_prop)) %>%
  pivot_wider(names_from = vir, values_from = pt) %>%
  mutate(n_h1 = n_h1 - n_rsv,
         n_b = n_b - n_rsv,
         n_h1_plus_b = n_h1_plus_b - n_rsv) %>%
  select(-n_rsv) %>%
  pivot_longer(n_b:n_h1_plus_b, names_to = 'vir', values_to = 'pt_diff')

# Calculate the number of weeks containing 75% of reported cases:
metrics_conc <- NULL
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
    
    metrics_conc <- rbind(metrics_conc,
                          c(vir,
                            yr,
                            length(which_weeks),
                            min(which_weeks) + length(which_weeks) - 1 == max(which_weeks)))
    
  }
  
  rm(case_counts_temp, target_sum, sum_cases, which_weeks, which_max, max_val)
  
}
rm(hk_dat_temp, yr, vir)

metrics_conc <- metrics_conc %>%
  as_tibble() %>%
  rename('vir' = 'V1',
         'season' = 'V2',
         'duration' = 'V3',
         'consecutive' = 'V4')

# Join all metrics:
metrics <- metrics_pt_ar %>%
  left_join(metrics_pt_diff, by = c('vir', 'season')) %>%
  left_join(metrics_conc, by = c('vir', 'season')) %>%
  select(vir:season, ar:ar_prop, pt, pt_diff, duration, consecutive) %>%
  mutate(duration = as.numeric(duration),
         consecutive = as.numeric(as.logical(consecutive))) %>%
  mutate(ar = ar * 100,
         ar_prop = ar_prop * 100)
rm(metrics_pt_ar, metrics_pt_diff, metrics_conc)

# Print values:
metrics <- metrics %>%
  pivot_longer(ar:consecutive, names_to = 'metric', values_to = 'val')

metrics %>%
  filter(vir == 'n_h1_plus_b') %>%
  group_by(metric) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()
metrics %>%
  filter(vir == 'n_h1_plus_b', metric == 'pt') %>%
  mutate(val = if_else(season == 's16-17', val - 53, val - 52)) %>%
  pull(val) %>%
  summary()

metrics %>%
  filter(vir == 'n_rsv') %>%
  group_by(metric) %>%
  summarise(mean = mean(val),
            median = median(val),
            min = min(val),
            max = max(val),
            sum = sum(val)) %>%
  print()
metrics %>%
  filter(vir == 'n_rsv', metric == 'pt') %>%
  mutate(val = if_else(season == 's16-17', val - 53, val - 52)) %>%
  pull(val) %>%
  summary()

# Plot "realistic" values:
p1 <- ggplot(data = metrics %>% filter(metric != 'consecutive'),
             aes(x = vir, y = val, fill = metric)) +
  geom_violin() +
  facet_wrap(~ metric, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value', fill = 'Metric') +
  scale_fill_brewer(palette = 'Set1')
print(p1)

# Clean up:
rm(list = ls())

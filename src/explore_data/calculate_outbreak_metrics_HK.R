# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates and peak timings for flu and RSV outbreaks (France)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak_ALT_leadNAs.rds')

# Compile relevant data:
hk_dat <- hk_dat$h1_rsv %>%
  rename(n_h1 = n_P1,
         n_rsv = n_P2) %>%
  select(time, season, n_h1, n_rsv, n_T, GOPC) %>%
  full_join(hk_dat$b_rsv %>%
              rename(n_b = n_P1,
                     n_rsv = n_P2) %>%
              select(time, season, n_b, n_rsv, n_T, GOPC),
            by = c('time', 'season')) %>%
  mutate(n_rsv = if_else(is.na(n_rsv.x), n_rsv.y, n_rsv.x),
         n_T = if_else(is.na(n_T.x), n_T.y, n_T.x),
         GOPC = if_else(is.na(GOPC.x), GOPC.y, GOPC.x)) %>%
  select(time:n_h1, n_b, n_rsv:GOPC) %>%
  arrange(season)

# Loop through seasons and compile:
metrics_list = ILI_metrics_list = vector('list', length = length(unique(hk_dat$season)))
for (i in 1:length(unique(hk_dat$season))) {
  yr <- c(unique(hk_dat$season))[i]
  
  metrics_list[[i]] <- hk_dat %>%
    filter(season == yr) %>%
    # mutate(n_h1 = n_h1 / n_T,
    #        n_b = n_b / n_T,
    #        n_rsv = n_rsv / n_T) %>%
    select(-c(n_T, GOPC)) %>%
    pivot_longer(n_h1:n_rsv, names_to = 'vir', values_to = 'val') %>%
    group_by(vir) %>%
    summarise(ar = sum(val, na.rm = TRUE),
              pt = which.max(val) + 45) %>%
    mutate(season = yr)
  
  ILI_metrics_list[[i]] <- hk_dat %>%
    filter(season == yr) %>%
    mutate(ili_h1 = GOPC * (n_h1 / n_T),
           ili_b = GOPC * (n_b / n_T),
           ili_rsv = GOPC * (n_rsv / n_T)) %>%
    select(-c(n_h1:n_T)) %>%
    pivot_longer(ili_h1:ili_rsv, names_to = 'vir', values_to = 'val') %>%
    group_by(vir) %>%
    summarise(ar_all = sum(GOPC, na.rm = TRUE),
              ar_spec = sum(val, na.rm = TRUE),
              pt_spec = which.max(val) + 45) %>%
    mutate(season = yr)
}
rm(i, yr)

metrics <- bind_rows(metrics_list)
ILI_metrics <- bind_rows(ILI_metrics_list)
rm(metrics_list, ILI_metrics_list)

# Calculate peak timing difference:
metrics <- metrics %>%
  select(-ar) %>%
  pivot_wider(names_from = vir, values_from = pt) %>%
  mutate(n_h1 = n_h1 - n_rsv,
         n_b = n_b - n_rsv) %>%
  select(-n_rsv) %>%
  pivot_longer(c(n_h1, n_b), names_to = 'vir', values_to = 'pt_diff') %>%
  right_join(metrics, by = c('vir', 'season')) %>%
  arrange(season) %>%
  select(season:vir, ar:pt, pt_diff)

# Plot "realistic" values:
metrics_plot <- metrics %>%
  pivot_longer(ar:pt_diff, names_to = 'metric', values_to = 'val')

p1 <- ggplot(data = metrics_plot, aes(x = vir, y = val, fill = metric)) + geom_violin() +
  facet_wrap(~ metric, scales = 'free_y') + theme_classic() +
  labs(x = 'Virus', y = 'Value') + scale_fill_brewer(palette = 'Set1')
print(p1)

metrics_plot %>%
  group_by(vir, metric) %>%
  summarise(min = min(val),
            max = max(val),
            mean = mean(val),
            median = median(val)) %>%
  print()
# for flu_H1: AR between 2000 and 16000, PT between 55 and 71
# for flu_B: AR between 2600 and 14500, PT between 60 and 73
# for RSV: AR between 3000 to 8000 (smaller mean but larger median than either flu), PT between 56 and 90
# H1 up to 26 weeks earlier and 1 later (mean/median -10.4/-7); B up to 26 weeks earlier and 17 later (mean/median -3.8/3)

# in terms of positivity rates, both flus have wider range (smaller mins and larger maxes) than RSV
# using positivity rates, both flus always before RSV by at least 4 weeks

# Same, but using ILI+:
metrics_plot <- ILI_metrics %>%
  pivot_longer(ar_all:pt_spec, names_to = 'metric', values_to = 'val')

p2 <- ggplot(data = metrics_plot, aes(x = vir, y = val, fill = metric)) + geom_violin() +
  facet_wrap(~ metric, scales = 'free_y') + theme_classic() +
  labs(x = 'Virus', y = 'Value') + scale_fill_brewer(palette = 'Set1')
print(p2)

# Compare peak timing values for viral vs. syn+ data:
metrics_comp <- metrics %>%
  mutate(vir = str_remove(vir, 'n_')) %>%
  left_join(ILI_metrics %>%
              mutate(vir = str_remove(vir, 'ili_')),
            by = c('vir', 'season')) %>%
  select(season:vir, pt, pt_spec) %>%
  mutate(pt_diff = pt - pt_spec)

summary(metrics_comp$pt)
summary(metrics_comp$pt_spec)
summary(metrics_comp$pt_diff)
# pretty similar; biggest differences are B in 14-15 and H1 in 16-17 (10 weeks earlier in vir than ILI+), but both are small/not very peaky

# Check syn+ attack rates:
ILI_metrics <- ILI_metrics %>%
  select(vir, season, ar_spec) %>%
  pivot_wider(names_from = vir, values_from = ar_spec)

summary(ILI_metrics$ili_h1) # 2-17%, mean/median ~9.8/9.4
summary(ILI_metrics$ili_b) # 2-17%, mean/meidan ~9.6/10.8
summary(ILI_metrics$ili_rsv) # 5-11%, mean/median ~7.8/8.1

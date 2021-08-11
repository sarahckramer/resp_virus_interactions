# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates and peak timings for flu and RSV outbreaks
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in necessary functions:
source('src/functions/functions_flu_RSV.R')

# Read in data:
fr_dat <- read_rds('data/formatted/GROG_pop_vir_ari_dat_2003-4_2013-14.rds')

# Loop through seasons and compile:
metrics_list = ARI_metrics_list = vector('list', length = length(2006:2014))

for (i in 1:length(2006:2014)) {
  yr <- c(2006:2014)[i]
  
  formatted_dat <- prepare_data('flu_A', 'rsv', yr, fr_dat)
  dat_full <- formatted_dat[[1]]
  
  metrics_list[[i]] <- dat_full %>%
    select(week_no, virus, n_samp:i_ARI, time) %>%
    group_by(virus) %>%
    summarise(ar = sum(n_pos),
              pt = which.max(n_pos) + 39) %>%
    mutate(season = yr)
  
  ARI_metrics_list[[i]] <- dat_full %>%
    select(week_no, virus, n_samp:i_ARI, time) %>%
    mutate(prop_pos = n_pos / n_samp,
           i_ARI_spec = i_ARI * prop_pos) %>%
    group_by(virus) %>%
    summarise(ar_all = sum(i_ARI) * 100,
              ar_spec = sum(i_ARI_spec) * 100,
              pt_spec = which.max(i_ARI_spec) + 39) %>%
    mutate(season = yr)
}

metrics <- bind_rows(metrics_list)
ARI_metrics <- bind_rows(ARI_metrics_list)

# Calculate peak timing difference:
metrics <- metrics %>%
  select(-ar) %>%
  pivot_wider(names_from = virus, values_from = pt) %>%
  mutate(flu_A = flu_A - rsv,
         flu_B = flu_B - rsv) %>%
  select(-rsv) %>%
  pivot_longer(flu_A:flu_B, names_to = 'virus', values_to = 'pt_diff') %>%
  right_join(metrics, by = c('virus', 'season')) %>%
  arrange(season) %>%
  select(season:virus, ar:pt, pt_diff)

# Remove pandemic:
metrics <- metrics %>% filter(season != 2010)
ARI_metrics <- ARI_metrics %>% filter(season != 2010)

# Remove relevant flu_B seasons:
metrics <- metrics %>% filter(!(virus == 'flu_B' & season %in% c(2007, 2012, 2014)))
ARI_metrics <- ARI_metrics %>% filter(!(virus == 'flu_B' & season %in% c(2007, 2012, 2014)))

# Plot "realistic" values:
metrics_plot <- metrics %>%
  pivot_longer(ar:pt_diff, names_to = 'metric', values_to = 'val')

p1 <- ggplot(data = metrics_plot, aes(x = virus, y = val, fill = metric)) + geom_violin() +
  facet_wrap(~ metric, scales = 'free_y') + theme_classic() +
  labs(x = 'Virus', y = 'Value') + scale_fill_brewer(palette = 'Set1')
print(p1)

metrics_plot %>%
  group_by(virus, metric) %>%
  summarise(min = min(val),
            max = max(val),
            mean = mean(val),
            median = median(val)) %>%
  print()
# for flu: AR between 200 and 1400, PT between 53-63
# for RSV: AR between 100 to 400, PT between 46 and 55
# flu should be anywhere from 0 to 13 weeks later

# Compare peak timing values for viral vs. syn+ data:
metrics_comp <- metrics %>%
  left_join(ARI_metrics, by = c('virus', 'season')) %>%
  select(season:virus, pt, pt_spec) %>%
  mutate(pt_diff = pt - pt_spec)

summary(metrics_comp$pt)
summary(metrics_comp$pt_spec)
summary(metrics_comp$pt_diff)
# overall very similar; flu_B in 2008 has a larger difference, but vir data aren't very peaky

# Check syn+ attack rates (to assess appropriateness of rho = 0.5):
ARI_metrics <- ARI_metrics %>%
  select(virus, season, ar_spec) %>%
  pivot_wider(names_from = virus, values_from = ar_spec)

summary(ARI_metrics$flu_A) # 3-10%, mean/median ~8/8.6
summary(ARI_metrics$flu_B) # 3-10%, mean/meidan 6/5
summary(ARI_metrics$rsv) # 1.5-4%, mean/median ~2.7
# flu is broadly consistent with a reporting rate of about 50%, since it seems to infect 5-20% each year
# unclear what actual rates of RSV are; however, if similar to influenza, then reporting rate of 0.2-0.3 might be
# more appropriate

# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# Plot ARI/flu/RSV incidence and positivity by age group
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in data:
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
  select(-c(pop_eff, n_ARI)) %>%
  filter(epiyear >= 2006)

# ---------------------------------------------------------------------------------------------------------------------

pdf('results/plots/patterns_by_age.pdf', width = 14, height = 9)

# Plots to explore dynamics by age:
fr_dat %>% group_by(agecat) %>% summarise(mean_pop = mean(pop_tot))
# ages 0-4 make up ~6% of the population
# so this is the age group that we expect reporting to come from at least, and maybe also where we expect the
# most transmissible cases to occur (assuming asymptomatic cases are not or less likely to be transmissible)

p1 <- ggplot(data = fr_dat, aes(x = week_no, y = i_ARI, group = virus)) + geom_line() +
  facet_grid(agecat ~ epiyear) + theme_classic() + labs(x = 'Week', y = 'ARI Incidence')
print(p1)
# ARI tends to have double-peak in younger age groups; still there to some extent in older but overall flatter
# appearance and much lower rates

p2 <- ggplot(data = fr_dat, aes(x = week_no, y = n_pos, group = virus, col = virus)) + geom_line() +
  facet_grid(agecat ~ epiyear, scales = 'free_y') + theme_classic() + labs(x = 'Week', y = 'Positive Samples') +
  scale_color_brewer(palette = 'Set1')
p3 <- ggplot(data = fr_dat %>% filter(virus == 'rsv'), aes(x = week_no, y = n_pos)) + geom_line() +
  facet_grid(agecat ~ epiyear, scales = 'free_y') + theme_classic() + labs(x = 'Week', y = 'Positive RSV Samples')
print(p2)
print(p3)
# rates of everything are really low in 2006 - should we remove this year?
# flus tend to hit all age groups, but most reported RSV cases are in the 0-4 group; does this implay substantial
# asymptomatic infection among older age groups? probably; but to what extent do these age groups transmit?

p4 <- ggplot(data = fr_dat, aes(x = week_no, y = prop_pos, group = virus, col = virus)) + geom_line(lwd = 0.75) +
  facet_grid(agecat ~ epiyear) + theme_classic() + labs(x = 'Week', y = 'Proportion Positive') +
  scale_color_brewer(palette = 'Set1')
print(p4)

fr_dat <- fr_dat %>%
  mutate(i_ARI_plus = i_ARI * prop_pos)
p5 <- ggplot(data = fr_dat, aes(x = week_no, y = i_ARI_plus, group = virus, col = virus)) + geom_line() +
  facet_grid(agecat ~ epiyear) + theme_classic() + labs(x = 'Week', y = 'ARI+ Incidence') +
  scale_color_brewer(palette = 'Set1')
print(p5)
# estimated incidence of (medically-attended) flu also decreases in adults compared to young children

dev.off()

# ---------------------------------------------------------------------------------------------------------------------

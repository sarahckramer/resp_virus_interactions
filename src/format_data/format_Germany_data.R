# ---------------------------------------------------------------------------------------------------------------------
# Read in and format German data
# Source: https://github.com/robert-koch-institut/ARE-Konsultationsinzidenz
# ---------------------------------------------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(gridExtra)
library(ISOweek)
library(testthat)

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format syndromic data

# Read in ARI data:
# Note: Seasons seem to start in week 40
de_ari <- read_tsv('data/raw/germany/ARE-Konsultationsinzidenz.tsv')

# Limit to data of interest:
de_ari <- de_ari %>%
  filter(Bundesland_ID == 0,
         Altersgruppe == '00+') %>%
  select(Saison:Kalenderwoche, ARE_Konsultationsinzidenz)

# Get dates and seasons:
de_ari <- de_ari %>%
  mutate(year = str_sub(Kalenderwoche, 1, 4),
         week = str_sub(Kalenderwoche, 7, 8)) %>%
  mutate(season = paste0('s', str_sub(Saison, 3, 4), '-', str_sub(Saison, 6, 7))) %>%
  select(year:season, ARE_Konsultationsinzidenz) %>%
  mutate(year = as.numeric(year),
         week = as.numeric(week))

# Rename ARI column:
de_ari <- de_ari %>%
  rename('ari_rate' = 'ARE_Konsultationsinzidenz')

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format virologic data

# Read in data:
de_vir <- read_csv('data/raw/germany/02_positivity_complete.csv')

# Clean up:
de_vir <- de_vir %>%
  mutate(PATHOGEN = paste0(PATHOGEN, PATH.TYPE, PATH.STRAIN),
         PATHOGEN = str_remove_all(PATHOGEN, 'NA')) %>%
  select(YEAR:NoSamples, PATHOGEN, VALUE) %>%
  rename_with(tolower) %>%
  pivot_wider(names_from = pathogen,
              values_from = value) %>%
  rowwise() %>%
  mutate(tot_a = sum(FluAH3N2, FluAH1N1, na.rm = TRUE),
         tot_b = if_else(is.na(FluB), sum(FluBYam, FluBVic, na.rm = TRUE), FluB),
         flu = sum(tot_a, tot_b, FluXnottyped, na.rm = TRUE)) %>%
  ungroup() %>%
  select(year:nosamples, flu, tot_a:tot_b, RSV, Rhino) %>%
  rename_with(tolower) %>%
  rename('n_test' = 'nosamples') %>%
  drop_na()

# Check dates:
de_vir %>%
  mutate(date_check = ISOweek2date(paste0(year, '-W', week, '-1'))) %>%
  filter(date != date_check) %>%
  nrow() %>%
  expect_equal(0)

# ---------------------------------------------------------------------------------------------------------------------

# Combine and format for model
de_dat <- de_ari %>%
  inner_join(de_vir %>%
               mutate(week = as.numeric(week)),
             by = c('year', 'week')) %>%
  select(year:week, date, season:ari_rate, n_test:rhino)
expect_true(nrow(de_dat) == nrow(de_vir))
rm(de_ari, de_vir)

# Remove pandemic seasons:
de_dat <- de_dat %>%
  filter(!(season %in% c('s19-20', 's20-21', 's21-22', 's22-23')))

# Get population sizes:
# Source: https://www.destatis.de/EN/Themes/Society-Environment/Population/Current-Population/Tables/liste-current-population.html
# for 2014-2015 season, want 2014 Q4
pop_dat <- read_csv2('data/raw/germany/pop_dat_germany.csv') %>%
  mutate(pop = pop * 1000)

de_dat <- de_dat %>%
  left_join(pop_dat, by = 'season')
rm(pop_dat)

# Add Time column:
de_dat <- de_dat %>%
  split(~ .$season) %>%
  lapply(function(ix) {
    ix %>% mutate(time = 1:nrow(ix))
  }) %>%
  bind_rows()

# Rename columns for use in model:
de_dat <- de_dat %>%
  select(time, date, year:week, season, pop, n_test, flu:rhino, ari_rate) %>%
  rename('n_T' = n_test,
         'n_P1' = flu,
         'n_P2' = rsv,
         'i_ARI' = ari_rate) %>%
  mutate(i_ARI = i_ARI / 100000)

# Replace dips in ARI around New Year's with NA:
de_dat <- de_dat %>%
  mutate(i_ARI = if_else(week %in% 52:53 | (week == 1 & year %in% c(2015, 2019)), NA, i_ARI))

# Write data to file:
write_csv(de_dat, 'data/formatted/dat_germany.csv')

# ---------------------------------------------------------------------------------------------------------------------

# Visualize data

p1 <- ggplot(de_dat, aes(x = date, y = i_ARI)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'ARI Rate', title = 'Germany ARI')

de_dat_long <- de_dat %>%
  mutate(perc_flu = n_P1 / n_T * 100,
         perc_rsv = n_P2 / n_T * 100,
         perc_rhino = rhino / n_T * 100,
         perc_fluA = tot_a / n_T * 100,
         perc_fluB = tot_b / n_T * 100) %>%
  select(time:date, season, perc_flu:perc_fluB) %>%
  pivot_longer(perc_flu:perc_fluB,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = if_else(virus == 'perc_flu', 'Influenza', virus),
         virus = if_else(virus == 'perc_rsv', 'RSV', virus),
         virus = if_else(virus == 'perc_rhino', 'Rhinovirus', virus),
         virus = if_else(virus == 'perc_fluA', 'Influenza (A)', virus),
         virus = if_else(virus == 'perc_fluB', 'Influenza (B)', virus))

p2 <- ggplot(de_dat_long %>% filter(virus %in% c('Influenza', 'RSV')),
             aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus', title = 'Germany Virologic') +
  scale_color_brewer(palette = 'Set1')

p3 <- ggplot(de_dat_long %>% filter(virus %in% c('Influenza (A)', 'Influenza (B)')),
       aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus', title = 'Germany Virologic') +
  scale_color_brewer(palette = 'Set2')

p4 <- ggplot(de_dat_long %>%
               filter(virus %in% c('Influenza', 'RSV', 'Rhinovirus')) %>%
               mutate(virus = factor(virus, levels = c('Influenza', 'RSV', 'Rhinovirus'))),
             aes(x = time, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ season) +
  labs(x = 'Time (Weeks)', y = 'Proportion Positive', col = 'Virus', title = 'Germany Virologic (By Season)')

p5 <- ggplot(de_dat_long %>%
               filter(virus %in% c('Influenza', 'RSV', 'Rhinovirus')) %>%
               mutate(virus = factor(virus, levels = c('Influenza', 'RSV', 'Rhinovirus'))),
             aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus', title = 'Germany Virologic') +
  scale_color_brewer(palette = 'Set1')

# Save plots to file:
pdf('results/plots/data_Germany.pdf', width = 10, height = 7)
grid.arrange(p1, p2, ncol = 1)
grid.arrange(p5, p3, ncol = 1)
print(p4)
dev.off()

# Clean up:
rm(list = ls())

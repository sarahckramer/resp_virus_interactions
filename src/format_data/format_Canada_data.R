# ---------------------------------------------------------------------------------------------------------------------
# Read in and format Canadian data
# Source: https://search.open.canada.ca/opendata/?od-search-portal=Open%20Data&search_text=fluwatch
# ---------------------------------------------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(gridExtra)
library(ISOweek)
library(testthat)

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format syndromic data

# Read in ILI data (Canada.ca):
# Note: Seasons seem to be defined as starting in week 35 here
file_list <- list.files(path = 'data/raw/canada/', pattern = 'consultation', full.names = TRUE)

dat_list <- vector('list', length = length(file_list))
for (i in 1:length(dat_list)) {
  dat_list[[i]] <- read_csv2(file_list[i])
}
rm(i)

# Format all in same way and join:
dat_list <- lapply(dat_list, function(ix) {
  ix <- ix[, 1:2]
  names(ix) <- c('season_week', 'ili_rate')
  ix
})
dat_list <- lapply(dat_list, drop_na)
can_ili <- bind_rows(dat_list)
rm(file_list, dat_list)

# Get dates:
can_ili <- can_ili %>%
  mutate(year = str_sub(season_week, 1, 4),
         week = str_sub(season_week, 5, 6),
         date = ISOweek2date(paste0(year, '-W', week, '-1'))) %>%
  select(year:date, ili_rate) %>%
  mutate(year = as.numeric(year),
         week = as.numeric(week))

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format virologic data

# Read in data:
file_list <- list.files(path = 'data/raw/canada/', pattern = 'detection-', full.names = TRUE)

dat_list <- vector('list', length = length(file_list))
for (i in 1:length(dat_list)) {
  dat_list[[i]] <- read_csv2(file_list[i], col_types = 'cddddddddddddddddddddd')
}
rm(i)

# Join all years:
can_vir <- bind_rows(dat_list)
rm(file_list, dat_list)

# Clean up:
can_vir <- can_vir %>%
  select(`Influenza Season and Epidemiological Week`, `Total Flu tests`:`RSV positive`, `Total Rhinovirus`:`Positive Rhinovirus`) %>%
  rename('season_week' = 1,
         'n_test_flu' = 2,
         'h1_09' = 3,
         'h1' = 4,
         'h3' = 5,
         'a_unsubtyped' = 6,
         'tot_a' = 7,
         'tot_b' = 8,
         'n_test_rsv' = 9,
         'rsv' = 10,
         'n_test_rhino' = 11,
         'rhino' = 12) %>%
  mutate(flu = tot_a + tot_b,
         h1 = h1_09 + h1) %>%
  select(season_week:n_test_flu, flu, tot_a:tot_b, n_test_rsv:rhino) %>%
  drop_na()

# Get dates:
can_vir <- can_vir %>%
  mutate(year = str_sub(season_week, 1, 4),
         week = str_sub(season_week, 5, 6),
         date = ISOweek2date(paste0(year, '-W', week, '-1'))) %>%
  select(year:date, n_test_flu:rhino) %>%
  mutate(year = as.numeric(year),
         week = as.numeric(week))

# ---------------------------------------------------------------------------------------------------------------------

# Combine and format for model

# Combine syndromic and virologic data:
can_dat <- can_ili %>%
  inner_join(can_vir, by = c('year', 'week', 'date'))
expect_true(nrow(can_dat) == nrow(can_vir))
rm(can_ili, can_vir)

# Remove 2009 pandemic:
can_dat <- can_dat %>%
  filter(year > 2010 | (year == 2010 & week >= 35))

# Label with northern-hemisphere seasons:
start_week <- 35
can_dat <- can_dat %>%
  mutate(season = NA) %>%
  mutate(season = ifelse((year == 2010 & week >= start_week) | (year == 2011 & week < start_week), 's10-11', season),
         season = ifelse((year == 2011 & week >= start_week) | (year == 2012 & week < start_week), 's11-12', season),
         season = ifelse((year == 2012 & week >= start_week) | (year == 2013 & week < start_week), 's12-13', season),
         season = ifelse((year == 2013 & week >= start_week) | (year == 2014 & week < start_week), 's13-14', season))
rm(start_week)

# Check for NAs:
can_dat %>% pull(season) %>% is.na() %>% any() %>% expect_false()

# Get population sizes:
# Source: https://www.statcan.gc.ca/en/subjects-start/population_and_demography
# for 2010-2011 season, want: Q1 2011
pop_dat <- read_csv('data/raw/canada/pop_dat_Canada.csv', skip = 8, n_max = 15) %>%
  filter(Geography == 'Canada') %>%
  select(contains('Q1')) %>%
  rename_with(~ str_sub(., 4)) %>%
  mutate(`2010` = str_remove_all(`2010`, ','),
         `2010` = as.numeric(`2010`))

pop_dat <- pop_dat %>%
  pivot_longer(cols = everything(), names_to = 'season', values_to = 'pop') %>%
  mutate(season = as.numeric(season),
         season = paste0('s', str_sub(season - 1, 3, 4), '-', str_sub(season, 3, 4)))

can_dat <- can_dat %>%
  left_join(pop_dat, by = 'season')
rm(pop_dat)

# Add Time column:
can_dat <- can_dat %>%
  group_by(season) %>%
  mutate(time = 1:52) %>%
  ungroup()

# Rename columns for use in model:
can_dat <- can_dat %>%
  select(time, date, year:week, season, pop, n_test_flu, n_test_rsv, n_test_rhino, flu:tot_b, rsv, rhino, ili_rate) %>%
  rename('n_T1' = n_test_flu,
         'n_T2' = n_test_rsv,
         'n_P1' = flu,
         'n_P2' = rsv,
         'i_ILI' = 'ili_rate') %>%
  mutate(i_ILI = i_ILI / 1000)

# Write data to file:
write_csv(can_dat, 'data/formatted/dat_canada.csv')

# ---------------------------------------------------------------------------------------------------------------------

# Visualize data

p1 <- ggplot(can_dat, aes(x = date, y = i_ILI)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'ILI Rate', title = 'Canada ILI')

can_dat_long <- can_dat %>%
  mutate(perc_flu = n_P1 / n_T1 * 100,
         perc_rsv = n_P2 / n_T2 * 100,
         perc_rhino = rhino / n_test_rhino * 100) %>%
  select(time:date, season, perc_flu:perc_rhino) %>%
  pivot_longer(perc_flu:perc_rhino,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = if_else(virus == 'perc_flu', 'Influenza', virus),
         virus = if_else(virus == 'perc_rsv', 'RSV', virus),
         virus = if_else(virus == 'perc_rhino', 'Rhinovirus', virus))

p2 <- ggplot(can_dat_long %>% filter(virus %in% c('Influenza', 'RSV')),
             aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus', title = 'Canada Virologic') +
  scale_color_brewer(palette = 'Set1')

p3 <- ggplot(can_dat_long, aes(x = time, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ season) +
  labs(x = 'Time (Weeks)', y = 'Proportion Positive', col = 'Virus', title = 'Canada Virologic (By Season)')

p4 <- ggplot(can_dat_long,
             aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  theme_classic() +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus', title = 'Canada Virologic') +
  scale_color_brewer(palette = 'Set1')

# Save plots to file:
pdf('results/plots/data_Canada.pdf', width = 10, height = 7)
grid.arrange(p1, p2, ncol = 1)
print(p3)
print(p4)
dev.off()

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Read in and format United States data
# Source (influenza): https://www.cdc.gov/flu/weekly/fluviewinteractive.htm
# Source (RSV): https://github.com/JianiC/RSV_flu/tree/main (originally NREVSS)
# ---------------------------------------------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(ISOweek)
library(testthat)

# ---------------------------------------------------------------------------------------------------------------------

# Read and format syndromic data

# Read in data:
us_ili <- read_csv('data/raw/us/FluViewPhase2Data/ILINet.csv', skip = 1)

# Keep only columns/seasons of interest:
us_ili <- us_ili %>%
  select(REGION:WEEK, ILITOTAL) %>%
  # select(REGION:WEEK, `%UNWEIGHTED ILI`) %>%
  # rename('perc_ili' = `%UNWEIGHTED ILI`) %>%
  rename_with(tolower) %>%
  filter(year > 2009 & (year > 2010 | (year == 2010 & week >= 40))) %>%
  filter(year < 2020 & (year < 2019 | (year == 2019 & week < 20)))

# Get dates:
us_ili <- us_ili %>%
  mutate(week_temp = str_pad(week, 2, 'left', 0)) %>%
  mutate(date = ISOweek2date(paste0(year, '-W', week_temp, '-1'))) %>%
  select(region:week, date, ilitotal)
  # select(region:week, date, perc_ili)

# # Get ILI per single consultation:
# us_ili <- us_ili %>%
#   mutate(i_ILI = perc_ili / 100) %>%
#   select(-perc_ili)

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format virologic data (influenza)

# Read in all data:
us_vir_pre2015 <- read_csv('data/raw/us/FluViewPhase2Data/WHO_NREVSS_Combined_prior_to_2015_16.csv', skip = 1)
us_vir_clin_post2015 <- read_csv('data/raw/us/FluViewPhase2Data/WHO_NREVSS_Clinical_Labs.csv', skip = 1)
# us_vir_ph_post2015 <- read_csv('data/raw/us/FluViewPhase2Data/WHO_NREVSS_Public_Health_Labs.csv', skip = 1)

# Prepare data for combining (pre 2015):
us_vir_pre2015 <- us_vir_pre2015 %>%
  select(REGION:H3N2v) %>%
  rename_with(tolower) %>%
  mutate(tot_a = `a (2009 h1n1)` + `a (h1)` + `a (h3)` + `a (subtyping not performed)` + `a (unable to subtype)` + `h3n2v`,
         tot_b = b,
         flu = tot_a + tot_b) %>%
  select(region:`percent positive`, flu, tot_a:tot_b) %>%
  rename('n_test_flu' = `total specimens`) %>%
  mutate(check = flu / n_test_flu * 100)

expect_true(
  us_vir_pre2015 %>%
    filter(round(`percent positive`, 1) != round(check, 1)) %>%
    nrow() == 0
)

us_vir_pre2015 <- us_vir_pre2015 %>%
  select(region:n_test_flu, flu:tot_b)

# # Prepare data for combining (post 2015, public health labs):
# us_vir_ph_post2015 <- us_vir_ph_post2015 %>%
#   select(REGION:H3N2v) %>%
#   rename_with(tolower) %>%
#   mutate(tot_a = `a (2009 h1n1)` + `a (h3)` + `a (subtyping not performed)` + h3n2v,
#          tot_b = b + bvic + byam,
#          flu = tot_a + tot_b) %>%
#   rename('n_test_flu' = `total specimens`) %>%
#   select(region:n_test_flu, flu, tot_a:tot_b)

# Prepare data for combining (post 2015, clinical labs):
us_vir_clin_post2015 <- us_vir_clin_post2015 %>%
  select(REGION:`PERCENT POSITIVE`) %>%
  rename_with(tolower) %>%
  rename('n_test_flu' = `total specimens`,
         'tot_a' = `total a`,
         'tot_b' = `total b`) %>%
  mutate(flu = tot_a + tot_b,
         check = flu / n_test_flu * 100)

expect_true(
  us_vir_clin_post2015 %>%
    filter(round(`percent positive`, 1) != round(check, 1)) %>%
    nrow() == 0
)

us_vir_clin_post2015 <- us_vir_clin_post2015 %>%
  select(region:n_test_flu, flu, tot_a:tot_b)

# # Combine post-2015 data:
# us_vir_post2015 <- us_vir_clin_post2015 %>%
#   inner_join(us_vir_ph_post2015,
#              by = c('region', 'year', 'week'))
# expect_true(nrow(us_vir_post2015) == nrow(us_vir_clin_post2015))
# 
# us_vir_post2015 <- us_vir_post2015 %>%
#   mutate(n_test_flu = n_test_flu.x + n_test_flu.y,
#          flu = flu.x + flu.y,
#          tot_a = tot_a.x + tot_a.y,
#          tot_b = tot_b.x + tot_b.y) %>%
#   select(region:week, n_test_flu:tot_b)
# rm(us_vir_clin_post2015)

# Combine data from before and after 2015:
us_vir_flu <- bind_rows(us_vir_pre2015,
                        us_vir_clin_post2015)
rm(us_vir_pre2015, us_vir_clin_post2015)

# Keep only seasons of interest:
us_vir_flu <- us_vir_flu %>%
  filter(year > 2009 & (year > 2010 | (year == 2010 & week >= 40))) %>%
  filter(year < 2020 & (year < 2019 | (year == 2019 & week < 20)))

# Get dates:
us_vir_flu <- us_vir_flu %>%
  mutate(week_temp = str_pad(week, 2, 'left', 0)) %>%
  mutate(date = ISOweek2date(paste0(year, '-W', week_temp, '-1'))) %>%
  select(region:week, date, n_test_flu:tot_b)

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format virologic data (RSV)

# Read in data:
us_vir_rsv <- read_csv('data/raw/us/flu_RSV_by_HHS.csv')

# Keep only columns/seasons of interest:
us_vir_rsv <- us_vir_rsv %>%
  select(HHS_REGION, YEAR:WEEK, RSVtest:RSVpos, flutest:fluBpos) %>%
  rename_with(tolower) %>%
  filter(year > 2009 & (year > 2010 | (year == 2010 & week >= 40))) %>%
  filter(year < 2020 & (year < 2019 | (year == 2019 & week < 20)))

# Format for joining with virologic data:
us_vir_rsv <- us_vir_rsv %>%
  rename('region' = 'hhs_region',
         'n_test_rsv' = 'rsvtest',
         'rsv' = 'rsvpos') %>%
  mutate(region = paste('Region', region))

# # Check flu data:
# us_vir_check <- us_vir_ph_post2015 %>%
#   inner_join(us_vir_rsv,
#              by = c('region', 'year', 'week'))
# 
# expect_true(us_vir_check %>% filter(n_test_flu != flutest) %>% nrow == 0)
# expect_true(us_vir_check %>% filter(flu != flupos) %>% nrow == 0)
# expect_true(us_vir_check %>% filter(tot_a != fluapos) %>% nrow == 0)
# expect_true(us_vir_check %>% filter(tot_b != flubpos) %>% nrow == 0)
# rm(us_vir_check, us_vir_ph_post2015)

# Join with flu virologic data:
us_vir <- us_vir_flu %>%
  full_join(us_vir_rsv %>%
              select(region:rsv),
            by = c('region', 'year', 'week'))
expect_true(nrow(us_vir) == nrow(us_vir_flu))
rm(us_vir_flu, us_vir_rsv)

# ---------------------------------------------------------------------------------------------------------------------

# Combine syndromic and virologic data
us_dat <- us_ili %>%
  inner_join(us_vir,
             by = c('region', 'year', 'week', 'date'))
expect_true(nrow(us_dat) == nrow(us_ili))
expect_true(nrow(us_dat) == nrow(us_vir))
rm(us_ili, us_vir)

# ---------------------------------------------------------------------------------------------------------------------

# Format for model input

# Label with northern-hemisphere seasons:
start_week <- 40
us_dat <- us_dat %>%
  mutate(season = NA) %>%
  mutate(season = ifelse((year == 2010 & week >= start_week) | (year == 2011 & week < start_week), 's10-11', season),
         season = ifelse((year == 2011 & week >= start_week) | (year == 2012 & week < start_week), 's11-12', season),
         season = ifelse((year == 2012 & week >= start_week) | (year == 2013 & week < start_week), 's12-13', season),
         season = ifelse((year == 2013 & week >= start_week) | (year == 2014 & week < start_week), 's13-14', season),
         season = ifelse((year == 2014 & week >= start_week) | (year == 2015 & week < start_week), 's14-15', season),
         season = ifelse((year == 2015 & week >= start_week) | (year == 2016 & week < start_week), 's15-16', season),
         season = ifelse((year == 2016 & week >= start_week) | (year == 2017 & week < start_week), 's16-17', season),
         season = ifelse((year == 2017 & week >= start_week) | (year == 2018 & week < start_week), 's17-18', season),
         season = ifelse((year == 2018 & week >= start_week) | (year == 2019 & week < start_week), 's18-19', season))
rm(start_week)

# Check for NAs:
us_dat %>% pull(season) %>% is.na() %>% any() %>% expect_false()

# Get population sizes:
# Source: https://www.census.gov/data/tables/time-series/demo/popest/2010s-counties-total.html
# States by HHS region: https://www.hhs.gov/about/agencies/iea/regional-offices/index.html
# for 2010-11 season, want average of 2010 and 2011 data
pop_dat <- read_csv('data/raw/us/pop_dat_us.csv')
pop_dat <- pop_dat %>%
  filter(STNAME == CTYNAME) %>%
  select(STNAME, POPESTIMATE2010:POPESTIMATE2019) %>%
  distinct() %>%
  pivot_longer(-STNAME,
               names_to = 'year',
               values_to = 'pop') %>%
  mutate(year = str_remove(year, 'POPESTIMATE')) %>%
  pivot_wider(names_from = STNAME,
              values_from = pop)

# Sum to get HHS region totals:
pop_dat <- pop_dat %>%
  mutate(Region_1 = Connecticut + Maine + Massachusetts + `New Hampshire` + `Rhode Island` + Vermont,
         Region_2 = `New Jersey` + `New York`,
         Region_3 = Delaware + `District of Columbia` + Maryland + Pennsylvania + Virginia + `West Virginia`,
         Region_4 = Alabama + Florida + Georgia + Kentucky + Mississippi + `North Carolina` + `South Carolina` + Tennessee,
         Region_5 = Illinois + Indiana + Michigan + Minnesota + Ohio + Wisconsin,
         Region_6 = Arkansas + Louisiana + `New Mexico` + Oklahoma + Texas,
         Region_7 = Iowa + Kansas + Missouri + Nebraska,
         Region_8 = Colorado + Montana + `North Dakota` + `South Dakota` + Utah + Wyoming,
         Region_9 = Arizona + California + Hawaii + Nevada,
         Region_10 = Alaska + Idaho + Oregon + Washington) %>%
  # select(-c(Connecticut, Maine, Massachusetts, `New Hampshire`, `Rhode Island`, Vermont,
  #           `New Jersey`, `New York`,
  #           Delaware, `District of Columbia`, Maryland, Pennsylvania, Virginia, `West Virginia`,
  #           Alabama, Florida, Georgia, Kentucky, Mississippi, `North Carolina`, `South Carolina`, Tennessee,
  #           Illinois, Indiana, Michigan, Minnesota, Ohio, Wisconsin,
  #           Arkansas, Louisiana, `New Mexico`, Oklahoma, Texas,
  #           Iowa, Kansas, Missouri, Nebraska,
  #           Colorado, Montana, `North Dakota`, `South Dakota`, Utah, Wyoming,
  #           Arizona, California, Hawaii, Nevada,
  #           Alaska, Idaho, Oregon, Washington)) %>%
  select(year, Region_1:Region_10)
# Note: Missing pop data for territories in Regions 2 (Puerto Rico, Virgin Islands) and 9 (American Samoa,
# Commonwealth of the Northern Mariana Islands, Federated States of Micronesia, Guam, Marshall Islands, and Republic of Palau)

# Get seasonal averages:
pop_dat <- pop_dat %>%
  pivot_longer(-year,
               names_to = 'region',
               values_to = 'pop') %>%
  pivot_wider(names_from = year,
              values_from = pop) %>%
  mutate(`s10-11` = round((`2010` + `2011`) / 2),
         `s11-12` = round((`2011` + `2012`) / 2),
         `s12-13` = round((`2012` + `2013`) / 2),
         `s13-14` = round((`2013` + `2014`) / 2),
         `s14-15` = round((`2014` + `2015`) / 2),
         `s15-16` = round((`2015` + `2016`) / 2),
         `s16-17` = round((`2016` + `2017`) / 2),
         `s17-18` = round((`2017` + `2018`) / 2),
         `s18-19` = round((`2018` + `2019`) / 2)) %>%
  select(region, `s10-11`:`s18-19`) %>%
  pivot_longer(-region,
               names_to = 'season',
               values_to = 'pop')

# Join to influenza data:
nrow_orig <- nrow(us_dat)
us_dat <- us_dat %>%
  inner_join(pop_dat %>%
               mutate(region = str_replace(region, '_', ' ')),
             by = c('region', 'season'))
expect_true(nrow(us_dat) == nrow_orig)
rm(nrow_orig, pop_dat)

# Calculate ILI rate per population:
us_dat <- us_dat %>%
  mutate(i_ILI = ilitotal / pop)

# Add time column:
us_dat <- us_dat %>%
  group_by(region, season) %>%
  mutate(time = 1:length(season)) %>%
  ungroup() %>%
  select(time, region, date, year:week, season, pop, n_test_flu, n_test_rsv, flu:tot_b, rsv, i_ILI)

# Rename columns for use in model:
us_dat_mod <- us_dat %>%
  select(time:flu, rsv:i_ILI) %>%
  rename('n_T1' = 'n_test_flu',
         'n_T2' = 'n_test_rsv',
         'n_P1' = 'flu',
         'n_P2' = 'rsv')

# Separate by region:
us_dat_mod <- us_dat_mod %>%
  mutate(region = factor(region, levels = c('Region 1', 'Region 2', 'Region 3', 'Region 4', 'Region 5',
                                            'Region 6', 'Region 7', 'Region 8', 'Region 9', 'Region 10'))) %>%
  split(.$region)

# Write data to file:
write_rds(us_dat_mod, file = 'data/formatted/dat_us_byRegion.rds')

# ---------------------------------------------------------------------------------------------------------------------

# Visualize data
p1 <- ggplot(us_dat, aes(x = date, y = i_ILI * 1000)) +
  geom_line() +
  # geom_vline(xintercept = us_dat$date[us_dat$week == 40]) +
  facet_wrap(~ region, ncol = 2) +
  theme_classic() +
  labs(x = 'Date', y = 'ILI per 1000 Population')

us_dat_long <- us_dat %>%
  mutate(perc_flu = flu / n_test_flu * 100,
         perc_a = tot_a / n_test_flu * 100,
         perc_b = tot_b / n_test_flu * 100,
         perc_rsv = rsv / n_test_rsv * 100) %>%
  select(region, date, perc_flu:perc_rsv) %>%
  pivot_longer(perc_flu:perc_rsv,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = if_else(virus == 'perc_flu', 'Influenza (A + B)', virus),
         virus = if_else(virus == 'perc_a', 'Influenza (A)', virus),
         virus = if_else(virus == 'perc_b', 'Influenza (B)', virus),
         virus = if_else(virus == 'perc_rsv', 'RSV', virus))

p2 <- ggplot(us_dat_long %>% filter(virus %in% c('Influenza (A + B)', 'RSV')),
       aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  # geom_vline(xintercept = us_dat$date[us_dat$week == 40]) +
  facet_wrap(~ region, ncol = 2) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus') +
  scale_color_brewer(palette = 'Set1')
p3 <- ggplot(us_dat_long %>% filter(virus %in% c('Influenza (A)', 'Influenza (B)', 'RSV')),
             aes(x = date, y = perc_pos, group = virus, col = virus)) +
  geom_line() +
  facet_wrap(~ region, ncol = 2) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(x = 'Date', y = 'Percent Positive', col = 'Virus') +
  scale_color_brewer(palette = 'Set1')
rm(us_dat_long)

# Save plots to file:
pdf('results/plots/data_US.pdf', width = 12, height = 7)
print(p1)
print(p2)
print(p3)
dev.off()

rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Get and format climate (temperature, humidity, other?) data
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(GSODR)
library(tidyverse)
library(testthat)

# # List all available country names:
# get_inventory() %>%
#   pull(COUNTRY_NAME) %>%
#   table()

# Get Hong Kong data:
dat_hk <- get_GSOD(2012:2019, country = 'HONG KONG SAR CHINA')

# Get station nearest/in Santiago, Chile:
stat_id_chl <- nearest_stations(LAT = -33.45, LON = -70.67, distance = 5)

# Get Chile data:
dat_chile <- get_GSOD(2006:2019, station = stat_id_chl)

# Get French data:
stat_id_fr <- nearest_stations(LAT = 48.85, LON = 2.35, distance = 13)
dat_fr <- get_GSOD(2012:2019, station = stat_id_fr)

# Format data from all regions and calculate AH in g/m3:
dat_hk <- dat_hk %>%
  select(NAME, ISO2C, LATITUDE, LONGITUDE, YEARMODA:TEMP, MIN, MAX, RH, PRCP) %>%
  mutate(ES_CC = 0.611 * exp((2256000 / 461.52) * (1 / 273.15 - 1 / (TEMP + 273.15))),
         ES_ARM = 0.61094 * exp((17.625 * (TEMP)) / (TEMP + 243.04)),
         EA = (RH / 100) * ES_ARM,
         AH = (EA * 1000) * 2.16679 / (TEMP + 273.15)) %>%
  select(NAME:MAX, AH, RH:PRCP)
dat_chile <- dat_chile %>%
  select(NAME, ISO2C, LATITUDE, LONGITUDE, YEARMODA:TEMP, MIN, MAX, RH, PRCP) %>%
  mutate(ES_CC = 0.611 * exp((2256000 / 461.52) * (1 / 273.15 - 1 / (TEMP + 273.15))),
         ES_ARM = 0.61094 * exp((17.625 * (TEMP)) / (TEMP + 243.04)),
         EA = (RH / 100) * ES_ARM,
         AH = (EA * 1000) * 2.16679 / (TEMP + 273.15)) %>%
  select(NAME:MAX, AH, RH:PRCP)
dat_fr <- dat_fr %>%
  select(NAME, ISO2C, LATITUDE, LONGITUDE, YEARMODA:TEMP, MIN, MAX, RH, PRCP) %>%
  mutate(ES_CC = 0.611 * exp((2256000 / 461.52) * (1 / 273.15 - 1 / (TEMP + 273.15))),
         ES_ARM = 0.61094 * exp((17.625 * (TEMP)) / (TEMP + 243.04)),
         EA = (RH / 100) * ES_ARM,
         AH = (EA * 1000) * 2.16679 / (TEMP + 273.15)) %>%
  select(NAME:MAX, AH, RH:PRCP)

# Get week numbers:
dat_hk <- dat_hk %>%
  mutate(WEEK = as.numeric(strftime(YEARMODA, format = '%U'))) %>%
  group_by(YEAR) %>%
  mutate(YEAR_NEW = ifelse(WEEK == max(WEEK) & YEAR != 2016, YEAR + 1, YEAR),
         WEEK = ifelse(WEEK == max(WEEK) & YEAR != 2016, 0, WEEK)) %>%
  ungroup() %>%
  mutate(YEAR = YEAR_NEW) %>%
  select(-YEAR_NEW) %>%
  mutate(WEEK = ifelse(YEAR == 2017, WEEK, WEEK + 1)) %>%
  filter(YEAR < 2020 & YEAR >= 2013)

dat_chile <- dat_chile %>%
  mutate(WEEK = as.numeric(strftime(YEARMODA, format = '%V')),
         YEAR = ifelse(WEEK == 1 & MONTH == 12, YEAR + 1, YEAR),
         YEAR = ifelse((WEEK == 52 | WEEK == 53) & MONTH == 1, YEAR - 1, YEAR))

dat_fr <- dat_fr %>%
  mutate(WEEK = as.numeric(strftime(YEARMODA, format = '%U'))) %>%
  group_by(YEAR) %>%
  mutate(YEAR_NEW = ifelse(WEEK == max(WEEK) & YEAR != 2016, YEAR + 1, YEAR),
         WEEK = ifelse(WEEK == max(WEEK) & YEAR != 2016, 0, WEEK)) %>%
  ungroup() %>%
  mutate(YEAR = YEAR_NEW) %>%
  select(-YEAR_NEW) %>%
  mutate(WEEK = ifelse(YEAR == 2017, WEEK, WEEK + 1)) %>%
  filter(YEAR < 2020 & YEAR >= 2013)

# Check whether data are available for full weeks or not:
dat_hk %>%
  group_by(YEAR, WEEK) %>%
  summarise(days_w_data = length(AH)) %>%
  filter(days_w_data != 7)

dat_chile %>%
  group_by(YEAR, WEEK) %>%
  summarise(days_w_data = length(AH)) %>%
  filter(days_w_data != 7)
# A lot more missingness here - take average of just the available days?

dat_fr %>%
  group_by(YEAR, WEEK) %>%
  summarise(days_w_data = length(AH)) %>%
  filter(days_w_data != 7)
# a few weeks each in 2015 and 2018, but all still have 4-6 days of data

# Check that all weeks are present in data:
dat_hk %>%
  mutate(yearweek = paste(YEAR, WEEK, sep = '_')) %>%
  select(YEAR, yearweek) %>%
  unique() %>%
  pull(YEAR) %>% table()

dat_chile %>%
  mutate(yearweek = paste(YEAR, WEEK, sep = '_')) %>%
  select(YEAR, yearweek) %>%
  unique() %>%
  pull(YEAR) %>% table()
# 2006 missing 6 weeks; 2009 missing 3; 2010 missing 1; 2016 missing 16 (weeks 36-51); 2017 missing 1

dat_fr %>%
  mutate(yearweek = paste(YEAR, WEEK, sep = '_')) %>%
  select(YEAR, yearweek) %>%
  unique() %>%
  pull(YEAR) %>% table()

# Convert data to weekly:
dat_hk <- dat_hk %>%
  group_by(NAME, ISO2C, YEAR, WEEK) %>%
  summarise(date = min(YEARMODA),
            lat = unique(LATITUDE),
            long = unique(LONGITUDE),
            temp = mean(TEMP),
            min_temp = mean(MIN),
            max_temp = mean(MAX),
            ah = mean(AH),
            rh = mean(RH),
            prcp = mean(PRCP)) %>%
  rename('stn' = 'NAME',
         'iso2c' = 'ISO2C',
         'year' = 'YEAR',
         'week' = 'WEEK') %>%
  ungroup()
# dat_chile %>%
#   group_by(NAME, ISO2C, YEAR, WEEK) %>%
#   summarise(date = min(YEARMODA),
#             lat = unique(LATITUDE),
#             long = unique(LONGITUDE),
#             temp = mean(TEMP),
#             min_temp = mean(MIN),
#             max_temp = mean(MAX),
#             ah = mean(AH),
#             rh = mean(RH),
#             prcp = mean(PRCP)) %>%
#   rename('stn' = 'NAME',
#          'iso2c' = 'ISO2C',
#          'year' = 'YEAR',
#          'week' = 'WEEK') %>%
#   ungroup()
dat_fr <- dat_fr %>%
  group_by(NAME, ISO2C, YEAR, WEEK) %>%
  summarise(date = min(YEARMODA),
            lat = unique(LATITUDE),
            long = unique(LONGITUDE),
            temp = mean(TEMP),
            min_temp = mean(MIN),
            max_temp = mean(MAX),
            ah = mean(AH),
            rh = mean(RH),
            prcp = mean(PRCP)) %>%
  rename('stn' = 'NAME',
         'iso2c' = 'ISO2C',
         'year' = 'YEAR',
         'week' = 'WEEK') %>%
  ungroup()

# Write data to file:
write_csv(dat_hk, file = 'data/formatted/clim_dat_hk.csv')
write_csv(dat_fr, file = 'data/formatted/clim_dat_fr.csv')

# Plot data:
dat_hk_plot <- dat_hk %>%
  select(date, temp, ah:prcp) %>%
  rename('Temperature' = 'temp',
         'Abs. Humidity' = 'ah',
         'Rel. Humidity' = 'rh',
         'Precipitation' = 'prcp') %>%
  pivot_longer(-date, names_to = 'measure', values_to = 'val') %>%
  mutate(measure = factor(measure))
dat_hk_plot$measure <- factor(dat_hk_plot$measure, levels = levels(dat_hk_plot$measure)[c(4, 1, 3, 2)])

p_hk <- ggplot(data = dat_hk_plot, aes(x = date, y = val)) + geom_point() + geom_line() +
  theme_classic() + facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  labs(x = 'Year', y = '')
print(p_hk)

dat_fr_plot <- dat_fr %>%
  select(date, temp, ah:prcp) %>%
  rename('Temperature' = 'temp',
         'Abs. Humidity' = 'ah',
         'Rel. Humidity' = 'rh',
         'Precipitation' = 'prcp') %>%
  pivot_longer(-date, names_to = 'measure', values_to = 'val') %>%
  mutate(measure = factor(measure))
dat_fr_plot$measure <- factor(dat_fr_plot$measure, levels = levels(dat_fr_plot$measure)[c(4, 1, 3, 2)])

p_fr <- ggplot(data = dat_fr_plot, aes(x = date, y = val)) + geom_point() + geom_line() +
  theme_classic() + facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  labs(x = 'Year', y = '')
print(p_fr)

# Normalize data (as in Yaari et al.):
dat_hk_norm <- dat_hk %>%
  mutate(across(temp:prcp, ~ (.x - mean(.x)) / sd(.x)))
dat_hk_norm_CHECK <- dat_hk %>%
  mutate(across(temp:prcp, ~ scale(.x)[,1 ]))
expect_true(all.equal(dat_hk_norm, dat_hk_norm_CHECK))
rm(dat_hk_norm_CHECK)

dat_fr_norm <- dat_fr %>%
  mutate(across(temp:prcp, ~ (.x - mean(.x)) / sd(.x)))
dat_fr_norm_CHECK <- dat_fr %>%
  mutate(across(temp:prcp, ~ scale(.x)[,1 ]))
expect_true(all.equal(dat_fr_norm, dat_fr_norm_CHECK))
rm(dat_fr_norm_CHECK)

# Check that means ~ 0 and sd ~ 1:
dat_hk_norm %>%
  summarise(across(temp:prcp, mean))
dat_hk_norm %>%
  summarise(across(temp:prcp, var))

dat_fr_norm %>%
  summarise(across(temp:prcp, mean))
dat_fr_norm %>%
  summarise(across(temp:prcp, var))

# Write data to file:
write_csv(dat_hk_norm, file = 'data/formatted/clim_dat_hk_NORM.csv')
write_csv(dat_fr_norm, file = 'data/formatted/clim_dat_fr_NORM.csv')

# Plot data:
dat_hk_plot <- dat_hk_norm %>%
  select(date, temp, ah:prcp) %>%
  rename('Temperature' = 'temp',
         'Abs. Humidity' = 'ah',
         'Rel. Humidity' = 'rh',
         'Precipitation' = 'prcp') %>%
  pivot_longer(-date, names_to = 'measure', values_to = 'val') %>%
  mutate(measure = factor(measure))
dat_hk_plot$measure <- factor(dat_hk_plot$measure, levels = levels(dat_hk_plot$measure)[c(4, 1, 3, 2)])

p_hk_norm <- ggplot(data = dat_hk_plot, aes(x = date, y = val)) + geom_point() + geom_line() +
  theme_classic() + facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  labs(x = 'Year', y = '')
print(p_hk_norm)

dat_fr_plot <- dat_fr_norm %>%
  select(date, temp, ah:prcp) %>%
  rename('Temperature' = 'temp',
         'Abs. Humidity' = 'ah',
         'Rel. Humidity' = 'rh',
         'Precipitation' = 'prcp') %>%
  pivot_longer(-date, names_to = 'measure', values_to = 'val') %>%
  mutate(measure = factor(measure))
dat_fr_plot$measure <- factor(dat_fr_plot$measure, levels = levels(dat_fr_plot$measure)[c(4, 1, 3, 2)])

p_fr_norm <- ggplot(data = dat_fr_plot, aes(x = date, y = val)) + geom_point() + geom_line() +
  theme_classic() + facet_wrap(~ measure, scales = 'free_y', ncol = 2) +
  labs(x = 'Year', y = '')
print(p_fr_norm)

# Clean up:
rm(list = ls())

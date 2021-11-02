# ---------------------------------------------------------------------------------------------------------------------
# Format WHO (RSV and flu) data
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(ISOweek)
library(gridExtra)

# Read in positivity data:
dat_pos <- read_csv('data/raw/WHO_country_RSV_Flu_trends_data.csv')

# Data repeats?:
expect_true(all.equal(dat_pos[1:2202, ], dat_pos[2203:4404, ]))
dat_pos <- dat_pos[1:2202, ]

# Format:
dat_pos <- dat_pos %>%
  filter(Year > 2016) %>%
  rename('Country' = 'Countryname',
         'prop_flu' = 'Avg. Flu % positivity',
         'prop_rsv' = 'Avg. RSV % positivity') %>%
  mutate(prop_flu = prop_flu / 100,
         prop_rsv = prop_rsv / 100) %>%
  mutate(Date = ISOweek2date(paste0(Year, '-W', str_pad(Epiweek, width=2, side = 'left', pad = '0'), '-1')))

# Get number of samples:
dat_samp <- read_csv('data/raw/WHO_RSV_specimens_tested_data.csv')

# Format:
dat_samp <- dat_samp %>%
  filter(Year > 2016) %>%
  rename('Country' = 'Countryname',
         'n_samp' = 'Total Spec Tested')

# Join and convert to absolute number of positive tests (RSV):
expect_true(nrow(dat_pos) == nrow(dat_samp))
dat_who <- dat_pos %>%
  inner_join(dat_samp, by = c('Year', 'Epiweek', 'Country')) %>%
  mutate(n_rsv = n_samp * prop_rsv) %>%
  select(Year:Epiweek, Date, Country:prop_rsv, n_rsv, n_samp) %>%
  filter(Country != 'Canada')

# Explore where flu occurs before RSV:
p1 <- ggplot(data = dat_who %>% select(-c(n_rsv:n_samp)) %>% pivot_longer(prop_flu:prop_rsv, names_to = 'vir', values_to = 'val'),
             aes(x = Date, y = val, color = vir)) +
  geom_line() + theme_classic() + labs(x = 'Date', y = 'Prop. Positive', color = 'Virus') +
  scale_color_brewer(palette = 'Set1') + facet_wrap(~ Country)
print(p1)
# Primarily Thailand, which has a lot of overlap; Chile also has overlap, but RSV tends to begin earlier

# p2 <- ggplot(data = dat_who) + geom_line(aes(x = Date, y = n_samp)) + geom_line(aes(x = Date, y = n_rsv), col = 'coral') +
#   theme_classic() + labs(x = 'Date', y = '# Samples (Total and RSV)') + facet_wrap(~ Country, scales = 'free_y')
# print(p2)

# Compare flu positivity to FluNet data:
dat_flunet <- read_csv('data/raw/FluNet_Thailand.csv')
dat_flunet <- dat_flunet %>%
  select(Year:Week, SDATE, Country, SPEC_RECEIVED_NB:SPEC_PROCESSED_NB, INF_A, INF_B, ALL_INF) %>%
  rename(Epiweek = Week)

dat_who <- dat_who %>%
  filter(Country == 'Thailand') %>%
  left_join(dat_flunet, by = c('Year', 'Epiweek', 'Country'))

expect_true(all.equal(dat_who$SPEC_RECEIVED_NB, dat_who$SPEC_PROCESSED_NB))
expect_true(all.equal(dat_who$SDATE, dat_who$Date))

dat_who <- dat_who %>%
  select(!c(SDATE, SPEC_RECEIVED_NB)) %>%
  rename(n_fluA = INF_A,
         n_fluB = INF_B,
         n_flu = ALL_INF,
         n_samp_rsv = n_samp,
         n_samp_flu = SPEC_PROCESSED_NB) %>%
  mutate(prop_flu_check = n_flu / n_samp_flu,
         prop_flu_A = n_fluA / n_samp_flu,
         prop_flu_B = n_fluB / n_samp_flu) %>%
  select(Year:Country, prop_rsv, prop_flu, prop_flu_check:prop_flu_B, n_rsv, n_flu, n_fluA:n_fluB, n_samp_rsv:n_samp_flu)

p3 <- ggplot(data = dat_who) + geom_line(aes(x = Date, y = prop_flu)) + geom_line(aes(x = Date, y = prop_flu_check), col = 'coral') +
  theme_classic() + labs(x = 'Date', y = 'Prop. Flu')
print(p3)
# Not exactly the same, but similar; and FluNet data allows us to separate into fluA and fluB

dat_who <- dat_who %>%
  select(!prop_flu) %>%
  rename(prop_flu = prop_flu_check)

# Read in and format syndromic data:
dat_thailand <- read_csv('data/raw/FluID_Thailand.csv')

dat_thailand <- dat_thailand %>%
  select(Country:Week, ILI_CASES:ILI_OUTPATIENTS, ILI_OUT_RATE) %>%
  rename(Epiweek = Week)

dat_who <- dat_who %>%
  inner_join(dat_thailand, by = c('Year', 'Epiweek', 'Country'))

# Plot ILI and positivity data:
p4 <- ggplot(data = dat_who %>% select(Year:Country, prop_rsv, prop_flu_A:prop_flu_B) %>% pivot_longer(prop_rsv:prop_flu_B, names_to = 'vir', values_to = 'val'),
       aes(x = Date, y = val, color = vir)) +
  geom_line() + theme_classic() + labs(x = 'Date', y = 'Prop. Positive', color = 'Virus') +
  scale_color_brewer(palette = 'Set1')
p5 <- ggplot(data = dat_who) + geom_line(aes(x = Date, y = ILI_OUT_RATE)) + theme_classic() + labs(x = 'Date', y = 'ILI Rate')
p6 <- ggplot(data = dat_who %>%
               mutate(ili_plus_fluA = ILI_OUT_RATE * prop_flu_A, ili_plus_fluB = ILI_OUT_RATE * prop_flu_B, ili_plus_rsv = ILI_OUT_RATE * prop_rsv) %>%
               select(Year:Country, ili_plus_fluA:ili_plus_rsv) %>%
               pivot_longer(ili_plus_fluA:ili_plus_rsv, names_to = 'vir', values_to = 'val'),
             aes(x = Date, y = val, color = vir)) +
  geom_line() + theme_classic() + labs(x = 'Date', y = 'ILI+', color = 'Virus') +
  scale_color_brewer(palette = 'Set1')
grid.arrange(p4, p5, p6, ncol = 1)

# ---------------------------------------------------------------------------------------------------------------------
# Format WHO (___) data
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(ISOweek)

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

# Compare flu positivity to FluNet data:





# Read in syndromic data:




# Format:




# Plot:
p1 <- ggplot(data = dat_who %>% select(-c(n_rsv:n_samp)) %>% pivot_longer(prop_flu:prop_rsv, names_to = 'vir', values_to = 'val'),
             aes(x = Date, y = val, color = vir)) +
  geom_line() + theme_classic() + labs(x = 'Date', y = 'Prop. Positive', color = 'Virus') +
  scale_color_brewer(palette = 'Set1') + facet_wrap(~ Country)
p2 <- ggplot(data = dat_who) + geom_line(aes(x = Date, y = n_samp)) + geom_line(aes(x = Date, y = n_rsv), col = 'coral') +
  theme_classic() + labs(x = 'Date', y = '# Samples (Total and RSV)') + facet_wrap(~ Country, scales = 'free_y')
print(p1)
print(p2)

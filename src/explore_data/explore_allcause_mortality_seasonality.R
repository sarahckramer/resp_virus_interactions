# ---------------------------------------------------------------------------------------------------------------------
# Check whether seasonal trends in all-cause mortality exist in Hong Kong
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)

# Read in data:
mortality_dat <- read_csv('data/raw/hk_all-cause_mortality.csv')

# Format data:
mortality_dat <- mortality_dat %>%
  mutate(Yr = as.numeric(Yr),
         Month = as.numeric(Month),
         No_of_deaths = as.numeric(No_of_deaths)) %>%
  filter(!is.na(No_of_deaths)) %>%
  drop_na()

# Get average over all years:
avg_dat <- mortality_dat %>%
  group_by(Month) %>%
  summarise(No_of_deaths = median(No_of_deaths))

# Plot data:
p1 <- ggplot(data = mortality_dat, aes(x = Month, y = No_of_deaths, group = Yr, col = Yr)) +
  geom_point() + geom_line() + theme_classic() + labs(y = 'Number of Deaths', color = 'Year') +
  scale_x_continuous(breaks = 1:12) + scale_color_distiller(palette = 'RdYlGn')
print(p1)

p2 <- ggplot(data = mortality_dat) + geom_line(aes(x = Month, y = No_of_deaths, group = Yr, col = Yr)) +
  geom_line(data = avg_dat, aes(x = Month, y = No_of_deaths), lwd = 1.2) +
  theme_classic() + labs(y = 'Number of Deaths', color = 'Year') +
  scale_x_continuous(breaks = 1:12) + scale_color_distiller(palette = 'RdYlGn')
print(p2)

# Clean up:
rm(list = ls())

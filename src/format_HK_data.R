# ---------------------------------------------------------------------------------------------------------------------
# Format Hong Kong data
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Enter years for which data are available:
years <- 2014:2021

# Read in data:
file_list_flu <- list.files(path = 'data/raw/', pattern = 'flu', full.names = TRUE)
file_list_rsv <- list.files(path = 'data/raw/', pattern = 'rsv', full.names = TRUE)

data_list_flu = data_list_rsv = vector('list', length = length(file_list_flu))
for (i in 1:length(file_list_flu)) {
  data_list_flu[[i]] <- read_csv2(file_list_flu[i]) %>%
    mutate(Year = years[i])
  data_list_rsv[[i]] <- read_csv2(file_list_rsv[i]) %>%
    mutate(Year = years[i])
}

rm(file_list_flu, file_list_rsv, i)

# Format flu data:
data_list_flu <- lapply(data_list_flu, function(ix) {
  ix[-1, ] %>%
    filter(!is.na(Week)) %>%
    rename('n_samp' = 'No. of specimen tested',
           'n_h1' = 'Type A subtype H1',
           'n_h3' = 'Type A subtype H3',
           'n_h5' = 'Type A subtype H5',
           'n_h7' = 'Type A subtype H7',
           'n_h9' = 'Type A subtype H9',
           'n_b' = 'Type B',
           'n_c' = 'Type C') %>%
    select(Year, Week, n_samp, n_h1, n_h3, n_h5, n_h7, n_h9, n_b, n_c) %>%
    # select(18, 1, 3:4, 6, 8, 10, 12, 14, 16) %>%
    mutate(across(.cols = everything(), as.numeric))
})

# Format RSV data:
data_list_rsv <- lapply(data_list_rsv, function(ix) {
  ix[-1, ] %>%
    rename('Week' = 'Week number',
           'n_samp' = 'No. of specimen tested',
           'n_rsv' = 'RSV') %>%
    # select(13, 1, 3, 6) %>%
    select(Year, Week, n_samp, n_rsv) %>%
    mutate(across(.cols = everything(), as.numeric))
})
# Starting in 2020 week 8, ages are only <18!

# Combine:
dat_flu <- bind_rows(data_list_flu)
dat_rsv <- bind_rows(data_list_rsv)
expect_true(nrow(dat_flu) == nrow(dat_rsv))

dat_flu <- dat_flu %>%
  select(Year:n_h3, n_b)

dat_hk <- dat_flu %>%
  inner_join(dat_rsv, by = c('Year', 'Week'))
# expect_true(all(dat_hk$n_samp.x == dat_hk$n_samp.y))
if (!all(dat_hk$n_samp.x == dat_hk$n_samp.y)) {
  dat_hk %>% filter(n_samp.x != n_samp.y)
}
# 2018 week 30 could easily be a typo; starting 2020 week 8 number changes b/c RSV is <18 only I assume
# 2020 week 8 is also around when we start seeing highly reduced case counts due to lockdowns - can probably remove

dat_hk <- dat_hk %>%
  select(-n_samp.y) %>%
  rename('n_samp' = 'n_samp.x') %>%
  mutate(Time = 1:nrow(dat_hk)) %>%
  filter(Time < Time[Year == 2020 & Week == 8])

# Plot:
p1 <- ggplot(data = dat_hk %>% pivot_longer(n_h1:n_rsv, names_to = 'vir', values_to = 'val'),
             aes(x = Time, y = val, color = vir)) + geom_line() + theme_classic() + #facet_wrap(~ vir) +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years)
p2 <- ggplot(data = dat_hk %>% pivot_longer(n_h1:n_rsv, names_to = 'vir', values_to = 'val'),
             aes(x = Time, y = val / n_samp, color = vir)) + geom_line() + theme_classic() +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years)
grid.arrange(p1, p2, ncol = 1)

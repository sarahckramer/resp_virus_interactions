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
file_list_flu <- list.files(path = 'data/raw/hk/', pattern = 'flu', full.names = TRUE)
file_list_rsv <- list.files(path = 'data/raw/hk/', pattern = 'rsv', full.names = TRUE)
file_list_ili <- list.files(path = 'data/raw/hk/', pattern = 'ili', full.names = TRUE)

data_list_flu = data_list_rsv = data_list_ili = vector('list', length = length(file_list_flu))
for (i in 1:length(file_list_flu)) {
  data_list_flu[[i]] <- read_csv2(file_list_flu[i]) %>%
    mutate(Year = years[i])
  data_list_rsv[[i]] <- read_csv2(file_list_rsv[i]) %>%
    mutate(Year = years[i])
  data_list_ili[[i]] <- read.table(file_list_ili[i], sep = '\t', skip = 1) %>%
    mutate(Year = years[i])
}

rm(file_list_flu, file_list_rsv, file_list_ili, i)

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
  
  names(ix) <- make.names(names(ix))
  if (!('Rhinovirus...Enterovirus' %in% names(ix))) {
    ix <- ix %>%
      rename('Rhinovirus...Enterovirus' = 'X.Rhinovirus...Enterovirus')
  }
  
  if ('Human.metapneumovirus' %in% names(ix)) {
    ix <- ix %>%
      rename('n_hmpv' = 'Human.metapneumovirus')
  } else if ('Human.metapneumovirus.' %in% names(ix)) {
    ix <- ix %>%
      rename('n_hmpv' = 'Human.metapneumovirus.')
  }
  
  locate_diff <- ix[1, ] %>% unlist() %>% as.vector()
  locate_diff <- which(locate_diff == 'Differentiation performed')
  select_cols <- paste0('X', seq(locate_diff, locate_diff + 2))
  
  ix <- ix[-1, ] %>%
    rename('Week' = 'Week.number',
           'n_samp' = 'No..of.specimen.tested',
           'n_rsv' = 'RSV',
           'n_adeno' = 'Adenovirus',
           'n_rhino_entero' = 'Rhinovirus...Enterovirus')
  
  if (!('n_hmpv' %in% names(ix))) {
    ix <- ix %>%
      mutate(n_hmpv = NA)
  }
  
  ix <- ix %>%
    select(Year, Week, n_samp, n_rsv, n_adeno, n_hmpv, n_rhino_entero, all_of(select_cols)) %>%
    rename('n_diff' = select_cols[1],
           'n_rhino' = select_cols[2],
           'n_entero' = select_cols[3]) %>%
    mutate(across(.cols = everything(), as.numeric))
  
  ix
})
# Starting in 2020 week 8, ages are only <18!

# Format ILI data:
data_list_ili <- lapply(data_list_ili, function(ix) {
  names(ix)[1:2] <- c('GOPC', 'PMP.Clinics')
  ix %>%
    as_tibble() %>%
    mutate(Week = 1:nrow(ix)) %>%
    select(Year, Week, GOPC, PMP.Clinics) %>%
    mutate(across(.cols = everything(), as.numeric))
})
# General Out-patient Clinics (GOPC) and Private Medical Practitioner (PMP) Clinics

# Combine:
dat_flu <- bind_rows(data_list_flu)
dat_rsv <- bind_rows(data_list_rsv)
dat_ili <- bind_rows(data_list_ili)
expect_true(nrow(dat_flu) == nrow(dat_rsv))
expect_true(nrow(dat_flu) == nrow(dat_ili))

dat_flu <- dat_flu %>%
  select(Year:n_h3, n_b)

dat_hk <- dat_flu %>%
  inner_join(dat_rsv, by = c('Year', 'Week')) %>%
  inner_join(dat_ili, by = c('Year', 'Week'))
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

# Estimate number of rhino/entero cases based on proportion of differentiated samples that are rhino vs. entero:
dat_hk <- dat_hk %>%
  mutate(n_rhino_est = (n_rhino / n_diff) * n_rhino_entero,
         n_entero_est = (n_entero / n_diff) * n_rhino_entero,
         n_rhino_est_rnd = round(n_rhino_est),
         n_entero_est_rnd = round(n_entero_est)) %>%
  select(Time, Year:n_rhino_entero, n_rhino_est:n_entero_est_rnd, GOPC:PMP.Clinics)

# Save data:
write_csv(dat_hk, file = 'data/formatted/dat_hk.csv')

# Plot:
p1 <- ggplot(data = dat_hk %>% pivot_longer(c(n_h1:n_rsv, n_rhino_est), names_to = 'vir', values_to = 'val'),
             aes(x = Time, y = val, color = vir)) + geom_line() + theme_classic() + #facet_wrap(~ vir) +
  labs(x = 'Time', y = '# Samples Pos', color = 'Virus') +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years)
p2 <- ggplot(data = dat_hk %>% pivot_longer(c(n_h1:n_rsv, n_rhino_est), names_to = 'vir', values_to = 'val'),
             aes(x = Time, y = val / n_samp, color = vir)) + geom_line() + theme_classic() +
  labs(x = 'Time', y = '% Samples Pos', color = 'Virus') +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years) +
  theme(legend.position = 'bottom')
p3 <- ggplot(data = dat_hk %>% pivot_longer(GOPC:PMP.Clinics, names_to = 'clinic_type', values_to = 'val'),
       aes(x = Time, y = val, color = clinic_type)) + geom_line() + theme_classic() +
  labs(x = 'Time', y = 'ILI Rate (per 1000 consultations)', color = 'Clinic') +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years) +
  theme(legend.position = 'bottom')
grid.arrange(p2, p3, ncol = 1)

p4 <- ggplot(data = dat_hk %>% pivot_longer(c(n_h1:n_hmpv, n_rhino_est:n_entero_est), names_to = 'vir', values_to = 'val'),
       aes(x = Time, y = val / n_samp, color = vir)) + geom_line() + theme_classic() +
  labs(x = 'Time', y = '% Samples Pos', color = 'Virus') + facet_wrap(~ vir) +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years) +
  theme(legend.position = 'none')
print(p4)

p5 <- ggplot(data = dat_hk, aes(x = Time, y = n_samp)) + geom_line() +
  theme_classic() + labs(x = 'Time', y = '# of Vir. Samples') +
  scale_x_continuous(breaks = c(1, 53, 105, 158, 210, 262, 314, 366), labels = years)
print(p5)

# Output to pdf:
pdf('results/plots/data_Hong_Kong.pdf', width = 10, height = 6)
grid.arrange(p2, p3, ncol = 1)
print(p4)
print(p5)
dev.off()

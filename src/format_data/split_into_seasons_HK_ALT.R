# ---------------------------------------------------------------------------------------------------------------------
# Try dividing the data into regular winter seasons, as an alternative
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Set season start week:
start_week <- 46

# Read in data:
dat_hk <- read_csv('data/formatted/dat_hk.csv')

# Limit to data of interest:
dat_hk <- dat_hk %>%
  select(Time:n_rsv, n_rhino_est_rnd, GOPC:PMP.Clinics) %>%
  rename('n_rhino' = 'n_rhino_est_rnd')

# Label with northern-hemisphere seasons:
dat_hk <- dat_hk %>%
  mutate(season = NA) %>%
  mutate(season = ifelse((Year == 2014 & Week >= start_week) | (Year == 2015 & Week < start_week), 's14-15', season),
         season = ifelse((Year == 2015 & Week >= start_week) | (Year == 2016 & Week < start_week), 's15-16', season),
         season = ifelse((Year == 2016 & Week >= start_week) | (Year == 2017 & Week < start_week), 's16-17', season),
         season = ifelse((Year == 2017 & Week >= start_week) | (Year == 2018 & Week < start_week), 's17-18', season),
         season = ifelse((Year == 2018 & Week >= start_week) | (Year == 2019 & Week < start_week), 's18-19', season),
         season = ifelse((Year == 2019 & Week >= start_week) | (Year == 2020 & Week < start_week), 's19-20', season),
         season = ifelse(Year == 2014 & Week < start_week, 's13-14', season))
# In this case, do we remove 13-14 and 19-20, since they are not full seasons?

# Check for NAs:
dat_hk %>% pull(season) %>% is.na() %>% any() %>% expect_false()

# Get population sizes:
pop_dat <- read_csv('data/raw/hk/pop_dat_hk.csv', skip = 4) %>%
  rename('X1' = 'Year') %>%
  filter(X1 == 'Both sexes',
         X2 == 'All age groups') %>%
  select(!contains('_1')) %>%
  select(all_of(as.character(unique(dat_hk$Year))))

pop_dat <- pop_dat %>%
  pivot_longer(cols = everything(), names_to = 'season', values_to = 'pop') %>%
  mutate(season = as.numeric(season),
         pop = as.numeric(pop)) %>%
  mutate(season = paste0('s', str_sub(season - 1, 3, 4), '-', str_sub(season, 3, 4))) %>%
  mutate(pop = pop * 1000)

dat_hk <- dat_hk %>%
  left_join(pop_dat, by = 'season')

# Ensure that Time begins at 1 for each season:
dat_hk <- dat_hk %>%
  group_by(season) %>%
  mutate(Time = Time - min(Time) + 1) %>%
  ungroup()

# Break down by virus pairs and rename columns:
vir_list <- c('h1', 'h3', 'b', 'rsv', 'rhino')
dat_hk_pomp <- vector('list', length = 6)

counter <- 1
for (vir1_nm in vir_list[1:3]) {
  for (vir2_nm in vir_list[4:5]) {
    
    # Rename list element:
    names(dat_hk_pomp)[counter] <- paste(vir1_nm, vir2_nm, sep = '_')
    
    # Limit to data of interest:
    dat_out <- dat_hk %>%
      select(Time:n_samp, contains(vir1_nm), contains(vir2_nm), GOPC:pop)
    
    # Remove 19-20 season:
    dat_out <- dat_out %>%
      filter(season != 's19-20')
    
    # Format column names:
    dat_out <- dat_out %>%
      rename('time' = 'Time',
             'n_T' = 'n_samp') %>%
      rename_with(~ 'n_P1', .cols = contains(vir1_nm)) %>%
      rename_with(~ 'n_P2', .cols = contains(vir2_nm))
    
    # Remove seasons where a given threshold for flu is not exceeded:
    dat_out <- dat_out %>%
      group_by(season) %>%
      mutate(tot_pos = sum(n_P1)) %>%
      ungroup() %>%
      filter(tot_pos > 2000) %>%
      select(-tot_pos)
    
    # Store data in list:
    dat_hk_pomp[[counter]] <- dat_out
    
    # Iterate counter:
    counter <- counter + 1
    
  }
}

# Clean up:
rm(counter, vir1_nm, vir2_nm, dat_out)

# Add leading NAs to 2013-14 season:
dat_hk_pomp <- lapply(dat_hk_pomp, function(ix) {
  ix %>%
    complete(Week, season) %>%
    filter(!(Week == 53 & is.na(n_P1))) %>%
    select(time:GOPC, Week:season, pop) %>%
    mutate(Year = if_else(is.na(Year), 2013, Year),
           time = if_else(is.na(time), Week - 52, time),
           pop = if_else(is.na(pop), min(dat_hk$pop), pop)) %>%
    arrange(Year, Week) %>%
    group_by(season) %>%
    mutate(time = time - min(time) + 1) %>%
    ungroup()
})

# Plot outbreaks for each combination of pathogens:
plot_list <- vector('list', length = length(dat_hk_pomp))
for (i in 1:length(dat_hk_pomp)) {
  dat_temp <- dat_hk_pomp[[i]]
  
  p_temp <- ggplot(data = dat_temp %>%
                     mutate(prop_P1 = n_P1 / n_T, prop_P2 = n_P2 / n_T) %>%
                     pivot_longer(prop_P1:prop_P2, names_to = 'vir', values_to = 'val')
  ) +
    geom_line(aes(x = time, y = val, group = vir, col = vir)) +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    facet_wrap(~ season) +
    labs(x = 'Time (Weeks)', y = 'Proportion Positive', col = 'Virus', title = names(dat_hk_pomp)[i])
  plot_list[[i]] <- p_temp
  
}
print(plot_list)

# Clean up:
rm(i, dat_temp)

# Save plots to file:
pdf('results/plots/data_Hong_Kong_byOutbreak.pdf', width = 8, height = 5)
print(plot_list)
dev.off()

# Save formatted data:
write_rds(dat_hk_pomp, file = 'data/formatted/dat_hk_byOutbreak.rds')

# Clean up:
rm(list = ls())

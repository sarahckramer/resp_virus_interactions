# ---------------------------------------------------------------------------------------------------------------------
# Get good-quality start values for the amplitude of sinusoidal forcing from fit climate forcing parameters
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(viridis)
library(gridExtra)

# Read in MLEs:
mle_h1 <- read_rds('results/MLEs_flu_h1.rds')

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')

# Set duration of infection:
gamma1 <- 7/5
gamma2 <- 7/10

# Get results by season:
mle_h1 <- mle_h1 %>%
  select(eta_temp1:eta_ah2,
         contains('R10'),
         contains('R20'),
         contains('R120'),
         contains('Ri')) %>%
  pivot_longer(-c(eta_temp1:eta_ah2),
               names_to = 'parameter',
               values_to = 'mle') %>%
  mutate(season = str_sub(parameter, 2, 6),
         parameter = str_sub(parameter, 8)) %>%
  pivot_wider(names_from = parameter,
              values_from = mle) %>%
  select(eta_temp1:eta_ah2, Ri1:Ri2, R10:R120, season)

# Get vector of seasons:
seasons <- unique(mle_h1$season)

# Calculate beta over time based on fit parameter values:
beta_t <- vector('list', length = length(seasons))

for (yr_index in 1:length(seasons)) {
  
  seas <- seasons[yr_index]
  
  dat_temp <- hk_dat[['h1_rsv']] %>%
    filter(season == paste0('s', seas)) %>%
    inner_join(dat_clim,
               by = c('Year' = 'year',
                      'Week' = 'week')) %>%
    select(time, temp, ah)
  
  mle_temp <- mle_h1 %>%
    filter(season == seas)
  
  # beta1_temp <- unlist(mle_temp['Ri1']) / (1.0 - (unlist(mle_temp['R10']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp) * gamma1
  # beta2_temp <- unlist(mle_temp['Ri2']) / (1.0 - (unlist(mle_temp['R20']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp) * gamma2
  
  beta_temp_list <- vector('list', length = nrow(mle_temp))
  
  for (i in 1:length(beta_temp_list)) {
    
    beta1_temp <- unlist(mle_temp[i, ]['Ri1']) / (1.0 - (unlist(mle_temp[i, ]['R10']) + unlist(mle_temp[i, ]['R120']))) * exp(unlist(mle_temp[i, ]['eta_ah1']) * dat_temp$ah + unlist(mle_temp[i, ]['eta_temp1']) * dat_temp$temp) * gamma1
    beta2_temp <- unlist(mle_temp[i, ]['Ri2']) / (1.0 - (unlist(mle_temp[i, ]['R20']) + unlist(mle_temp[i, ]['R120']))) * exp(unlist(mle_temp[i, ]['eta_ah2']) * dat_temp$ah + unlist(mle_temp[i, ]['eta_temp2']) * dat_temp$temp) * gamma2
    
    beta1_temp <- bind_cols(1:max(dat_temp$time), beta1_temp, seas, i)
    beta2_temp <- bind_cols(1:max(dat_temp$time), beta2_temp, seas, i)
    
    names(beta1_temp) <- c('time', 'beta1', 'season', 'result')
    names(beta2_temp) <- c('time', 'beta2', 'season', 'result')
    
    beta_temp <- beta1_temp %>%
      inner_join(beta2_temp,
                 by = c('time', 'season', 'result')) %>%
      select(time, beta1, beta2, season, result)
    
    beta_temp_list[[i]] <- beta_temp
    
  }
  
  beta_t[[yr_index]] <- bind_rows(beta_temp_list)
  
}
rm(yr_index, seas, dat_temp, mle_temp, beta_temp_list, beta1_temp, beta2_temp, beta_temp, i)

# Combine results into tibble:
beta_t <- bind_rows(beta_t)

# Get minimum and maximum values for plotting:
min_val <- min(c(min(beta_t$beta1), min(beta_t$beta2)))
max_val <- max(c(max(beta_t$beta1), max(beta_t$beta2)))

# Plot climate forcing over time for all seasons:
p_a <- ggplot(data = beta_t, aes(x = time, y = beta1, col = season, group = paste(season, result))) +
  geom_line() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.97)) +
  scale_x_continuous(breaks = seq(1, 52, by = 5),
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  # scale_y_continuous(limits = c(min_val, max_val)) +
  scale_y_log10(limits = c(min_val, max_val)) +
  scale_color_manual(values = viridis(6)) +
  labs(x = 'Week Number', y = expression(beta[1]))
p_b <- ggplot(data = beta_t, aes(x = time, y = beta2, col = season, group = paste(season, result))) +
  geom_line() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.97)) +
  scale_x_continuous(breaks = seq(1, 52, by = 5),
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_y_log10(limits = c(min_val, max_val)) +
  scale_color_manual(values = viridis(6)) +
  labs(x = 'Week Number', y = expression(beta[2]))

p_legend <- ggplot(data = beta_t, aes(x = time, y = beta1, col = season)) +
  geom_line() +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = 'Season')
p_legend <- ggplotGrob(p_legend)$grobs[[which(sapply(ggplotGrob(p_legend)$grobs, function(x) x$name) == 'guide-box')]]

p_beta <- arrangeGrob(arrangeGrob(p_a, p_b, ncol = 1), p_legend, nrow = 2, heights = c(15, 1))
plot(p_beta)

# Fit sine wave to each seasonal beta and get range of amplitudes:

# Clean up:
rm(list = ls())

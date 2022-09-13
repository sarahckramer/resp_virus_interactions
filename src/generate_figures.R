# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript (main text)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(lemon)
library(gridExtra)
library(patchwork)
library(RColorBrewer)

# ---------------------------------------------------------------------------------------------------------------------

# Figure 1: Visualize Hong Kong data
dat_hk <- read_csv('data/formatted/dat_hk.csv')

dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_h1, n_b, n_rsv, GOPC) %>%
  mutate(n_h1 = n_h1 / n_samp * 100,
         n_b = n_b / n_samp * 100,
         n_rsv = n_rsv / n_samp * 100)

dat_pos <- dat_hk %>%
  select(-c(n_samp, GOPC)) %>%
  pivot_longer(n_h1:n_rsv,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_h1', 'n_b', 'n_rsv'))) %>%
  mutate(virus = recode(virus, n_h1 = 'Influenza (H1)', n_b = 'Influenza (B)', n_rsv = 'RSV'))

x_lab_breaks <- dat_hk %>% filter(Week == 1) %>% pull(Time)

p1a <- ggplot(data = dat_pos, aes(x = Time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Year', y = '\n% Positive', col = 'Virus', tag = 'A')
p1a <- reposition_legend(p1a, position = 'top left', plot = FALSE)

p1b <- ggplot(data = dat_hk, aes(x = Time, y = n_samp)) + geom_line() +
  theme_classic() + theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          plot.tag = element_text(size = 22),
                          plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  labs(x = 'Year', y = 'Total # of Tests', tag = 'B')

p1c <- ggplot(data = dat_hk, aes(x = Time, y = GOPC)) + geom_line() +
  theme_classic() + theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          plot.tag = element_text(size = 22),
                          plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  labs(x = 'Year', y = 'ILI Cases per 1000\nConsultations', tag = 'C')

fig1 <- arrangeGrob(p1a, p1b, p1c, ncol = 1)
plot(fig1)

ggsave('results/plots/figures_for_manuscript/Figure1.svg', width = 9.5, height = 8.5, fig1)
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Figure 2: Model schematic
# Not generated in R

# ---------------------------------------------------------------------------------------------------------------------

# Figure 3: Observed vs. simulated cases
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

dat_h1 <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_rsv
dat_b <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$b_rsv

mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
mle_b <- read_rds('results/MLEs_flu_b.rds')

vir1 <- 'flu_h1'
prof_lik <- FALSE; lag_val <- 0
source('src/functions/setup_global_likelilhood.R')

sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Setup
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- mle_h1[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  # Run deterministic simulation and calculate mean/median/IQR of binomial distribution:
  traj_temp <- trajectory(resp_mod, format = 'data.frame') %>%
    mutate(season = yr) %>%
    left_join(dat_h1, by = c('time', 'season')) %>%
    rename('i_ILI' = 'GOPC') %>%
    mutate(i_ILI = i_ILI / 1000)
  
  rho1 <- as.numeric(pars_temp['rho1'])
  rho2 <- as.numeric(pars_temp['rho2'])
  alpha <- as.numeric(pars_temp['alpha'])
  phi <- as.numeric(pars_temp['phi'])
  
  rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H1 / traj_temp$i_ILI
  rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H2 / traj_temp$i_ILI
  
  rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
  rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
  
  expect_equal(nrow(traj_temp), length(rho1_w))
  expect_equal(nrow(traj_temp), length(rho2_w))
  
  traj_temp$rho1_w <- rho1_w
  traj_temp$rho2_w <- rho2_w
  
  traj_temp <- traj_temp %>%
    as_tibble() %>%
    mutate(mean1 = rho1_w * n_T,
           mean2 = rho2_w * n_T,
           # median1 = qbinom(p = 0.50, size = n_T, prob = rho1_w),
           # median2 = qbinom(p = 0.50, size = n_T, prob = rho2_w),
           lower1 = qbinom(p = 0.025, size = n_T, prob = rho1_w),
           lower2 = qbinom(p = 0.025, size = n_T, prob = rho2_w),
           upper1 = qbinom(p = 0.975, size = n_T, prob = rho1_w),
           upper2 = qbinom(p = 0.975, size = n_T, prob = rho2_w)) %>%
    select(time, season, Year, Week, mean1:upper2)
  
  # Run several stochastic simulations and calculate mean/median/IQR:
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    as_tibble() %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              # median1_stoch = median(n_P1),
              # median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.025, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.025, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.975, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.975, na.rm = TRUE))
  
  # Combine and store:
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  sim_temp <- traj_temp %>%
    inner_join(sim_temp, by = 'time')
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  
  sim_list[[i]] <- sim_temp
  
}

res_h1 <- bind_rows(sim_list) %>%
  inner_join(dat_h1,
             by = c('time', 'season', 'Year', 'Week')) %>%
  select(-c(n_T, GOPC, pop)) %>%
  mutate(vir1 = 'H1N1')

vir1 <- 'flu_b'
prof_lik <- FALSE; lag_val <- 0
source('src/functions/setup_global_likelilhood.R')

sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Setup
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- mle_b[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  # Run deterministic simulation and calculate mean/median/IQR of binomial distribution:
  traj_temp <- trajectory(resp_mod, format = 'data.frame') %>%
    mutate(season = yr) %>%
    left_join(dat_b, by = c('time', 'season')) %>%
    rename('i_ILI' = 'GOPC') %>%
    mutate(i_ILI = i_ILI / 1000)
  
  rho1 <- as.numeric(pars_temp['rho1'])
  rho2 <- as.numeric(pars_temp['rho2'])
  alpha <- as.numeric(pars_temp['alpha'])
  phi <- as.numeric(pars_temp['phi'])
  
  rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H1 / traj_temp$i_ILI
  rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H2 / traj_temp$i_ILI
  
  rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
  rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
  
  expect_equal(nrow(traj_temp), length(rho1_w))
  expect_equal(nrow(traj_temp), length(rho2_w))
  
  traj_temp$rho1_w <- rho1_w
  traj_temp$rho2_w <- rho2_w
  
  traj_temp <- traj_temp %>%
    as_tibble() %>%
    mutate(mean1 = rho1_w * n_T,
           mean2 = rho2_w * n_T,
           # median1 = qbinom(p = 0.50, size = n_T, prob = rho1_w),
           # median2 = qbinom(p = 0.50, size = n_T, prob = rho2_w),
           lower1 = qbinom(p = 0.025, size = n_T, prob = rho1_w),
           lower2 = qbinom(p = 0.025, size = n_T, prob = rho2_w),
           upper1 = qbinom(p = 0.975, size = n_T, prob = rho1_w),
           upper2 = qbinom(p = 0.975, size = n_T, prob = rho2_w)) %>%
    select(time, season, Year, Week, mean1:upper2)
  
  # Run several stochastic simulations and calculate mean/median/IQR:
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    as_tibble() %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              # median1_stoch = median(n_P1),
              # median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.025, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.025, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.975, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.975, na.rm = TRUE))
  
  # Combine and store:
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  sim_temp <- traj_temp %>%
    inner_join(sim_temp, by = 'time')
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  
  sim_list[[i]] <- sim_temp
  
}

res_b <- bind_rows(sim_list) %>%
  inner_join(dat_b,
             by = c('time', 'season', 'Year', 'Week')) %>%
  select(-c(n_T, GOPC, pop)) %>%
  mutate(vir1 = 'B')

# Plot means and 95% CIs from binomial distribution:
res <- bind_rows(res_h1, res_b)
res <- res %>%
  select(time:Week, vir1, mean1:mean2, lower1:upper2, n_P1:n_P2) %>%
  pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'obs') %>%
  pivot_longer(mean1:mean2, names_to = 'virus1', values_to = 'mean') %>%
  pivot_longer(lower1:lower2, names_to = 'virus2', values_to = 'lower') %>%
  pivot_longer(upper1:upper2, names_to = 'virus3', values_to = 'upper') %>%
  mutate(virus = str_sub(virus, 4),
         virus1 = str_sub(virus1, 5),
         virus2 = str_sub(virus2, 6),
         virus3 = str_sub(virus3, 6)) %>%
  filter(virus == virus1,
         virus == virus2,
         virus == virus3) %>%
  select(-c(virus1, virus2, virus3)) %>%
  mutate(virus = if_else(virus == 1, 'Influenza', 'RSV'))

p_legend <- ggplot(data = res, aes(x = mean, y = obs, color = season)) +
  geom_point() +
  facet_grid(vir1 ~ virus) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = 'Season')
p_legend <- ggplotGrob(p_legend)$grobs[[which(sapply(ggplotGrob(p_legend)$grobs, function(x) x$name) == 'guide-box')]]

p3a <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  # geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
  #             inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(title = 'Influenza (A(H1N1)-RSV)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(title = 'RSV (A(H1N1)-RSV)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(title = 'Influenza (B-RSV)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(title = 'RSV (B-RSV)', x = 'Simulated Cases', y = 'Observed Cases',color = 'Season', tag = 'D')

# p_column_title1 <- tableGrob('Influenza', theme = ttheme_minimal(core = list(fg_params = list(fontsize = 16))))
# p_column_title2 <- tableGrob('RSV', theme = ttheme_minimal(core = list(fg_params = list(fontsize = 16))))

# fig3 <- arrangeGrob(arrangeGrob(p_column_title1, p_column_title2, ncol = 2), arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 3, heights = c(1, 15, 1))
fig3 <- arrangeGrob(arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3.svg', fig3, width = 14, height = 8)

# For each virus, virus-virus pair, and season, calculate coefficient of efficiency:
# https://stats.stackexchange.com/questions/185898/difference-between-nash-sutcliffe-efficiency-and-coefficient-of-determination
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
r2a_list = r2b_list = r2c_list = r2d_list = c()
for (yr in seasons) {
  
  res_temp <- res %>%
    filter(season == yr)
  
  if (yr != 's14-15') {
    
    r2_a <- res_temp %>%
      filter(vir1 == 'H1N1' & virus == 'Influenza') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total) %>%
      pull(r2)
    
    r2_b <- res_temp %>%
      filter(vir1 == 'H1N1' & virus == 'RSV') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total) %>%
      pull(r2)
    
    r2a_list <- c(r2a_list, r2_a)
    r2b_list <- c(r2b_list, r2_b)
    
  } else {
    
    r2a_list <- c(r2a_list, NA)
    r2b_list <- c(r2b_list, NA)
    
  }
  
  
  if (yr != 's16-17') {
    
    r2_c <- res_temp %>%
      filter(vir1 == 'B' & virus == 'Influenza') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total) %>%
      pull(r2)
    
    r2_d <- res_temp %>%
      filter(vir1 == 'B' & virus == 'RSV') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total) %>%
      pull(r2)
    
    r2c_list <- c(r2c_list, r2_c)
    r2d_list <- c(r2d_list, r2_d)
    
  } else {
    
    r2c_list <- c(r2c_list, NA)
    r2d_list <- c(r2d_list, NA)
    
  }
  
}

print(summary(r2a_list))
print(summary(r2b_list))
print(summary(r2c_list))
print(summary(r2d_list))

# For each virus, virus-virus pair, and season, calculate the proportion of observations falling within confidence intervals:
prop_a_list = prop_b_list = prop_c_list = prop_d_list = c()
for (yr in seasons) {
  
  prop_a <- res %>%
    filter(season == yr,
           vir1 == 'H1N1',
           virus == 'Influenza') %>%
    mutate(within_ci = obs >= lower & obs <= upper) %>%
    summarise(prop = sum(within_ci, na.rm = TRUE) / length(within_ci[!is.na(within_ci)])) %>%
    pull(prop)
  
  prop_b <- res %>%
    filter(season == yr,
           vir1 == 'H1N1',
           virus == 'RSV') %>%
    mutate(within_ci = obs >= lower & obs <= upper) %>%
    summarise(prop = sum(within_ci, na.rm = TRUE) / length(within_ci[!is.na(within_ci)])) %>%
    pull(prop)
  
  prop_c <- res %>%
    filter(season == yr,
           vir1 == 'B',
           virus == 'Influenza') %>%
    mutate(within_ci = obs >= lower & obs <= upper) %>%
    summarise(prop = sum(within_ci, na.rm = TRUE) / length(within_ci[!is.na(within_ci)])) %>%
    pull(prop)
  
  prop_d <- res %>%
    filter(season == yr,
           vir1 == 'B',
           virus == 'RSV') %>%
    mutate(within_ci = obs >= lower & obs <= upper) %>%
    summarise(prop = sum(within_ci, na.rm = TRUE) / length(within_ci[!is.na(within_ci)])) %>%
    pull(prop)
  
  prop_a_list <- c(prop_a_list, prop_a)
  prop_b_list <- c(prop_b_list, prop_b)
  prop_c_list <- c(prop_c_list, prop_c)
  prop_d_list <- c(prop_d_list, prop_d)
  
}

print(summary(prop_a_list))
print(summary(prop_b_list))
print(summary(prop_c_list))
print(summary(prop_d_list))

rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Figure 4: Plot heatmaps from simulation study
file_list_hk <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'TEMPERATE', full.names = TRUE)

file_list_hk_red <- list.files('results/vaccine_simulation_study/simulations/thetalambdavacc0.50/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp_red <- list.files('results/vaccine_simulation_study/simulations/thetalambdavacc0.50/', pattern = 'TEMPERATE', full.names = TRUE)

res_list <- vector('list', length(file_list_hk))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk[i]) %>% mutate(season = str_split(file_list_hk[i], '_')[[1]][5])
}
res_hk <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop',
         scenario = 'natural')

res_list <- vector('list', length(file_list_temp))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_temp[i]) %>% mutate(season = str_split(file_list_temp[i], '_')[[1]][5])
}
res_temp <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'temp',
         scenario = 'natural')

res_list <- vector('list', length(file_list_hk_red))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk_red[i]) %>% mutate(season = str_split(file_list_hk_red[i], '_')[[1]][5])
}
res_hk_red <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop',
         scenario = 'half')

res_list <- vector('list', length(file_list_temp_red))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_temp_red[i]) %>% mutate(season = str_split(file_list_temp_red[i], '_')[[1]][5])
}
res_temp_red <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'temp',
         scenario = 'half')

rm(i, file_list_hk, file_list_temp, file_list_hk_red, file_list_temp_red, res_list)

res <- res_hk %>%
  bind_rows(res_temp) %>%
  bind_rows(res_hk_red) %>%
  bind_rows(res_temp_red)
rm(res_hk, res_temp, res_hk_red, res_temp_red)

res_metrics <- res %>%
  group_by(climate, scenario, season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(climate, scenario, season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

res_metrics <- res_metrics %>%
  filter(vacc_cov <= 0.60) %>%
  mutate(vacc_cov = vacc_cov * 100)

# res_metrics_AVG <- res_metrics %>%
#   group_by(climate, scenario, vacc_cov, vacc_time) %>%
#   summarise(ar2_impact = median(ar2_impact))

res <- res %>%
  filter((climate == 'subtrop' & season == 's18-19') |
           (climate == 'temp' & season == 's17-18'))
res_metrics <- res_metrics %>%
  filter((climate == 'subtrop' & season == 's18-19') |
           (climate == 'temp' & season == 's17-18'))

upper_bound_ar <- max(res_metrics$ar2_impact)

p_legend1 <- ggplot(data= res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
                    aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.5, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  labs(x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')
p_legend1 <- ggplotGrob(p_legend1)$grobs[[which(sapply(ggplotGrob(p_legend1)$grobs, function(x) x$name) == 'guide-box')]]

# res_metrics %>%
#   filter(climate == 'temp' & scenario == 'natural') %>%
#   filter(ar2_impact == min(ar2_impact))

res_simA <- res %>%
  filter(climate == 'temp',
         scenario == 'natural',
         vacc_cov == '0.6',
         vacc_time == '0') %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV'))

p_legend2 <- ggplot(data = res_simA, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  scale_linetype(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')
p_legend2 <- ggplotGrob(p_legend2)$grobs[[which(sapply(ggplotGrob(p_legend2)$grobs, function(x) x$name) == 'guide-box')]]

p4a <- ggplot(data= res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 0, y = 60, shape = 16, size = 3, color = 'black', fill = 'black') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.2, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], '=', theta[lambda*1], ')')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')

p4a_sim <- ggplot(data = res_simA, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(x = 52, y = 1.13, shape = 16, size = 3, col = 'black') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none') +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4b <- ggplot(data= res_metrics %>% filter(climate == 'subtrop' & scenario == 'natural'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 15, y = 60, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 60, shape = 17, size = 3, color = 'black', fill = 'black') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], '=', theta[lambda*1], ')')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'B')

# res_metrics %>%
#   filter(climate == 'subtrop' & scenario == 'natural') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simB <- res %>%
  filter(climate == 'subtrop',
         scenario == 'natural',
         vacc_cov == 0.6,
         vacc_time %in% c('0', '15')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_time = factor(vacc_time, levels = c('15', '0')))

p4b_sim <- ggplot(data = res_simB, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = as.numeric(as.character(vacc_time))), lty = 2) +
  geom_point(x = 52, y = 3.8, aes(shape = vacc_time), size = 3, col = 'black', fill = 'black') +
  facet_wrap(~ vacc_time, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4c <- ggplot(data= res_metrics %>% filter(climate == 'temp' & scenario == 'half'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 0, y = 60, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 20, shape = 17, size = 3, color = 'black', fill = 'black') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], '= 0.5)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'C')

# res_metrics %>%
#   filter(climate == 'temp' & scenario == 'half') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simC <- res %>%
  filter(climate == 'temp',
         scenario == 'half',
         vacc_cov %in% c('0.2', '0.6'),
         vacc_time == '0') %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_cov = factor(vacc_cov, levels = c('0.6', '0.2')))

p4c_sim <- ggplot(data = res_simC, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = vacc_time), lty = 2) +
  geom_point(x = 52, y = 1.13, aes(shape = vacc_cov), size = 3, col = 'black') +
  facet_wrap(~ vacc_cov, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4d <- ggplot(data= res_metrics %>% filter(climate == 'subtrop' & scenario == 'half'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 19, y = 60, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 40, shape = 17, size = 3, color = 'black', fill = 'black') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], '= 0.5)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'D')

# res_metrics %>%
#   filter(climate == 'subtrop' & scenario == 'half') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simD <- res %>%
  filter(climate == 'subtrop',
         scenario == 'half') %>%
  filter((vacc_cov == '0.4' & vacc_time == '0') | (vacc_cov == '0.6' & vacc_time == '19')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_cov = factor(vacc_cov, levels = c('0.6', '0.4')))

p4d_sim <- ggplot(data = res_simD, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = vacc_time), lty = 2) +
  geom_point(x = 52, y = 3.8, aes(shape = vacc_cov), size = 3, col = 'black') +
  facet_wrap(~ vacc_cov, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

fig4 <- arrangeGrob(arrangeGrob(arrangeGrob(p4a, p4a_sim, layout_matrix = rbind(c(1, 2), c(1, NA)), heights = c(5.86, 4.14), widths = c(5.5, 4.5)),
                                arrangeGrob(p4b, p4b_sim, nrow = 1, widths = c(5.5, 4.5)), nrow = 1),
                    arrangeGrob(arrangeGrob(p4c, p4c_sim, nrow = 1, widths = c(5.5, 4.5)),
                                arrangeGrob(p4d, p4d_sim, nrow = 1, widths = c(5.5, 4.5)), nrow = 1),
                    arrangeGrob(p_legend1, p_legend2, layout_matrix = rbind(c(NA, 1, 2, NA)), widths = c(3, 2, 2, 3)),
                    nrow = 3, heights = c(12, 12, 2.25))
plot(fig4)

ggsave('results/plots/figures_for_manuscript/Figure4.svg', fig4, width = 20, height = 10)
rm(list = ls())

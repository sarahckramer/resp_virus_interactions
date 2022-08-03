# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript and supplement
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(lemon)
library(gridExtra)
library(patchwork)

# ---------------------------------------------------------------------------------------------------------------------

# Figure 1: Model schematic
# Not generated in R

# ---------------------------------------------------------------------------------------------------------------------

# Figure 2: Visualize Hong Kong data
dat_hk <- read_csv('data/formatted/dat_hk.csv')

dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_h1, n_b, n_rsv, GOPC) %>%
  mutate(n_h1 = n_h1 / n_samp * 100,
         n_b = n_b / n_samp * 100,
         n_rsv = n_rsv / n_samp * 100) %>%
  select(-n_samp)

dat_pos <- dat_hk %>%
  select(-GOPC) %>%
  pivot_longer(n_h1:n_rsv,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_h1', 'n_b', 'n_rsv'))) %>%
  mutate(virus = recode(virus, n_h1 = 'Influenza (H1)', n_b = 'Influenza (B)', n_rsv = 'RSV'))

x_lab_breaks <- dat_hk %>% filter(Week == 1) %>% pull(Time)

p2a <- ggplot(data = dat_pos, aes(x = Time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.05, 0.96)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Year', y = '\n% Positive', col = 'Virus', tag = 'A')
p2a <- reposition_legend(p2a, position = 'top left', plot = FALSE)
p2b <- ggplot(data = dat_hk, aes(x = Time, y = GOPC)) + geom_line() +
  theme_classic() + theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          plot.tag = element_text(size = 22),
                          plot.tag.position = c(0.05, 0.96)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  labs(x = 'Year', y = 'ILI Rate per 1000\nConsultations', tag = 'B')

# fig2 = p2a + p2b + plot_layout(ncol = 1)
# print(fig2)

fig2 <- arrangeGrob(p2a, p2b, ncol = 1)
plot(fig2)

ggsave('results/plots/figures_for_manuscript/Figure2.svg', width = 9.5, height = 6, fig2)
rm(dat_hk, dat_pos, p2a, p2b, x_lab_breaks)

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
           median1 = qbinom(p = 0.50, size = n_T, prob = rho1_w),
           median2 = qbinom(p = 0.50, size = n_T, prob = rho2_w),
           lower1 = qbinom(p = 0.25, size = n_T, prob = rho1_w),
           lower2 = qbinom(p = 0.25, size = n_T, prob = rho2_w),
           upper1 = qbinom(p = 0.75, size = n_T, prob = rho1_w),
           upper2 = qbinom(p = 0.75, size = n_T, prob = rho2_w)) %>%
    select(time, season, Year, Week, mean1:upper2)
  
  # Run several stochastic simulations and calculate mean/median/IQR:
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    as_tibble() %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              median1_stoch = median(n_P1),
              median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.25, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.25, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.75, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.75, na.rm = TRUE))
  
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
           median1 = qbinom(p = 0.50, size = n_T, prob = rho1_w),
           median2 = qbinom(p = 0.50, size = n_T, prob = rho2_w),
           lower1 = qbinom(p = 0.25, size = n_T, prob = rho1_w),
           lower2 = qbinom(p = 0.25, size = n_T, prob = rho2_w),
           upper1 = qbinom(p = 0.75, size = n_T, prob = rho1_w),
           upper2 = qbinom(p = 0.75, size = n_T, prob = rho2_w)) %>%
    select(time, season, Year, Week, mean1:upper2)
  
  # Run several stochastic simulations and calculate mean/median/IQR:
  sim_temp <- simulate(resp_mod, nsim = 100, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    as_tibble() %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              median1_stoch = median(n_P1),
              median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.25, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.25, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.75, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.75, na.rm = TRUE))
  
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

# Plot means from binomial distribution (with or without error bars):
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
       aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'RSV'),
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500), limits = c(9, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'Influenza'),
            aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'RSV'),
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'D')

fig3 <- arrangeGrob(arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3_DETERM.svg', fig3, width = 14, height = 8)

p3a <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500), limits = c(9, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'D')

fig3 <- arrangeGrob(arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3_DETERM_errorbars.svg', fig3, width = 14, height = 8)

# For each virus and virus-virus pair, calculate correlation coefficient:
cor.test(unlist(res[res$vir1 == 'H1N1' & res$virus == 'Influenza', 'mean']), unlist(res[res$vir1 == 'H1N1' & res$virus == 'Influenza', 'obs']))
cor.test(unlist(res[res$vir1 == 'H1N1' & res$virus == 'RSV', 'mean']), unlist(res[res$vir1 == 'H1N1' & res$virus == 'RSV', 'obs']))
cor.test(unlist(res[res$vir1 == 'B' & res$virus == 'Influenza', 'mean']), unlist(res[res$vir1 == 'B' & res$virus == 'Influenza', 'obs']))
cor.test(unlist(res[res$vir1 == 'B' & res$virus == 'RSV', 'mean']), unlist(res[res$vir1 == 'B' & res$virus == 'RSV', 'obs']))

# Plot means from fully stochastic simulations (with or without error bars):
res <- bind_rows(res_h1, res_b)
res <- res %>%
  select(time:Week, vir1, mean1_stoch:mean2_stoch, lower1_stoch:upper2_stoch, n_P1:n_P2) %>%
  pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'obs') %>%
  pivot_longer(mean1_stoch:mean2_stoch, names_to = 'virus1', values_to = 'mean') %>%
  pivot_longer(lower1_stoch:lower2_stoch, names_to = 'virus2', values_to = 'lower') %>%
  pivot_longer(upper1_stoch:upper2_stoch, names_to = 'virus3', values_to = 'upper') %>%
  mutate(virus = str_sub(virus, 4),
         virus1 = str_sub(virus1, 5, 5),
         virus2 = str_sub(virus2, 6, 6),
         virus3 = str_sub(virus3, 6, 6)) %>%
  filter(virus == virus1,
         virus == virus2,
         virus == virus3) %>%
  select(-c(virus1, virus2, virus3)) %>%
  mutate(virus = if_else(virus == 1, 'Influenza', 'RSV'))

# Plot means from binomial distribution (with or without error bars):
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
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'RSV'),
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500), limits = c(9, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'Influenza'),
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'RSV'),
              aes(x = mean, y = obs, color = season)) +
  geom_point() +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'D')

fig3 <- arrangeGrob(arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3_STOCH.svg', fig3, width = 14, height = 8)

p3a <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(vir1 == 'H1N1', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500), limits = c(9, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(vir1 == 'B', virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_pointrange(size = 0.2) +
  geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
              inherit.aes = FALSE, color = 'black', se = FALSE) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Simulated Cases', y = 'Observed Cases', color = 'Season', tag = 'D')

fig3 <- arrangeGrob(arrangeGrob(p3a, p3b, p3c, p3d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3_STOCH_errorbars.svg', fig3, width = 14, height = 8)

# For each virus and virus-virus pair, calculate correlation coefficient:
cor.test(unlist(res[res$vir1 == 'H1N1' & res$virus == 'Influenza', 'mean']), unlist(res[res$vir1 == 'H1N1' & res$virus == 'Influenza', 'obs']))
cor.test(unlist(res[res$vir1 == 'H1N1' & res$virus == 'RSV', 'mean']), unlist(res[res$vir1 == 'H1N1' & res$virus == 'RSV', 'obs']))
cor.test(unlist(res[res$vir1 == 'B' & res$virus == 'Influenza', 'mean']), unlist(res[res$vir1 == 'B' & res$virus == 'Influenza', 'obs']))
cor.test(unlist(res[res$vir1 == 'B' & res$virus == 'RSV', 'mean']), unlist(res[res$vir1 == 'B' & res$virus == 'RSV', 'obs']))

# ---------------------------------------------------------------------------------------------------------------------

# Figure 4: Plot heatmaps from simulation study
file_list_hk <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'TEMPERATE', full.names = TRUE)

res_list <- vector('list', length(file_list_hk))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk[i]) %>% mutate(season = str_sub(file_list_hk[i], 62, 67))
}
res_hk <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop')

res_list <- vector('list', length(file_list_temp))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_temp[i]) %>% mutate(season = str_sub(file_list_temp[i], 62, 67))
}

res_temp <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'temp')
rm(i, file_list_hk, file_list_temp, res_list)

res <- res_hk %>%
  bind_rows(res_temp)
rm(res_hk, res_temp)

res_metrics <- res %>%
  group_by(climate, season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(climate, season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

res_metrics <- res_metrics %>%
  filter(vacc_cov <= 0.60)

res_metrics_AVG <- res_metrics %>%
  group_by(climate, vacc_cov, vacc_time) %>%
  summarise(ar2_impact = median(ar2_impact))

upper_bound_ar <- max(res_metrics_AVG$ar2_impact)

# p4 <- ggplot(data= res_metrics_AVG,
#              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
#   geom_tile() + facet_wrap(~ climate, ncol = 2) +
#   theme_classic() + theme(strip.text = element_blank()) +
#   scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1),
#                        breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
#   scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
#   labs(x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'Impact')
# p4

p4a <- ggplot(data= res_metrics_AVG %>% filter(climate == 'temp'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.2, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'Impact', tag = 'A')
p4b <- ggplot(data= res_metrics_AVG %>% filter(climate == 'subtrop'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'Impact', tag = 'B')

fig4 <- grid_arrange_shared_legend(p4a, p4b, ncol = 2, plot = FALSE)
ggsave('results/plots/figures_for_manuscript/Figure4.svg', fig4, width = 10, height = 5)

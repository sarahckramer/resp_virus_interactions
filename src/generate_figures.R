# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript (main text)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(lemon)
library(gridExtra)
library(RColorBrewer)

# ---------------------------------------------------------------------------------------------------------------------

# Figure 1: Visualize data
dat_hk <- read_csv('data/formatted/dat_hk.csv')

dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_samp, n_h1_b, n_rsv, GOPC) %>%
  mutate(n_h1_b = n_h1_b / n_samp * 100,
         n_rsv = n_rsv / n_samp * 100)

dat_pos_hk <- dat_hk %>%
  select(-c(n_samp, GOPC)) %>%
  pivot_longer(n_h1_b:n_rsv,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_h1_b', 'n_rsv'))) %>%
  mutate(virus = recode(virus, n_h1_b = 'Influenza A(H1N1) + B    ', n_rsv = 'RSV'))

dat_can <- read_csv('data/formatted/dat_canada.csv')

dat_can <- dat_can %>%
  mutate(n_P1 = n_P1 / n_T1 * 100,
         n_P2 = n_P2 / n_T2 * 100,
         i_ILI = i_ILI * 1000) %>%
  select(time, year:week, n_P1:n_P2, i_ILI) %>%
  mutate(time = 1:nrow(dat_can))

dat_pos_can <- dat_can %>%
  select(-c(i_ILI, tot_a, tot_b)) %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_P1', 'n_P2'))) %>%
  mutate(virus = recode(virus, n_P1 = 'Influenza    ', n_P2 = 'RSV'))

x_lab_breaks_hk <- dat_hk %>% filter(Week == 1) %>% pull(Time)
x_lab_breaks_can <- dat_can %>% filter(week == 1) %>% pull(time)

p1a <- ggplot(data = dat_pos_hk, aes(x = Time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks_hk, labels = 2014:2019) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Year', y = '\n% Positive', col = 'Virus', tag = 'A')
p1a <- reposition_legend(p1a, position = 'top left', plot = FALSE)

p1b <- ggplot(data = dat_pos_can, aes(x = time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks_can, labels = 2011:2014) +
  scale_y_continuous(limits = c(0, 35)) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Year', y = '\n% Positive', col = 'Virus', tag = 'B')
p1b <- reposition_legend(p1b, position = 'top left', plot = FALSE)

p1c <- ggplot(data = dat_hk, aes(x = Time, y = GOPC)) + geom_line() +
  theme_classic() + theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          plot.tag = element_text(size = 22),
                          plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks_hk, labels = 2014:2019) +
  labs(x = 'Year', y = 'ILI Cases per 1000\nConsultations', tag = 'C')
p1d <- ggplot(data = dat_can, aes(x = time, y = i_ILI)) + geom_line() +
  theme_classic() + theme(axis.title = element_text(size = 14),
                          axis.text = element_text(size = 12),
                          plot.tag = element_text(size = 22),
                          plot.tag.position = c(0.005, 0.98)) +
  scale_x_continuous(breaks = x_lab_breaks_can, labels = 2011:2014) +
  labs(x = 'Year', y = 'ILI Cases per 1000\nConsultations', tag = 'D')

fig1 <- arrangeGrob(p1a, p1b, p1c, p1d, ncol = 2, widths = c(1, 0.666))
plot(fig1)

ggsave('results/plots/figures_for_manuscript/Figure1.svg', width = 15, height = 7, fig1)
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Figure 2: Model schematic + influence of key model parameters
# Model schematic (Figure 2a) not generated in R

# Read in MLEs:
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')[1, ]

# Set shared and unit parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Load pomp models:
prof_lik <- FALSE
vir1 <- 'flu_h1_plus_b'
fit_canada <- FALSE

source('src/functions/setup_global_likelilhood.R')

# Use 2018-19:
yr <- seasons[6]

# Get pomp object:
resp_mod <- po_list[[6]]

# Get parameter values:
pars_temp <- mle %>%
  select(all_of(shared_estpars),
         contains(yr))
names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars

# Set coefficient values:
coef(resp_mod, c(shared_estpars, unit_estpars)) <- pars_temp

# Run model using various theta_lambda1 values:
param_mat <- parmat(coef(resp_mod), nrep = 6)
param_mat['theta_lambda1', ] <- seq(0, 1, by = 0.2)
param_mat['theta_lambda2', ] <- 1.0

traj_temp <- trajectory(resp_mod, params = param_mat, format = 'data.frame') %>%
  as_tibble() %>%
  mutate(season = yr,
         theta_lambda1 = as.character(seq(0, 1, by = 0.2)[.id])) %>%
  select(time, season, H1, H2, theta_lambda1)

p2b <- ggplot(data = traj_temp) +
  geom_line(aes(x = time, y = H1, group = theta_lambda1), linetype = 2) +
  geom_line(aes(x = time, y = H2, group = theta_lambda1, col = theta_lambda1)) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.95, 0.97)) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = 'Time (Weeks)', y = 'RSV Incidence', col = expr(theta[lambda*1]), tag = 'B')
p2b <- reposition_legend(p2b, x = 0.7, y = 0.56, just = 0, plot = FALSE)

# Run model using various delta1 values:
param_mat <- parmat(coef(resp_mod), nrep = 5)
param_mat['delta1', ] <- 7 / c(7, 15, 30, 60, 182.5)
param_mat['theta_lambda1', ][1] <- 1.0
param_mat['theta_lambda2', ] <- 1.0

traj_temp <- trajectory(resp_mod, params = param_mat, format = 'data.frame') %>%
  as_tibble() %>%
  mutate(season = yr,
         delta1 = as.character(round(c(7 / c(7, 15, 30, 60, 182.5)), 3)[.id]),
         delta1 = if_else(delta1 == 1, 'No interaction', delta1)) %>%
  select(time, season, H1, H2, delta1)

traj_temp <- traj_temp %>%
  mutate(delta1 = if_else(delta1 == '0.467', '0.467 (15 days)', delta1),
         delta1 = if_else(delta1 == '0.233', '0.233 (1 month)', delta1),
         delta1 = if_else(delta1 == '0.117', '0.117 (2 months)', delta1),
         delta1 = if_else(delta1 == '0.038', '0.038 (6 months)', delta1))

p2c <- ggplot(data = traj_temp) +
  geom_line(aes(x = time, y = H1, group = delta1), linetype = 2) +
  geom_line(aes(x = time, y = H2, group = delta1, col = delta1)) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.95, 0.97)) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 10) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = 'Time (Weeks)', y = 'RSV Incidence', col = expr(delta[1]), tag = 'C')
p2c <- reposition_legend(p2c, x = 0.58, y = 0.6, just = 0, plot = FALSE)

# Combine plots and save:
fig2 <- arrangeGrob(p2b, p2c, ncol = 1)
plot(fig2)

ggsave('results/plots/figures_for_manuscript/Figure2_bc.svg', width = 5.0, height = 6.0, fig2)
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Figure 3: Observed vs. simulated cases

# Set up:
source('src/functions/functions_evalutate_res.R')

shared_estpars_hk <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                       'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
shared_estpars_can <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                        'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars_hk <- c(shared_estpars_hk, unit_estpars)
true_estpars_can <- c(shared_estpars_can, unit_estpars)

mle_hk <- read_rds('results/MLEs_flu_h1_plus_b.rds')
mle_can <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')

dat_hk <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv
dat_can <- read_csv('data/formatted/dat_canada.csv')

# Run simulations for Hong Kong:
prof_lik <- FALSE
fit_canada <- FALSE
fit_us <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk

source('src/functions/setup_global_likelilhood.R')

sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  traj_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'deterministic', return_obs = TRUE, analysis = 'iqr') %>%
    mutate(mean1 = obs1,
           mean2 = obs2,
           lower1 = qbinom(p = 0.025, size = n_T, prob = rho1_w),
           lower2 = qbinom(p = 0.025, size = n_T, prob = rho2_w),
           upper1 = qbinom(p = 0.975, size = n_T, prob = rho1_w),
           upper2 = qbinom(p = 0.975, size = n_T, prob = rho2_w)) %>%
    select(time, season, mean1:upper2)
  
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'stochastic', n_sim = 100, analysis = 'basic') %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              # median1_stoch = median(n_P1),
              # median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.025, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.025, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.975, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.975, na.rm = TRUE))
  
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  sim_temp <- traj_temp %>%
    inner_join(sim_temp, by = 'time')
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  
  sim_list[[i]] <- sim_temp
  
}

res_hk <- bind_rows(sim_list) %>%
  inner_join(dat_hk,
             by = c('time', 'season')) %>%
  select(time:Year, Week, n_P1:n_P2) %>%
  mutate(loc = 'hk')

# Run simulations for Canada:
fit_canada <- TRUE
vir1 <- 'flu'
true_estpars <- true_estpars_can

source('src/functions/setup_global_likelilhood.R')

sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  traj_temp <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'deterministic', return_obs = TRUE, analysis = 'iqr') %>%
    mutate(mean1 = obs1,
           mean2 = obs2,
           lower1 = qbinom(p = 0.025, size = n_T1, prob = rho1_w),
           lower2 = qbinom(p = 0.025, size = n_T2, prob = rho2_w),
           upper1 = qbinom(p = 0.975, size = n_T1, prob = rho1_w),
           upper2 = qbinom(p = 0.975, size = n_T2, prob = rho2_w)) %>%
    select(time, season, mean1:upper2)
  
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'stochastic', n_sim = 100, analysis = 'basic') %>%
    group_by(time) %>%
    summarise(mean1_stoch = mean(n_P1),
              mean2_stoch = mean(n_P2),
              # median1_stoch = median(n_P1),
              # median2_stoch = median(n_P2),
              lower1_stoch = quantile(n_P1, probs = 0.025, na.rm = TRUE),
              lower2_stoch = quantile(n_P2, probs = 0.025, na.rm = TRUE),
              upper1_stoch = quantile(n_P1, probs = 0.975, na.rm = TRUE),
              upper2_stoch = quantile(n_P2, probs = 0.975, na.rm = TRUE))
  
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  sim_temp <- traj_temp %>%
    inner_join(sim_temp, by = 'time')
  expect_equal(nrow(traj_temp), nrow(sim_temp))
  
  sim_list[[i]] <- sim_temp
  
}

res_can <- bind_rows(sim_list) %>%
  inner_join(dat_can,
             by = c('time', 'season')) %>%
  mutate(Year = year,
         Week = week) %>%
  select(time:upper2_stoch, Year:Week, n_P1:n_P2) %>%
  mutate(loc = 'canada')

# Combine:
res <- bind_rows(res_hk, res_can)

# Format tibble for analyses/plotting:
res <- res %>%
  select(time:season, Year:Week, loc, mean1:upper2, n_P1:n_P2) %>%
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

# For each location, virus, and season, calculate coefficient of efficiency:
# https://stats.stackexchange.com/questions/185898/difference-between-nash-sutcliffe-efficiency-and-coefficient-of-determination
r2as = r2bs = r2a_list = r2b_list = vector('list', 2)
for (i in 1:2) {
  
  location <- c('hk', 'canada')[i]
  
  r2a_list_temp = r2b_list_temp = c()
  r2a_num = r2b_num = 0
  r2a_denom = r2b_denom = 0
  
  res_loc <- res %>%
    filter(loc == location)
  
  for (yr in unique(res_loc$season)) {
    
    res_temp <- res_loc %>%
      filter(season == yr)
    
    r2_a <- res_temp %>%
      filter(virus == 'Influenza') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total)
    r2a_num <- r2a_num + r2_a$ss_error
    r2a_denom <- r2a_denom + r2_a$ss_total
    r2_a <- r2_a %>%
      pull(r2)
    
    r2_b <- res_temp %>%
      filter(virus == 'RSV') %>%
      mutate(resid_sq = (obs - mean) ** 2,
             total_sq = (obs - mean(obs, na.rm = TRUE)) ** 2) %>%
      summarise(ss_error = sum(resid_sq, na.rm = TRUE),
                ss_total = sum(total_sq, na.rm = TRUE)) %>%
      mutate(r2 = 1 - ss_error / ss_total)
    r2b_num <- r2b_num + r2_b$ss_error
    r2b_denom <- r2b_denom + r2_b$ss_total
    r2_b <- r2_b %>%
      pull(r2)
    
    r2a_list_temp <- c(r2a_list_temp, r2_a)
    r2b_list_temp <- c(r2b_list_temp, r2_b)
    
  }
  
  print(summary(r2a_list_temp))
  print(summary(r2b_list_temp))
  
  r2a <- 1 - (r2a_num / r2a_denom)
  r2b <- 1 - (r2b_num / r2b_denom)
  
  r2a_list[[i]] <- r2a_list_temp
  r2b_list[[i]] <- r2b_list_temp
  r2as[[i]] <- r2a
  r2bs[[i]] <- r2b
  
}

print(r2as)
print(r2bs)

# # Do the same, but instead of comparing to observed mean, fit sine wave to data:
# r2as = r2bs = r2a_list = r2b_list = vector('list', 2)
# for (i in 1:2) {
#   
#   location <- c('hk', 'canada')[i]
#   
#   r2a_list_temp = r2b_list_temp = c()
#   r2a_num = r2b_num = 0
#   r2a_denom = r2b_denom = 0
#   
#   res_loc <- res %>%
#     filter(loc == location)
#   
#   omega <- 2 * pi / 52.25
#   
#   for (yr in unique(res_loc$season)) {
#     
#     res_temp <- res_loc %>%
#       filter(season == yr)
#     
#     res_temp1 <- res_temp %>%
#       filter(virus == 'Influenza')
#     res_temp2 <- res_temp %>%
#       filter(virus == 'RSV')
#     
#     m1 <- lm(log(obs) ~ sin(omega * time) + cos(omega * time), data = res_temp1)
#     m2 <- lm(log(obs) ~ sin(omega * time) + cos(omega * time), data = res_temp2)
#     
#     par(mfrow = c(1, 2))
#     plot(res_temp1$obs, pch = 20, xlab = 'Time', ylab = 'Influenza Cases')
#     lines(exp(m1$coefficients[1] + m1$coefficients[2] * sin(omega * res_temp1$time) + m1$coefficients[3] * cos(omega * res_temp1$time)))
#     lines(res_temp1$mean, col = 'blue')
#     plot(res_temp2$obs, pch = 20, xlab = 'Time', ylab = 'RSV Cases')
#     lines(exp(m2$coefficients[1] + m2$coefficients[2] * sin(omega * res_temp1$time) + m2$coefficients[3] * cos(omega * res_temp1$time)))
#     lines(res_temp2$mean, col = 'blue')
#     
#     res_fit1 <- bind_cols(time = names(exp(m1$fitted.values)), sin_fit = exp(m1$fitted.values)) %>%
#       mutate(time = as.numeric(time))
#     res_fit2 <- bind_cols(time = names(exp(m2$fitted.values)), sin_fit = exp(m2$fitted.values)) %>%
#       mutate(time = as.numeric(time))
#     
#     r2_a <- res_temp1 %>%
#       left_join(res_fit1, by = 'time') %>%
#       mutate(resid_sq = (obs - mean) ** 2,
#              total_sq = (obs - sin_fit) ** 2) %>%
#       summarise(ss_error = sum(resid_sq, na.rm = TRUE),
#                 ss_total = sum(total_sq, na.rm = TRUE)) %>%
#       mutate(r2 = 1 - ss_error / ss_total)
#     r2a_num <- r2a_num + r2_a$ss_error
#     r2a_denom <- r2a_denom + r2_a$ss_total
#     r2_a <- r2_a %>%
#       pull(r2)
#     
#     r2_b <- res_temp2 %>%
#       left_join(res_fit2, by = 'time') %>%
#       mutate(resid_sq = (obs - mean) ** 2,
#              total_sq = (obs - sin_fit) ** 2) %>%
#       summarise(ss_error = sum(resid_sq, na.rm = TRUE),
#                 ss_total = sum(total_sq, na.rm = TRUE)) %>%
#       mutate(r2 = 1 - ss_error / ss_total)
#     r2b_num <- r2b_num + r2_b$ss_error
#     r2b_denom <- r2b_denom + r2_b$ss_total
#     r2_b <- r2_b %>%
#       pull(r2)
#     
#     r2a_list_temp <- c(r2a_list_temp, r2_a)
#     r2b_list_temp <- c(r2b_list_temp, r2_b)
#     
#   }
#   
#   print(summary(r2a_list_temp))
#   print(summary(r2b_list_temp))
#   
#   r2a <- 1 - (r2a_num / r2a_denom)
#   r2b <- 1 - (r2b_num / r2b_denom)
#   
#   r2a_list[[i]] <- r2a_list_temp
#   r2b_list[[i]] <- r2b_list_temp
#   r2as[[i]] <- r2a
#   r2bs[[i]] <- r2b
#   
# }
# 
# print(r2as)
# print(r2bs)

# Plot means and 95% CIs from binomial distribution:
p_legend1 <- ggplot(data = res %>% filter(loc == 'hk'),
                    aes(x = mean, y = obs, color = season)) +
  geom_point() +
  facet_wrap(~ virus) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = 'Season')
p_legend1 <- ggplotGrob(p_legend1)$grobs[[which(sapply(ggplotGrob(p_legend1)$grobs, function(x) x$name) == 'guide-box')]]

p_legend2 <- ggplot(data = res %>% filter(loc == 'canada'),
                    aes(x = mean, y = obs, color = season)) +
  geom_point() +
  facet_wrap(~ virus) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = 'Season')
p_legend2 <- ggplotGrob(p_legend2)$grobs[[which(sapply(ggplotGrob(p_legend2)$grobs, function(x) x$name) == 'guide-box')]]

label_a <- paste0(' = ', format(round(r2as[[1]], 2), nsmall = 2), ' (', round(min(r2a_list[[1]], na.rm = TRUE), 2), '-', round(max(r2a_list[[1]], na.rm = TRUE), 2), ')')
label_b <- paste0(' = ', round(r2bs[[1]], 2), ' (', round(min(r2b_list[[1]], na.rm = TRUE), 2), '-', round(max(r2b_list[[1]], na.rm = TRUE), 2), ')')
label_c <- paste0(' = ', round(r2as[[2]], 2), ' (', format(round(min(r2a_list[[2]], na.rm = TRUE), 2), nsmall = 2), '-', round(max(r2a_list[[2]], na.rm = TRUE), 2), ')')
label_d <- paste0(' = ', format(round(r2bs[[2]], 2), nsmall = 2), ' (', round(min(r2b_list[[2]], na.rm = TRUE), 2), '-', round(max(r2b_list[[2]], na.rm = TRUE), 2), ')')

p3a <- ggplot(data = res %>% filter(loc == 'hk' & virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  # geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
  #             inherit.aes = FALSE, color = 'black', se = FALSE) +
  annotate (
    geom = 'text', x = 500, y = 100, hjust = 0, vjust = 1, size = 5.5,
    label = as.expression(bquote(paste(R^2, .(label_a))))
  ) +
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
  scale_color_manual(values = viridis(6)) +
  labs(title = 'Influenza (Hong Kong)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'A')

p3b <- ggplot(data = res %>% filter(loc == 'hk' & virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  annotate (
    geom = 'text', x = 190, y = 18, hjust = 0, vjust = 1, size = 5.5,
    label = as.expression(bquote(paste(R^2, .(label_b))))
  ) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 200, 300, 400, 500)) +
  scale_color_manual(values = viridis(6)) +
  labs(title = 'RSV (Hong Kong)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'B')

p3c <- ggplot(data = res %>% filter(loc == 'canada' & virus == 'Influenza'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  # geom_smooth(aes(x = mean, y = obs), method = 'lm', formula = y ~ x,
  #             inherit.aes = FALSE, color = 'black', se = FALSE) +
  annotate (
    geom = 'text', x = 900, y = 150, hjust = 0, vjust = 1, size = 5.5,
    label = as.expression(bquote(paste(R^2, .(label_c))))
  ) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000, 3000, 4000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000, 1500, 2000, 3000, 4000)) +
  scale_color_manual(values = viridis(4)) +
  labs(title = 'Influenza (Canada)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'C')

p3d <- ggplot(data = res %>% filter(loc == 'canada' & virus == 'RSV'),
              aes(x = mean, y = obs, xmin = lower, xmax = upper, color = season)) +
  geom_abline(slope = 1, intercept = 0, col = 'gray80') +
  geom_pointrange(size = 0.2, alpha = 0.6) +
  annotate (
    geom = 'text', x = 350, y = 50, hjust = 0, vjust = 1, size = 5.5,
    label = as.expression(bquote(paste(R^2, .(label_d))))
  ) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_sqrt(breaks = c(10, 50, 100, 250, 500, 1000)) +
  scale_y_sqrt(breaks = c(10, 50, 100, 250, 500, 1000)) +
  scale_color_manual(values = viridis(4)) +
  labs(title = 'RSV (Canada)', x = 'Simulated Cases', y = 'Observed Cases',
       color = 'Season', tag = 'D')

fig3 <- arrangeGrob(arrangeGrob(arrangeGrob(p3a, p3b, ncol = 2), p_legend1, nrow = 2, heights = c(8, 1)),
                    arrangeGrob(arrangeGrob(p3c, p3d, ncol = 2), p_legend2, nrow = 2, heights = c(8, 1)),
                    nrow = 2)
plot(fig3)

ggsave('results/plots/figures_for_manuscript/Figure3.svg', fig3, width = 14, height = 10)
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Figure 4: Plot heatmaps from simulation study
file_list_hk <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_can <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'TEMPERATE', full.names = TRUE)

file_list_hk_red <- list.files('results/vaccine_simulation_study/simulations/fitcan/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_can_red <- list.files('results/vaccine_simulation_study/simulations/fitcan/', pattern = 'TEMPERATE', full.names = TRUE)

res_list <- vector('list', length(file_list_hk))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk[i]) %>% mutate(season = str_split(file_list_hk[i], '_')[[1]][5])
}
res_hk <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop',
         scenario = 'hk')

res_list <- vector('list', length(file_list_can))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_can[i]) %>% mutate(season = str_split(file_list_can[i], '_')[[1]][5])
}
res_can <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'temp',
         scenario = 'hk')

res_list <- vector('list', length(file_list_hk_red))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk_red[i]) %>% mutate(season = str_split(file_list_hk_red[i], '_')[[1]][5])
}
res_hk_red <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop',
         scenario = 'canada')

res_list <- vector('list', length(file_list_can_red))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_can_red[i]) %>% mutate(season = str_split(file_list_can_red[i], '_')[[1]][5])
}
res_can_red <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'temp',
         scenario = 'canada')

rm(i, file_list_hk, file_list_can, file_list_hk_red, file_list_can_red, res_list)

res <- res_hk %>%
  bind_rows(res_can) %>%
  bind_rows(res_hk_red) %>%
  bind_rows(res_can_red)
rm(res_hk, res_can, res_hk_red, res_can_red)

res_metrics <- res %>%
  group_by(climate, scenario, season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(climate, scenario, season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

res_metrics <- res_metrics %>%
  filter(vacc_cov <= 0.70) %>%
  mutate(vacc_cov = vacc_cov * 100)

res <- res %>%
  filter((climate == 'subtrop' & season == 's15-16') |
           (climate == 'temp' & season == 's12-13'))
res_metrics <- res_metrics %>%
  filter((climate == 'subtrop' & season == 's15-16') |
           (climate == 'temp' & season == 's12-13'))

upper_bound_ar <- max(res_metrics$ar2_impact)

p_legend1 <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
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
#   filter(climate == 'temp' & scenario == 'hk') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simA <- res %>%
  filter(climate == 'temp',
         scenario == 'hk') %>%
  filter((vacc_cov == '0.3' & vacc_time == '0') | (vacc_cov == '0.7' & vacc_time == '19')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_time = factor(vacc_time, levels = c('19', '0')))

p_legend2 <- ggplot(data = res_simA, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')
p_legend2 <- ggplotGrob(p_legend2)$grobs[[which(sapply(ggplotGrob(p_legend2)$grobs, function(x) x$name) == 'guide-box')]]

p4a <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'hk'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 19, y = 70, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 30, shape = 17, size = 3, color = 'black', fill = 'black') +
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
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0), breaks = seq(10, 70, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], ', ', delta[vacc], ' from Hong Kong)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')

p4a_sim <- ggplot(data = res_simA, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = as.numeric(as.character(vacc_time))), lty = 2) +
  geom_point(x = 52, y = 8.0, aes(shape = vacc_time), size = 3, col = 'black') +
  facet_wrap(~ vacc_time, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4b <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'hk'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 15, y = 70, shape = 16, size = 3, color = 'black', fill = 'black') +
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
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], ', ', delta[vacc], ' from Hong Kong)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'B')

# res_metrics %>%
#   filter(climate == 'subtrop' & scenario == 'hk') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simB <- res %>%
  filter(climate == 'subtrop',
         scenario == 'hk') %>%
  filter((vacc_cov == '0.4' & vacc_time == '0') | (vacc_cov == '0.7' & vacc_time == '15')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_time = factor(vacc_time, levels = c('15', '0')))

p4b_sim <- ggplot(data = res_simB, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = as.numeric(as.character(vacc_time))), lty = 2) +
  geom_point(x = 52, y = 4.25, aes(shape = vacc_time), size = 3, col = 'black', fill = 'black') +
  facet_wrap(~ vacc_time, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4c <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'canada'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 19, y = 70, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 70, shape = 17, size = 3, color = 'black', fill = 'black') +
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
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], ', ', delta[vacc], ' from Canada)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'C')

# res_metrics %>%
#   filter(climate == 'temp' & scenario == 'canada') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simC <- res %>%
  filter(climate == 'temp',
         scenario == 'canada',
         vacc_cov == 0.7,
         vacc_time %in% c('19', '0')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_time = factor(vacc_time, levels = c('19', '0')))

p4c_sim <- ggplot(data = res_simC, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = as.numeric(as.character(vacc_time))), lty = 2) +
  geom_point(x = 52, y = 8.75, aes(shape = vacc_time), size = 3, col = 'black') +
  facet_wrap(~ vacc_time, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

p4d <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'canada'),
              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  geom_point(x = 20, y = 70, shape = 16, size = 3, color = 'black', fill = 'black') +
  geom_point(x = 0, y = 70, shape = 17, size = 3, color = 'black', fill = 'black') +
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
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], ', ', delta[vacc], ' from Canada)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'D')

# res_metrics %>%
#   filter(climate == 'subtrop' & scenario == 'canada') %>%
#   filter(ar2_impact == min(ar2_impact) |
#            ar2_impact == max(ar2_impact))

res_simD <- res %>%
  filter(climate == 'subtrop',
         scenario == 'canada',
         vacc_cov == 0.7,
         vacc_time %in% c('20', '0')) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV')) %>%
  mutate(vacc_time = factor(vacc_time, levels = c('20', '0')))

p4d_sim <- ggplot(data = res_simD, aes(x = time, y = val, col = Virus, lty = .id)) +
  geom_line() +
  geom_vline(aes(xintercept = as.numeric(as.character(vacc_time))), lty = 2) +
  geom_point(x = 52, y = 3.8, aes(shape = vacc_time), size = 3, col = 'black') +
  facet_wrap(~ vacc_time, ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        strip.text = element_blank()) +
  scale_color_brewer(palette = 'Dark2') +
  scale_linetype(guide = 'none') +
  scale_shape_discrete(guide = 'none') +
  labs(title = '', x = 'Time (Weeks)', y = 'Incidence (%)')

fig4 <- arrangeGrob(arrangeGrob(arrangeGrob(p4a, p4a_sim, nrow = 1, widths = c(5.5, 4.5)),
                                arrangeGrob(p4b, p4b_sim, nrow = 1, widths = c(5.5, 4.5)), nrow = 1),
                    arrangeGrob(arrangeGrob(p4c, p4c_sim, nrow = 1, widths = c(5.5, 4.5)),
                                arrangeGrob(p4d, p4d_sim, nrow = 1, widths = c(5.5, 4.5)), nrow = 1),
                    arrangeGrob(p_legend1, p_legend2, layout_matrix = rbind(c(NA, 1, 2, NA)), widths = c(3, 2, 2, 3)),
                    nrow = 3, heights = c(12, 12, 2.25))
plot(fig4)

ggsave('results/plots/figures_for_manuscript/Figure4.svg', fig4, width = 20, height = 10)
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript (supplement)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(lemon)
library(gridExtra)
library(patchwork)
library(GGally)
library(RColorBrewer)
library(viridis)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 1: Circulation of influenza (sub)types over time

dat_hk <- read_csv('data/formatted/dat_hk.csv')

dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_b, n_rsv) %>%
  mutate(n_h1 = n_h1 / n_samp * 100,
         n_h3 = n_h3 / n_samp * 100,
         n_b = n_b / n_samp * 100,
         n_rsv = n_rsv / n_samp * 100)

dat_pos <- dat_hk %>%
  select(-n_samp, -n_rsv) %>%
  pivot_longer(n_h1:n_b,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_h1', 'n_h3', 'n_b'))) %>%
  mutate(virus = recode(virus, n_h1 = 'Influenza A(H1N1)', n_h3 = 'Influenza A(H3N2)', n_b = 'Influenza (B)'))

x_lab_breaks <- dat_hk %>% filter(Week == 1) %>% pull(Time)

fig1s <- ggplot(data = dat_pos, aes(x = Time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  scale_y_continuous(limits = c(0, 42)) +
  scale_color_brewer(palette = 'Paired') +
  labs(x = 'Year', y = '\n% Positive', col = '(Sub)type')
fig1s <- reposition_legend(fig1s, position = 'top left', plot = FALSE)

ggsave('results/plots/figures_for_manuscript/supp/FigureS1.svg', width = 9.5, height = 4, fig1s)

rm(dat_pos, fig1s)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 2: Plot of climate data and correlations

dat_clim <- read_csv('data/formatted/clim_dat_hk.csv')

dat_clim <- dat_clim %>%
  inner_join(dat_hk,
             by = c('year' = 'Year',
                    'week' = 'Week')) %>%
  select(Time, temp, ah, rh)

p2a <- ggplot(data = dat_clim, aes(x = Time, y = temp)) +
  geom_line() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.031, 0.96)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  scale_y_continuous(limits = c(9, 33)) +
  labs(x = 'Year', y = 'Mean Temperature (°C)', tag = 'A')
p2b <- ggplot(data = dat_clim, aes(x = Time, y = ah)) +
  geom_line() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.04, 0.96)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  labs(x = 'Year', y = expression(paste('Mean Absolute Humidity ', (g/m^{3}))), tag = 'B')
p2c <- ggplot(data = dat_clim, aes(x = temp, y = ah)) +
  geom_point() + theme_classic() +
  annotate (
    geom = 'text', x = 13, y = 15, label = 'r = 0.944', hjust = 0, vjust = 1, size = 6
  ) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.04, 0.96)) +
  labs(x = 'Mean Temperature (°C)', y = expression(paste('Mean Absolute Humidity ', (g/m^{3}))), tag = 'C')

fig2s <- arrangeGrob(p2a, p2b, p2c, ncol = 1)
ggsave('results/plots/figures_for_manuscript/supp/FigureS2.svg', width = 10.5, height = 9, fig2s)

rm(dat_hk, dat_clim, p2a, p2b, p2c, fig2s, x_lab_breaks)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 3: Season-specific parameters and CIs

shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

mle <- read_csv('results/MLE_plus_99CI_from_boostrapping_HPDI.csv')
res <- mle %>%
  filter(!(parameter %in% shared_estpars) & parameter != 'delta2') %>%
  mutate(season = str_sub(parameter, 2, 6),
         parameter = str_sub(parameter, 8)) %>%
  mutate(outside = mle < lower | mle > upper)
res$parameter <- factor(res$parameter, levels = c(unit_estpars, 'R10 + R120', 'R20 + R120'))

p3a <- ggplot(data = res %>% filter(vir1 == 'flu_h1')) +
  geom_errorbar(aes(x = season, ymin = lower, ymax = upper)) +
  geom_point(aes(x = season, y = mle, col = outside)) +
  facet_wrap(~ parameter, scales = 'free_y', ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.05, 0.995)) +
  scale_color_manual(values = c('black', 'darkred'), guide = 'none') +
  labs(title = 'A(H1N1)-RSV', x = 'Season',
       y = 'Maximum Likelihood Estimate (99% CI)', tag = 'A')
p3b <- ggplot(data = res %>% filter(vir1 == 'flu_b')) +
  geom_errorbar(aes(x = season, ymin = lower, ymax = upper)) +
  geom_point(aes(x = season, y = mle, col = outside)) +
  facet_wrap(~ parameter, scales = 'free_y', ncol = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.05, 0.995)) +
  scale_color_manual(values = c('black', 'darkred'), guide = 'none') +
  labs(title = 'B-RSV', x = 'Season',
       y = 'Maximum Likelihood Estimate (99% CI)', tag = 'B')

fig3s <- arrangeGrob(p3a, p3b, ncol = 2)
ggsave('results/plots/figures_for_manuscript/supp/FigureS3.svg', width = 7.2, height = 21.5, fig3s)

rm(p3a, p3b, fig3s, res, mle)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 4: Correlations between estimated parameters

res_dir_h1 <- 'results/round2_4_fluH1_FULL/'
res_dir_b <- 'results/round2_3_fluB_FULL/'

load_and_format_mega_results <- function(filepath, shared_estpars, unit_estpars, run_name) {
  
  # Compile estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full = list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get parameter estimates and log-likelihoods:
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  expect_true(nrow(pars_df) == length(res_files))
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Check whether estimations stopped due to tolerance or time:
  print(table(pars_df$message))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  expect_equal(df_use, 47)
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Return formatted results:
  return(pars_top)
  
}

res_h1 <- load_and_format_mega_results(filepath = res_dir_h1,
                                       shared_estpars = shared_estpars,
                                       unit_estpars = unit_estpars,
                                       run_name = 'H1_FULL')
res_b <- load_and_format_mega_results(filepath = res_dir_b,
                                      shared_estpars = shared_estpars,
                                      unit_estpars = unit_estpars,
                                      run_name = 'B_FULL')

res_LIST <- list(res_h1, res_b)
pars_top_LIST <- vector('list', length = length(res_LIST))
for (i in 1:length(res_LIST)) {
  pars_top_LIST[[i]] <- res_LIST[[i]]
}
rm(i)
names(pars_top_LIST) <- c('flu_h1_FULL', 'flu_b_FULL')

rm(res_LIST, res_h1, res_b, res_dir_h1, res_dir_b)

pars_top_LIST_temp <- pars_top_LIST

names(pars_top_LIST_temp[[1]])[1:12] <- c('rho[1]', 'rho[2]', 'theta[lambda*1]', 'theta[lambda*2]', 'delta[1]', 'd[2]', 'alpha', 'phi', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]')
names(pars_top_LIST_temp[[2]])[1:12] <- c('rho[1]', 'rho[2]', 'theta[lambda*1]', 'theta[lambda*2]', 'delta[1]', 'd[2]', 'alpha', 'phi', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]')
shared_estpars_temp <- c('rho[1]', 'rho[2]', 'theta[lambda*1]', 'theta[lambda*2]', 'delta[1]', 'd[2]', 'alpha', 'phi', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]')

fig4sa <- ggpairs(pars_top_LIST_temp[[1]] %>% select(all_of(shared_estpars_temp)),
                  upper = list(continuous = wrap(ggally_cor, size = 4.5, method = 'kendall', digits = 2, display_grid = FALSE)),
                  lower = list(continuous = wrap('points', size = 1.1)),
                  labeller = 'label_parsed') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 55, vjust = 0.6),
        strip.text = element_text(size = 15),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.005, 0.988),
        panel.border = element_rect(size = 0.5, fill = NA)) +
  labs(tag = 'A')

fig4sb <- ggpairs(pars_top_LIST_temp[[2]] %>% select(all_of(shared_estpars_temp)),
                  upper = list(continuous = wrap(ggally_cor, size = 4.5, method = 'kendall', digits = 2, display_grid = FALSE)),
                  labeller = 'label_parsed') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 55, vjust = 0.6),
        strip.text = element_text(size = 15),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.005, 0.988),
        panel.border = element_rect(size = 0.5, fill = NA)) +
  labs(tag = 'B')

ggsave('results/plots/figures_for_manuscript/supp/FigureS4a.svg', fig4sa, width = 18, height = 12)
ggsave('results/plots/figures_for_manuscript/supp/FigureS4b.svg', fig4sb, width = 18, height = 12)

rm(fig4sa, fig4sb, pars_top_LIST_temp, shared_estpars_temp)

# ---------------------------------------------------------------------------------------------------------------------

# # EXTRA Supplementary Figure: Beta due to climate over time
# 
# gamma1 <- 7/5
# gamma2 <- 7/10
#
# mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
# mle_b <- read_rds('results/MLEs_flu_b.rds')
# 
# mles <- bind_rows(mle_h1[1, ], mle_b[1, ]) %>%
#   mutate(vir1 = c('flu_h1', 'flu_b')) %>%
#   select(vir1,
#          eta_temp1:eta_ah2,
#          contains('R10'),
#          contains('R20'),
#          contains('R120'),
#          contains('Ri')) %>%
#   pivot_longer(-c(vir1, eta_temp1, eta_temp2, eta_ah1, eta_ah2),
#                names_to = 'parameter',
#                values_to = 'mle') %>%
#   mutate(season = str_sub(parameter, 2, 6),
#          parameter = str_sub(parameter, 8)) %>%
#   pivot_wider(names_from = parameter, values_from = mle) %>%
#   select(vir1, eta_temp1:eta_ah2, R10:Ri2, season) %>%
#   filter(!is.na(R10))
# 
# seasons_h1 <- c('13-14', '15-16', '16-17', '17-18', '18-19')
# seasons_b <- c('13-14', '14-15', '15-16', '17-18', '18-19')
# 
# hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
# dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
# 
# beta_h1 <- vector('list', length = length(seasons_h1))
# for (seas_index in 1:length(seasons_h1)) {
#   
#   seas <- seasons_h1[seas_index]
#   
#   dat_temp <- hk_dat[['h1_rsv']] %>%
#     filter(season == paste0('s', seas)) %>%
#     inner_join(dat_clim,
#                by = c('Year' = 'year',
#                       'Week' = 'week')) %>%
#     filter(Week != 53) %>%
#     select(time, temp, ah)
#   
#   mle_temp <- mles %>%
#     filter(vir1 == 'flu_h1',
#            season == seas)
#   
#   beta1_temp <- unlist(mle_temp['Ri1']) / (1.0 - (unlist(mle_temp['R10']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp) * gamma1
#   beta2_temp <- unlist(mle_temp['Ri2']) / (1.0 - (unlist(mle_temp['R20']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp) * gamma2
#   
#   # beta1_temp <- bind_cols(1:52, (beta1_temp - mean(beta1_temp)) / sd(beta1_temp), seas)
#   # beta2_temp <- bind_cols(1:52, (beta2_temp - mean(beta2_temp)) / sd(beta2_temp), seas)
#   beta1_temp <- bind_cols(1:52, beta1_temp / mean(beta1_temp), seas)
#   beta2_temp <- bind_cols(1:52, beta2_temp / mean(beta2_temp), seas)
#   
#   names(beta1_temp) <- c('time', 'beta1', 'season')
#   names(beta2_temp) <- c('time', 'beta2', 'season')
#   
#   beta_temp <- beta1_temp %>%
#     inner_join(beta2_temp,
#                by = c('time', 'season')) %>%
#     select(time, beta1, beta2, season)
#   
#   beta_h1[[seas_index]] <- beta_temp
#   
# }
# rm(seas_index, seas, dat_temp, mle_temp, beta1_temp, beta2_temp)
# 
# beta_b <- vector('list', length = length(seasons_b))
# for (seas_index in 1:length(seasons_b)) {
#   
#   seas <- seasons_b[seas_index]
#   
#   dat_temp <- hk_dat[['b_rsv']] %>%
#     filter(season == paste0('s', seas)) %>%
#     inner_join(dat_clim,
#                by = c('Year' = 'year',
#                       'Week' = 'week')) %>%
#     filter(Week != 53) %>%
#     select(time, temp, ah)
#   
#   mle_temp <- mles %>%
#     filter(vir1 == 'flu_b',
#            season == seas)
#   
#   beta1_temp <- unlist(mle_temp['Ri1']) / (1.0 - (unlist(mle_temp['R10']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp) * gamma1
#   beta2_temp <- unlist(mle_temp['Ri2']) / (1.0 - (unlist(mle_temp['R20']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp) * gamma2
#   
#   # beta1_temp <- bind_cols(1:52, (beta1_temp - mean(beta1_temp)) / sd(beta1_temp), seas)
#   # beta2_temp <- bind_cols(1:52, (beta2_temp - mean(beta2_temp)) / sd(beta2_temp), seas)
#   beta1_temp <- bind_cols(1:52, beta1_temp / mean(beta1_temp), seas)
#   beta2_temp <- bind_cols(1:52, beta2_temp / mean(beta2_temp), seas)
#   
#   names(beta1_temp) <- c('time', 'beta1', 'season')
#   names(beta2_temp) <- c('time', 'beta2', 'season')
#   
#   beta_temp <- beta1_temp %>%
#     inner_join(beta2_temp,
#                by = c('time', 'season')) %>%
#     select(time, beta1, beta2, season)
#   
#   beta_b[[seas_index]] <- beta_temp
#   
# }
# rm(seas_index, seas, dat_temp, mle_temp, beta1_temp, beta2_temp)
# 
# beta_h1 <- bind_rows(beta_h1)
# beta_b <- bind_rows(beta_b)
# 
# beta_temp <- beta_h1 %>%
#   bind_rows(beta_b)
# 
# min_val <- min(c(min(beta_temp$beta1), min(beta_temp$beta2)))
# max_val <- max(c(max(beta_temp$beta1), max(beta_temp$beta2)))
# 
# p_legend <- ggplot(data = beta_temp, aes(x = time, y = beta1, col = season)) +
#   geom_line() +
#   theme_classic() +
#   theme(legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12),
#         legend.position = 'bottom') +
#   guides(color = guide_legend(nrow = 1)) +
#   scale_color_viridis(discrete = TRUE) +
#   labs(color = 'Season')
# p_legend <- ggplotGrob(p_legend)$grobs[[which(sapply(ggplotGrob(p_legend)$grobs, function(x) x$name) == 'guide-box')]]
# 
# p_a <- ggplot(data = beta_h1, aes(x = time, y = beta1, col = season)) +
#   geom_line() +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.position = 'none',
#         plot.tag = element_text(size = 22),
#         plot.tag.position = c(0.02, 0.97)) +
#   scale_x_continuous(breaks = seq(1, 52, by = 5),
#                      labels = c(46, 51, seq(4, 45, by = 5))) +
#   # scale_y_continuous(limits = c(min_val, max_val)) +
#   scale_y_log10(limits = c(min_val, max_val), breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
#   scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
#   labs(x = 'Week Number', y = expression(beta[1] * ' (Relative to Mean)'), tag = 'A')
# p_b <- ggplot(data = beta_h1, aes(x = time, y = beta2, col = season)) +
#   geom_line() +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.position = 'none',
#         plot.tag = element_text(size = 22),
#         plot.tag.position = c(0.02, 0.97)) +
#   scale_x_continuous(breaks = seq(1, 52, by = 5),
#                      labels = c(46, 51, seq(4, 45, by = 5))) +
#   scale_y_log10(limits = c(min_val, max_val), breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
#   scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
#   labs(x = 'Week Number', y = expression(beta[2] * ' (Relative to Mean)'), tag = 'B')
# p_c <- ggplot(data = beta_b, aes(x = time, y = beta1, col = season)) +
#   geom_line() +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.position = 'none',
#         plot.tag = element_text(size = 22),
#         plot.tag.position = c(0.02, 0.97)) +
#   scale_x_continuous(breaks = seq(1, 52, by = 5),
#                      labels = c(46, 51, seq(4, 45, by = 5))) +
#   scale_y_log10(limits = c(min_val, max_val), breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
#   scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
#   labs(x = 'Week Number', y = expression(beta[1] * ' (Relative to Mean)'), tag = 'C')
# p_d <- ggplot(data = beta_b, aes(x = time, y = beta2, col = season)) +
#   geom_line() +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         legend.position = 'none',
#         plot.tag = element_text(size = 22),
#         plot.tag.position = c(0.02, 0.97)) +
#   scale_x_continuous(breaks = seq(1, 52, by = 5),
#                      labels = c(46, 51, seq(4, 45, by = 5))) +
#   scale_y_log10(limits = c(min_val, max_val), breaks = c(0.6, 0.8, 1.0, 1.2, 1.4)) +
#   scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
#   labs(x = 'Week Number', y = expression(beta[2] * ' (Relative to Mean)'), tag = 'D')
# 
# fig_s_extra <- arrangeGrob(arrangeGrob(p_a, p_b, p_c, p_d, ncol = 2), p_legend, nrow = 2, heights = c(15, 1))
# # ggsave('results/plots/figures_for_manuscript/supp/FigureS_EXTRA.svg', fig6s, width = 15, height = 8.5)
# 
# rm(mles, beta_h1, beta_b, beta_temp, hk_dat, dat_clim, p_a, p_b, p_c, p_d, p_legend, fig_s_extra,
#    seasons_h1, seasons_b, min_val, max_val, gamma1, gamma2)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 5: Simulations at the MLE

true_estpars <- c(shared_estpars, unit_estpars)
prof_lik <- FALSE

vir1 <- 'flu_h1'
source('src/functions/setup_global_likelilhood.R')

set.seed(12075)
sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- pars_top_LIST[[1]] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp[1, ]
  sim_temp <- simulate(resp_mod, nsim = 10, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    select(time:.id, n_P1:n_P2) %>%
    arrange(.id) %>%
    cbind(t(resp_mod@data))
  names(sim_temp)[5:6] <- c('obs1', 'obs2')
  sim_temp <- sim_temp %>%
    mutate(season = yr) %>%
    as_tibble()
  
  sim_list[[i]] <- sim_temp
  
}

sim_h1 <- bind_rows(sim_list)

vir1 <- 'flu_b'
source('src/functions/setup_global_likelilhood.R')

set.seed(12075)
sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- pars_top_LIST[[2]] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp[1, ]
  sim_temp <- simulate(resp_mod, nsim = 10, format = 'data.frame')
  
  sim_temp <- sim_temp %>%
    select(time:.id, n_P1:n_P2) %>%
    arrange(.id) %>%
    cbind(t(resp_mod@data))
  names(sim_temp)[5:6] <- c('obs1', 'obs2')
  sim_temp <- sim_temp %>%
    mutate(season = yr) %>%
    as_tibble()
  
  sim_list[[i]] <- sim_temp
  
}

sim_b <- bind_rows(sim_list)

rm(sim_list, po_list, vir1, yr, resp_mod, pars_temp, sim_temp, i,
   dat_pomp, hk_dat, nrow_check, yr_index)

sim_h1 <- sim_h1 %>%
  pivot_longer(n_P1:n_P2, names_to = 'vir1', values_to = 'sim') %>%
  pivot_longer(obs1:obs2, names_to = 'vir2', values_to = 'obs') %>%
  mutate(vir1 = if_else(vir1 == 'n_P1', 'Influenza (H1)', 'RSV')) %>%
  mutate(vir2 = if_else(vir2 == 'obs1', 'Influenza (H1)', 'RSV')) %>%
  filter(vir1 == vir2) %>%
  mutate(vir = vir1) %>%
  select(time:.id, season, vir1, sim, obs) %>%
  mutate(season = str_sub(season, 2))
sim_b <- sim_b %>%
  pivot_longer(n_P1:n_P2, names_to = 'vir1', values_to = 'sim') %>%
  pivot_longer(obs1:obs2, names_to = 'vir2', values_to = 'obs') %>%
  mutate(vir1 = if_else(vir1 == 'n_P1', 'Influenza (B)', 'RSV')) %>%
  mutate(vir2 = if_else(vir2 == 'obs1', 'Influenza (B)', 'RSV')) %>%
  filter(vir1 == vir2) %>%
  mutate(vir = vir1) %>%
  select(time:.id, season, vir1, sim, obs) %>%
  mutate(season = str_sub(season, 2))

breaks_fxn <- function(time) {
  if (max(time) > 55) {
    c(1, 6, seq(12, 52, by = 5))
    # seq(1, 52, by = 10)
  } else {
    seq(1, 52, by = 5)
  }
}
# https://coolbutuseless.github.io/2019/03/07/custom-axis-breaks-on-facetted-ggplot/

p5a <- ggplot(data = sim_h1) +
  geom_line(aes(x = time, y = sim, group = paste(.id, vir1), col = vir1), lwd = 0.3) +
  geom_point(aes(x = time, y = obs, group = paste(.id, vir1), col = vir1), size = 0.75) +
  facet_wrap(~ season, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.0035, 0.97)) +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  labs(x = 'Week #', y = '# of Cases', col = 'Virus', tag = 'A')

p5b <- ggplot(data = sim_b) +
  geom_line(aes(x = time, y = sim, group = paste(.id, vir1), col = vir1), lwd = 0.3) +
  geom_point(aes(x = time, y = obs, group = paste(.id, vir1), col = vir1), size = 0.75) +
  facet_wrap(~ season, scales = 'free', nrow = 1) +
  theme_classic() +
  theme(legend.position = 'right',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.tag.position = c(0.0035, 0.97)) +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(2, 3)]) +
  labs(x = 'Week #', y = '# of Cases', col = 'Virus', tag = 'B')

fig5s <- arrangeGrob(p5a, p5b, ncol = 1)
ggsave('results/plots/figures_for_manuscript/supp/FigureS5.svg', width = 21, height = 7.5, fig5s)

rm(fig5s, p5a, p5b, sim_h1, sim_b)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 6: Profile likelihood of theta_lambda1

load_and_format_proflik_results <- function(filepath, prof_par, shared_estpars) {
  
  # Remove prof_par from shared_estpars:
  shared_estpars <- shared_estpars[!(shared_estpars == prof_par)]
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full <- list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get estimated parameter and log-likelihood values:
  res_temp <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    select(all_of(shared_estpars)) %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist()) %>%
    bind_cols(map_chr(str_split(res_files, '_'), 5),
              map_chr(str_split(res_files, '_'), 7),
              map_chr(str_split(map_chr(str_split(res_files, '_'), 8), fixed('.')), 1)) %>%
    rename(vir1 = '...14',
           profpar = '...15',
           run = '...16') %>%
    mutate(profpar = as.numeric(profpar),
           run = as.numeric(run),
           vir1 = if_else(vir1 == 'b', 'B', 'H1'),
           vir1 = factor(vir1),
           vir1 = relevel(vir1, ref = 'H1')) %>%
    arrange(vir1, profpar, run)
  expect_true(nrow(res_temp) == length(res_files))
  expect_true(all(is.finite(res_temp$loglik)))
  
  res_temp <- res_temp %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # Set profpar to correct values:
  res_temp <- res_temp %>%
    mutate(profpar = seq(0.0, 0.2, by = 0.01)[profpar])
  
  # Return formatted results:
  return(res_temp)
  
}

res_proflik <- load_and_format_proflik_results(filepath = 'results/prof_lik_thetalambda1/',
                                               prof_par = 'theta_lambda1',
                                               shared_estpars = shared_estpars)

maxloglik <- res_proflik %>%
  group_by(vir1) %>%
  summarise(loglik = max(loglik)) %>%
  pull(loglik) %>%
  unlist()

ci_cutoff <- maxloglik - 0.5 * qchisq(df = 1, p = 0.99)
res_proflik <- res_proflik %>%
  mutate(ci = if_else(vir1 == 'H1', ci_cutoff[1], ci_cutoff[2]))

res_proflik <- res_proflik %>%
  group_by(vir1, profpar) %>%
  filter(rank(-loglik) == 1) %>%
  ungroup()

p6a <- ggplot(res_proflik %>%
                filter(vir1 == 'H1'),# %>%
              # filter(loglik > (max(loglik) - 100)), # allows for zooming
              aes(x = profpar, y = loglik)) +
  geom_point() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.975)) +
  geom_smooth(method = 'loess', span = 0.75, color = 'black') +
  geom_hline(color = 'black', aes(yintercept = ci), size = 1, lty = 2) +
  labs(x = bquote(theta[lambda*1]), y = 'Log-Likelihood', tag = 'A') +
  scale_x_continuous(n.breaks = 10)
p6b <- ggplot(res_proflik %>%
                filter(vir1 == 'B'),# %>%
              # filter(loglik > (max(loglik) - 100)), # allows for zooming
              aes(x = profpar, y = loglik)) +
  geom_point() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.975)) +
  geom_smooth(method = 'loess', span = 0.75, color = 'black') +
  geom_hline(color = 'black', aes(yintercept = ci), size = 1, lty = 2) +
  labs(x = bquote(theta[lambda*1]), y = 'Log-Likelihood', tag = 'B') +
  scale_x_continuous(n.breaks = 10) + scale_y_continuous(breaks = seq(-3975, -3935, by = 10))

fig6s <- arrangeGrob(p6a, p6b, nrow = 1)
ggsave('results/plots/figures_for_manuscript/supp/FigureS6.svg', width = 11.5, height = 4.2, fig6s)

rm(fig6s, p6a, p6b, res_proflik, maxloglik, ci_cutoff)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 7: Simulated attack rates for flu and RSV at the MLE

mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
mle_b <- read_rds('results/MLEs_flu_b.rds')

vir1 <- 'flu_h1'
source('src/functions/setup_global_likelilhood.R')
dat_temp <- hk_dat$h1_rsv

traj_list_H1N1 = ar_list_seas_h1 = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- mle_h1[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  param_mat_temp <- parmat(coef(resp_mod), nrep = 3)
  param_mat_temp['I10', 2] <- 0
  param_mat_temp['I20', 3] <- 0
  
  traj_temp <- trajectory(resp_mod, params = param_mat_temp, format = 'data.frame') %>%
    mutate(season = yr) %>%
    select(time:season, H1:H2)
  
  traj_list_H1N1[[i]] <- traj_temp
  
  traj_temp <- traj_temp %>%
    filter(.id == 1) %>%
    select(-.id) %>%
    left_join(dat_temp, by = c('time', 'season')) %>%
    rename('i_ILI' = 'GOPC') %>%
    mutate(i_ILI = i_ILI / 1000)
  
  ar_tot <- traj_temp %>%
    select(time, H1:H2) %>%
    summarise(H1 = sum(H1),
              H2 = sum(H2)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
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
  
  ar_obs <- traj_temp %>%
    mutate(obs1 = rho1_w * n_T,
           obs2 = rho2_w * n_T) %>%
    summarise(obs1 = sum(obs1, na.rm = TRUE),
              obs2 = sum(obs2, na.rm = TRUE)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
  ar_list_seas_h1[[i]] <- inner_join(ar_tot, ar_obs, by = c('virus1', 'season')) %>%
    select(virus1:season, H1:H2, obs1:obs2)
  
}

vir1 <- 'flu_b'
source('src/functions/setup_global_likelilhood.R')
dat_temp <- hk_dat$b_rsv

traj_list_B = ar_list_seas_b = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- mle_b[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  param_mat_temp <- parmat(coef(resp_mod), nrep = 3)
  param_mat_temp['I10', 2] <- 0
  param_mat_temp['I20', 3] <- 0
  
  traj_temp <- trajectory(resp_mod, params = param_mat_temp, format = 'data.frame') %>%
    mutate(season = yr) %>%
    select(time:season, H1:H2)
  
  traj_list_B[[i]] <- traj_temp
  
  traj_temp <- traj_temp %>%
    filter(.id == 1) %>%
    select(-.id) %>%
    left_join(dat_temp, by = c('time', 'season')) %>%
    rename('i_ILI' = 'GOPC') %>%
    mutate(i_ILI = i_ILI / 1000)
  
  ar_tot <- traj_temp %>%
    select(time, H1:H2) %>%
    summarise(H1 = sum(H1),
              H2 = sum(H2)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
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
  
  ar_obs <- traj_temp %>%
    mutate(obs1 = rho1_w * n_T,
           obs2 = rho2_w * n_T) %>%
    summarise(obs1 = sum(obs1, na.rm = TRUE),
              obs2 = sum(obs2, na.rm = TRUE)) %>%
    mutate(virus1 = vir1,
           season = yr)
  
  ar_list_seas_b[[i]] <- inner_join(ar_tot, ar_obs, by = c('virus1', 'season')) %>%
    select(virus1:season, H1:H2, obs1:obs2)
  
}

rm(i, yr, resp_mod, pars_temp, param_mat_temp, traj_temp, ar_tot, ar_obs, rho1, rho2, alpha, phi, rho1_w, rho2_w,
   dat_temp, hk_dat, dat_pomp, po_list, nrow_check, yr_index, obj_fun_list, mle_h1, mle_b)

ar_h1 <- bind_rows(ar_list_seas_h1)
ar_b <- bind_rows(ar_list_seas_b)

ar_df <- bind_rows(ar_h1, ar_b)
rm(ar_h1, ar_b, ar_list_seas_h1, ar_list_seas_b)

ar_df <- ar_df %>%
  select(-c(obs1:obs2)) %>%
  pivot_longer(H1:H2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'H1', 'Influenza', 'RSV'),
         virus1 = if_else(virus1 == 'flu_h1', 'H1', 'B'),
         virus = if_else(virus == 'Influenza' & virus1 == 'H1', 'Influenza (H1)', virus),
         virus = if_else(virus == 'Influenza' & virus1 == 'B', 'Influenza (B)', virus),
         attack_rate = attack_rate * 100)

p_legend <- ggplot(data = ar_df, aes(x = virus, y = attack_rate, group = virus)) +
  geom_violin(fill = 'gray90') +
  geom_point(aes(col = season), size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'bottom',
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_viridis(discrete = TRUE) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(col = 'Season')
p_legend <- ggplotGrob(p_legend)$grobs[[which(sapply(ggplotGrob(p_legend)$grobs, function(x) x$name) == 'guide-box')]]

p7a <- ggplot(data = ar_df %>% filter(virus1 == 'H1'),
              aes(x = virus, y = attack_rate, group = virus)) +
  geom_violin(fill = 'gray90') +
  geom_point(aes(col = season), size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.position = 'none',
        plot.tag.position = c(0.015, 0.97)) +
  scale_color_manual(values = viridis(6)[c(1, 3:6)]) +
  labs(x = 'Virus', y = 'Attack Rate (%)', tag = 'A')
p7b <- ggplot(data = ar_df %>% filter(virus1 == 'B'),
              aes(x = virus, y = attack_rate, group = virus)) +
  geom_violin(fill = 'gray90') +
  geom_point(aes(col = season), size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        legend.position = 'none',
        plot.tag.position = c(0.015, 0.97)) +
  scale_color_manual(values = viridis(6)[c(1:3, 5:6)]) +
  labs(x = 'Virus', y = 'Attack Rate (%)', tag = 'B')

fig7s <- arrangeGrob(arrangeGrob(p7a, p7b, nrow = 1), p_legend, nrow = 2, heights = c(15, 1))
ggsave('results/plots/figures_for_manuscript/supp/FigureS7.svg', width = 10, height = 4, fig7s)

rm(fig7s, p7a, p7b, ar_df, shared_estpars, unit_estpars, true_estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 8: Simulations at the MLE, with and without interaction effect

res_H1N1 <- bind_rows(traj_list_H1N1) %>%
  mutate(virus_pair = 'A(H1N1)-RSV') %>%
  as_tibble()
res_B <- bind_rows(traj_list_B) %>%
  mutate(virus_pair = 'B-RSV') %>%
  as_tibble()
res <- bind_rows(res_H1N1, res_B)

rm(res_H1N1, res_B, traj_list_H1N1, traj_list_B)

res <- res %>%
  mutate(virus_pair = factor(virus_pair, levels = c('A(H1N1)-RSV', 'B-RSV'))) %>%
  mutate(season = str_sub(season, 2))

p8a <- ggplot() +
  geom_line(data = res %>%
              filter(.id %in% c(1, 2)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H2, col = .id)) +
  geom_line(data = res %>%
              filter(.id == 1),
            aes(x = time, y = H1), lty = 2) +
  facet_grid(virus_pair ~ season, scales = 'free_x') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.004, 0.975)) +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_manual(values = c('#3182bd', '#9ecae1')) +
  labs(x = 'Week #', y = 'RSV Incidence', color = '', tag = 'A')
p8b <- ggplot() +
  geom_line(data = res %>%
              filter(virus_pair == 'A(H1N1)-RSV') %>%
              filter(.id %in% c(1, 3)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H1, col = .id)) +
  geom_line(data = res %>%
              filter(virus_pair == 'A(H1N1)-RSV') %>%
              filter(.id == 1),
            aes(x = time, y = H2), lty = 2) +
  facet_grid(virus_pair ~ season, scales = 'free_x') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.004, 0.975)) +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_manual(values = c('#de2d26', '#fc9272')) +
  labs(x = 'Week #', y = 'Influenza Incidence', color = '', tag = 'B')

fig8s <- arrangeGrob(p8a, p8b, ncol = 1, heights = c(1.6, 1))
ggsave('results/plots/figures_for_manuscript/supp/FigureS8.svg', width = 16, height = 9.5, fig8s)

rm(p8a, p8b, fig8s, res, pars_top_LIST, seasons, vir1, vir2, age_structured,
   d2_max, debug_bool, lag_val, prof_lik, Ri_max1, Ri_max2)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 9: Model schematic for vaccine simulation study
# Not generated in R

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 10: Impact of LAIV by season for all scenarios

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

res_metrics <- res_metrics %>%
  filter(!(climate == 'subtrop' & season == 's18-19') &
           !(climate == 'temp' & season == 's17-18'))

upper_bound_ar <- max(res_metrics$ar2_impact)

p10a <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 2) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.75, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], ' = ', theta[lambda*1], ')')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')

p10b <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'natural'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 2) +
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
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], ' = ', theta[lambda*1], ')')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'B')

p10c <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'half'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 2) +
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
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], ' = 0.5)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'C')

p10d <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'half'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 2) +
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
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', theta[lambda[vacc]], ' = 0.5)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'D')

fig10s <- (p10a + p10b) / (p10c + p10d) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave('results/plots/figures_for_manuscript/supp/FigureS10.svg', fig10s, width = 14, height = 10)

rm(fig10s, p10a, p10b, p10c, p10d, res, res_metrics, upper_bound_ar)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 11: Sensitivity analyses for the vaccine simulation study

file_list_hk_deltaShort <- list.files('results/vaccine_simulation_study/simulations/deltavacc1month/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp_deltaShort <- list.files('results/vaccine_simulation_study/simulations/deltavacc1month/', pattern = 'TEMPERATE', full.names = TRUE)

file_list_hk_deltaLong <- list.files('results/vaccine_simulation_study/simulations/deltavacc6months/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp_deltaLong <- list.files('results/vaccine_simulation_study/simulations/deltavacc6months/', pattern = 'TEMPERATE', full.names = TRUE)

file_list_hk_effLow <- list.files('results/vaccine_simulation_study/simulations/vacceff60/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp_effLow <- list.files('results/vaccine_simulation_study/simulations/vacceff60/', pattern = 'TEMPERATE', full.names = TRUE)

file_list_hk_effHigh <- list.files('results/vaccine_simulation_study/simulations/vacceff95/', pattern = 'SUBTROPICAL', full.names = TRUE)
file_list_temp_effHigh <- list.files('results/vaccine_simulation_study/simulations/vacceff95/', pattern = 'TEMPERATE', full.names = TRUE)

all_files_LIST <- list(file_list_hk_deltaShort, file_list_temp_deltaShort,
                       file_list_hk_deltaLong, file_list_temp_deltaLong,
                       file_list_hk_effLow, file_list_temp_effLow,
                       file_list_hk_effHigh, file_list_temp_effHigh)
res_df_LIST <- vector('list', length = length(all_files_LIST))
climates <- c('subtrop', 'temp', 'subtrop', 'temp', 'subtrop', 'temp', 'subtrop', 'temp')
scenarios <- c('deltaShort', 'deltaShort', 'deltaLong', 'deltaLong', 'effLow', 'effLow', 'effHigh', 'effHigh')

for (i in 1:length(res_df_LIST)) {
  res_list_temp <- vector('list', length = length(all_files_LIST[[i]]))
  
  for (j in 1:length(res_list_temp)) {
    res_list_temp[[j]] <- read_rds(all_files_LIST[[i]][j]) %>% mutate(season = str_split(all_files_LIST[[i]][j], '_')[[1]][5])
  }
  
  res_temp <- bind_rows(res_list_temp) %>%
    as_tibble() %>%
    mutate(climate = climates[i],
           scenario = scenarios[i])
  
  res_df_LIST[[i]] <- res_temp
}

res <- bind_rows(res_df_LIST)

rm(file_list_hk_deltaShort, file_list_temp_deltaShort, file_list_hk_deltaLong, file_list_temp_deltaLong,
   file_list_hk_effLow, file_list_temp_effLow, file_list_hk_effHigh, file_list_temp_effHigh, all_files_LIST,
   climates, scenarios, res_list_temp, res_temp, i, j, res_df_LIST)

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

res_metrics <- res_metrics %>%
  filter((climate == 'subtrop' & season == 's18-19') |
           (climate == 'temp' & season == 's17-18'))

upper_bound_ar <- max(res_metrics$ar2_impact)

p11a <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'deltaShort'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.75, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', delta[vacc], ' = 1 month)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')

p11b <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'deltaShort'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', delta[vacc], ' = 1 month)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'B')

p11c <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'deltaLong'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', delta[vacc], ' = 6 months)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'C')

p11d <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'deltaLong'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Subtropical (', delta[vacc], ' = 6 months)')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'D')

p11e <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'effLow'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = 'Temperate (VE = 60%)', x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'E')

p11f <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'effLow'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = 'Subtropical (VE = 60%)', x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'F')

p11g <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'effHigh'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = 'Temperate (VE = 95%)', x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'G')

p11h <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'effHigh'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.985)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25)),
                       guide = 'none') +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = 'Subtropical (VE = 95%)', x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'H')

fig11s <- (p11a + p11b) / (p11c + p11d) / (p11e + p11f) / (p11g + p11h) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
ggsave('results/plots/figures_for_manuscript/supp/FigureS11.svg', fig11s, width = 10, height = 15)

rm(fig11s, p11a, p11b, p11c, p11d, p11e, p11f, p11g, p11h, res, res_metrics, upper_bound_ar)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 12: Age-structured synthetic data

res_all_ages <- read_csv('data/age_structured_SA/synthetic_obs_by_age.csv')
covar_all_ages <- read_csv('data/age_structured_SA/synthetic_covariate_data.csv')

covar_all_ages <- covar_all_ages %>%
  select(time, season, n_T1:n_T5) %>%
  pivot_longer(n_T1:n_T5, names_to = 'age', values_to = 'n_samp') %>%
  mutate(age = as.numeric(str_sub(age, -1)))

res_all_ages <- res_all_ages %>%
  inner_join(covar_all_ages, by = c('time', 'season', 'age'))

check <- res_all_ages %>%
  group_by(time, season) %>%
  summarise(n_T = unique(n_T), n_samp = sum(n_samp)) %>%
  mutate(check = n_T - n_samp) %>%
  drop_na() %>%
  pull(check)
expect_true(all(check == 0))
rm(check)

res_all_ages <- res_all_ages %>%
  mutate(obs1 = obs1_s1b / n_samp * 100,
         obs2 = obs2_s1b / n_samp * 100) %>%
  select(time, season, age, obs1:obs2) %>%
  pivot_longer(obs1:obs2, names_to = 'virus', values_to = 'val') %>%
  mutate(virus = if_else(virus == 'obs1', 'Influenza', 'RSV'),
         age = as.character(age),
         age = recode(age,
                      '1' = '<1 year',
                      '2' = '1-4 years',
                      '3' = '5-15 years',
                      '4' = '16-64 years',
                      '5' = '65+ years'),
         age = factor(age, levels = c('<1 year', '1-4 years', '5-15 years', '16-64 years', '65+ years'))) %>%
  filter(!(season == 's16-17' & time == 8)) %>%
  mutate(time = if_else(season == 's16-17' & time > 7, time - 1, time))

fig12s <- ggplot(data = res_all_ages, aes(x = time, y = val, col = age)) +
  geom_line() +
  facet_grid(season ~ virus, scales = 'free') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = seq(1, 52, by = 5),
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_viridis(discrete = TRUE) +
  labs(x = 'Week #', y = '% Positive', color = 'Age')
ggsave('results/plots/figures_for_manuscript/supp/FigureS12.svg', fig12s, width = 9.5, height = 9)

rm(fig12s, res_all_ages, covar_all_ages)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 13: Comparison between observed data and synthetic data generated by the age-structured model

res_combined <- read_csv('data/age_structured_SA/synthetic_obs_combined.csv')
seasons <- unique(res_combined$season)

hk_dat <- NULL
for (yr in seasons) {
  
  hk_dat_temp <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_rsv %>%
    filter(season == yr) %>%
    select(time, season, n_P1:n_P2)
  hk_dat <- bind_rows(hk_dat, hk_dat_temp)
  
}

hk_dat_long <- hk_dat %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'obs') %>%
  mutate(virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'))

res_combined_long <- res_combined %>%
  select(time, season, obs1_s1b, obs2_s1b) %>%
  pivot_longer(obs1_s1b:obs2_s1b,
               names_to = 'virus',
               values_to = 'synth') %>%
  mutate(virus = if_else(virus == 'obs1_s1b', 'Influenza', 'RSV'))

res_combined_long <- res_combined_long %>%
  inner_join(hk_dat_long,
             by = c('time', 'season', 'virus'))

fig13s <- ggplot(data = res_combined_long,
                 aes(x = time, color = virus)) +
  geom_point(aes(y = obs)) +
  geom_line(aes(y = synth)) +
  facet_wrap(~ season, scale = 'free', nrow = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_manual(values = brewer.pal(3, 'Dark2')[c(1, 3)]) +
  labs(x = 'Week #', y = '# of Cases', color = 'Virus')
ggsave('results/plots/figures_for_manuscript/supp/FigureS13.svg', width = 21, height = 4.5, fig13s)

rm(list = ls())

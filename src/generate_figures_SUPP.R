# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript (supplement)
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(lemon)
library(gridExtra)
library(grid)
library(patchwork)
library(cowplot)
library(GGally)
library(RColorBrewer)
library(viridis)
library(purrr)

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
season_breaks <- dat_hk %>% filter(Week == 46) %>% pull(Time)

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
  labs(x = 'Year', y = '\n% Positive', col = '(Sub)type') +
  geom_vline(xintercept = season_breaks, linetype = 'dashed')

fig1s <- reposition_legend(fig1s, position = 'top left', plot = FALSE)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS1.svg', width = 9.5, height = 4, fig1s)

rm(dat_pos, fig1s, season_breaks)

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
# ggsave('results/plots/figures_for_manuscript/supp/FigureS2.svg', width = 10.5, height = 9, fig2s)

rm(dat_hk, dat_clim, p2a, p2b, p2c, fig2s, x_lab_breaks)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 3: Season-specific parameters and CIs

shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

mle <- read_csv('results/MLE_plus_95CI_from_boostrapping_HPDI.csv')
res <- mle %>%
  filter(!(parameter %in% shared_estpars) & parameter != 'delta2') %>%
  mutate(season = str_sub(parameter, 2, 6),
         parameter = str_sub(parameter, 8)) %>%
  mutate(outside = mle < lower | mle > upper)

res <- res %>%
  mutate(parameter = if_else(parameter == 'Ri1', 'Ri[1]', parameter),
         parameter = if_else(parameter == 'Ri2', 'Ri[2]', parameter),
         parameter = if_else(parameter == 'I10', 'I[10]', parameter),
         parameter = if_else(parameter == 'I20', 'I[20]', parameter),
         parameter = if_else(parameter == 'R10', 'R[10]', parameter),
         parameter = if_else(parameter == 'R20', 'R[20]', parameter),
         parameter = if_else(parameter == 'R120', 'R[120]', parameter),
         parameter = if_else(parameter == 'R10 + R120', 'R[10] + R[120]', parameter),
         parameter = if_else(parameter == 'R20 + R120', 'R[20] + R[120]', parameter))

res$parameter <- factor(res$parameter, levels = c('Ri[1]', 'Ri[2]', 'I[10]', 'I[20]', 'R[10]', 'R[20]', 'R[120]', 'R[10] + R[120]', 'R[20] + R[120]'))

fig3s <- ggplot(data = res %>% filter(!is.na(parameter))) +
  geom_errorbar(aes(x = season, ymin = lower, ymax = upper)) +
  geom_point(aes(x = season, y = mle, col = outside)) +
  facet_wrap(~ parameter, scales = 'free_y', ncol = 1, labeller = 'label_parsed') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  scale_color_manual(values = c('black', 'darkred'), guide = 'none') +
  labs(x = 'Season', y = 'Maximum Likelihood Estimate (95% CI)')

# ggsave('results/plots/figures_for_manuscript/supp/FigureS3.svg', width = 4, height = 20.0, fig3s)

rm(fig3s, res, mle)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 4: Relationship between proportion immune and previous season's attack rate

ar_df <- read_rds('results/simulated_ar.rds')
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')
dat_temp <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv

prop_imm_df <- mle[1, ] %>%
  select(contains('R1') | contains('R2')) %>%
  pivot_longer(everything()) %>%
  mutate(season = str_sub(name, 1, 6),
         parameter = str_sub(name, 8)) %>%
  select(-name) %>%
  pivot_wider(names_from = parameter,
              values_from = value) %>%
  mutate(r_flu = R10 + R120, r_rsv = R20 + R120) %>%
  select(season, r_flu:r_rsv) %>%
  pivot_longer(r_flu:r_rsv,
               names_to = 'virus',
               values_to = 'prop_imm') %>%
  mutate(virus = if_else(virus == 'r_flu', 'Influenza', 'RSV'))

tot_cases_df <- dat_temp %>%
  select(n_P1:n_P2, season) %>%
  group_by(season) %>%
  summarise(tot_cases1 = sum(n_P1, na.rm = TRUE),
            tot_cases2 = sum(n_P2, na.rm = TRUE)) %>%
  pivot_longer(tot_cases1:tot_cases2,
               names_to = 'virus',
               values_to = 'tot_cases') %>%
  mutate(virus = if_else(virus == 'tot_cases1', 'Influenza', 'RSV'))

ar_df <- ar_df %>%
  inner_join(tot_cases_df,
             by = c('season', 'virus'))
rm(tot_cases_df)

ar_df <- ar_df %>%
  mutate(season = paste0('s', as.numeric(str_sub(season, 2, 3)) + 1, '-', as.numeric(str_sub(season, 5, 6)) + 1))

prop_imm_df <- prop_imm_df %>%
  inner_join(ar_df,
             by = c('season', 'virus'))

p4a <- ggplot(data = prop_imm_df, aes(x = attack_rate, y = prop_imm)) +
  geom_point(size = 2.5) +
  facet_wrap(~ virus) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.015, 0.95)) +
  labs(x = 'Estimated Total Attack Rate', y = 'Proportion Immune', tag = 'A')
p4b <- ggplot(data = prop_imm_df, aes(x = tot_cases, y = prop_imm)) +
  geom_point(size = 2.5) +
  facet_wrap(~ virus) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.015, 0.95)) +
  labs(x = 'Total Observed Cases', y = 'Proportion Immune', tag = 'B')

fig4s <- arrangeGrob(p4a, p4b, ncol = 1)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS4.svg', width = 9, height = 6.5, fig4s)

rm(fig4s, p4a, p4b, ar_df, prop_imm_df, mle, dat_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 5: Correlations between estimated parameters (shared parameters)

res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'

load_and_format_mega_results <- function(filepath, run_name) {
  
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
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime'))
  pars_top <- pars_top %>%
    select(-message)
  
  # Add run name to tibble:
  pars_top <- pars_top %>%
    mutate(condition = run_name)
  
  # Return formatted results:
  return(pars_top)
  
}

pars_top <- load_and_format_mega_results(filepath = res_dir,
                                         run_name = 'H1_plus_B_FULL')
rm(res_dir)

pars_top_temp <- pars_top
names(pars_top_temp)[1:12] <- c('rho[1]', 'rho[2]', 'theta[lambda*1]', 'theta[lambda*2]', 'delta[1]', 'd[2]', 'alpha', 'phi', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]')
shared_estpars_temp <- c('rho[1]', 'rho[2]', 'theta[lambda*1]', 'theta[lambda*2]', 'delta[1]', 'd[2]', 'alpha', 'phi', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]')

fig5s <- ggpairs(pars_top_temp %>% select(all_of(shared_estpars_temp)),
                 upper = list(continuous = wrap(ggally_cor, size = 4.5, method = 'kendall', digits = 2, display_grid = FALSE)),
                 lower = list(continuous = wrap('points', size = 1.1)),
                 labeller = 'label_parsed') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 55, vjust = 0.6),
        strip.text = element_text(size = 15),
        panel.border = element_rect(linewidth = 0.5, fill = NA))

# ggsave('results/plots/figures_for_manuscript/supp/FigureS5.svg', fig5s, width = 18, height = 12)

rm(fig5s, pars_top_temp, shared_estpars_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 6: Correlations between estimated parameters (season-specific R_eff and % immune)

pars_top_unit <- pars_top %>%
  select(contains('Ri') | contains('R1') | contains('R20')) %>%
  mutate(`s13-14_R10 + R120` = `s13-14_R10` + `s13-14_R120`,
         `s14-15_R10 + R120` = `s14-15_R10` + `s14-15_R120`,
         `s15-16_R10 + R120` = `s15-16_R10` + `s15-16_R120`,
         `s16-17_R10 + R120` = `s16-17_R10` + `s16-17_R120`,
         `s17-18_R10 + R120` = `s17-18_R10` + `s17-18_R120`,
         `s18-19_R10 + R120` = `s18-19_R10` + `s18-19_R120`,
         `s13-14_R20 + R120` = `s13-14_R20` + `s13-14_R120`,
         `s14-15_R20 + R120` = `s14-15_R20` + `s14-15_R120`,
         `s15-16_R20 + R120` = `s15-16_R20` + `s15-16_R120`,
         `s16-17_R20 + R120` = `s16-17_R20` + `s16-17_R120`,
         `s17-18_R20 + R120` = `s17-18_R20` + `s17-18_R120`,
         `s18-19_R20 + R120` = `s18-19_R20` + `s18-19_R120`) %>%
  select(contains('Ri') | contains('+')) %>%
  mutate(fit = 1:nrow(pars_top)) %>%
  pivot_longer(cols = !fit,
               names_to = 'parameter',
               values_to = 'value') %>%
  mutate(season = str_sub(parameter, 1, 6)) %>%
  mutate(parameter = str_sub(parameter, 8))

tags <- c('A', 'B', 'C', 'D', 'E', 'F')
seasons <- unique(pars_top_unit$season)
plot_list <- vector('list', length = length(seasons))

for (seas_index in 1:length(seasons)) {
  
  seas <- seasons[seas_index]
  
  pars_top_unit_temp <- pars_top_unit %>%
    filter(season == seas) %>%
    select(!season) %>%
    pivot_wider(names_from = 'parameter',
                values_from = 'value')
  names(pars_top_unit_temp)[2:5] <- c('Ri[1]', 'Ri[2]', 'R[10] + R[120]', 'R[20] + R[120]')
  
  p_temp <- ggpairs(pars_top_unit_temp %>% select(2, 4, 3, 5),
                    upper = list(continuous = wrap(ggally_cor, size = 4.5, method = 'kendall', digits = 2, display_grid = FALSE)),
                    lower = list(continuous = wrap('points', size = 1.1)),
                    labeller = 'label_parsed') +
    theme_classic() +
    theme(title = element_text(size = 14),
          axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 55, vjust = 0.6),
          strip.text = element_text(size = 12),
          plot.tag = element_text(size = 22),
          plot.tag.position = c(0.005, 0.988),
          panel.border = element_rect(linewidth = 0.5, fill = NA)) +
    labs(title = str_flatten(c('20', str_sub(seas, 2))),
         tag = tags[seas_index])
  
  plot_list[[seas_index]] <- grid.grabExpr(print(p_temp))
  
}
rm(seas_index, seas, pars_top_unit_temp, p_temp, tags)

fig6s <- arrangeGrob(grobs = plot_list, ncol = 2)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS6.svg', fig6s, width = 18, height = 16.5)

rm(fig6s, pars_top_unit, plot_list)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 7: Climate forcing over time

gamma1 <- 7/5
gamma2 <- 7/10

mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')[1, ]

mle <- mle %>%
  select(eta_temp1:eta_ah2,
         contains('R10'),
         contains('R20'),
         contains('R120'),
         contains('Ri')) %>%
  pivot_longer(-c(eta_temp1:eta_ah2),
               names_to = 'parameter',
               values_to = 'mle') %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_sub(parameter, 8)) %>%
  pivot_wider(names_from = parameter,
              values_from = mle) %>%
  select(eta_temp1:eta_ah2, Ri1:Ri2, R10:R120, season)

hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')

force_t <- vector('list', length = length(seasons))
for (yr_index in 1:length(seasons)) {
  
  seas <- seasons[yr_index]
  
  dat_temp <- hk_dat[['h1_plus_b_rsv']] %>%
    filter(season == seas) %>%
    inner_join(dat_clim,
               by = c('Year' = 'year',
                      'Week' = 'week')) %>%
    select(time, temp, ah)
  
  mle_temp <- mle %>%
    filter(season == seas)
  
  force1_temp <- exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp)
  force2_temp <- exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp)
  
  force1_temp <- bind_cols(1:max(dat_temp$time), force1_temp, seas)
  force2_temp <- bind_cols(1:max(dat_temp$time), force2_temp, seas)
  
  names(force1_temp) <- c('time', 'force1', 'season')
  names(force2_temp) <- c('time', 'force2', 'season')
  
  force_temp <- force1_temp %>%
    inner_join(force2_temp,
               by = c('time', 'season')) %>%
    select(time, force1, force2, season)
  
  force_t[[yr_index]] <- force_temp
  
}
rm(yr_index, seas, dat_temp, mle_temp, force1_temp, force2_temp, force_temp, hk_dat, dat_clim)

force_t <- bind_rows(force_t)

min_val <- min(c(min(force_t$force1), min(force_t$force2)))
max_val <- max(c(max(force_t$force1), max(force_t$force2)))

p7a <- ggplot(data = force_t, aes(x = time, y = force1, col = season)) +
  geom_line() +
  theme_classic() +
  theme(title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.97)) +
  scale_x_continuous(breaks = seq(1, 52, by = 5),
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_y_continuous(limits = c(min_val, max_val)) +
  # scale_y_log10(limits = c(min_val, max_val)) +
  scale_color_manual(values = viridis(6)) +
  labs(x = 'Week Number', y = 'Climate Forcing', title = 'Influenza')
p7b <- ggplot(data = force_t, aes(x = time, y = force2, col = season)) +
  geom_line() +
  theme_classic() +
  theme(title = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.97)) +
  scale_x_continuous(breaks = seq(1, 52, by = 5),
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_y_continuous(limits = c(min_val, max_val)) +
  # scale_y_log10(limits = c(min_val, max_val)) +
  scale_color_manual(values = viridis(6)) +
  labs(x = 'Week Number', y = 'Climate Forcing', title = 'RSV')

p_legend <- ggplot(data = force_t, aes(x = time, y = force1, col = season)) +
  geom_line() +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1)) +
  scale_color_viridis(discrete = TRUE) +
  labs(color = 'Season')
p_legend <- ggplotGrob(p_legend)$grobs[[which(sapply(ggplotGrob(p_legend)$grobs, function(x) x$name) == 'guide-box')]]

fig7s <- arrangeGrob(arrangeGrob(p7a, p7b, ncol = 1), p_legend, nrow = 2, heights = c(15, 1))
# ggsave('results/plots/figures_for_manuscript/supp/FigureS7.svg', fig7s, width = 10, height = 6.25)

rm(mle, force_t, p7a, p7b, p_legend, fig7s, min_val, max_val, gamma1, gamma2)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 8: Simulations at the MLE

true_estpars <- c(shared_estpars, unit_estpars)
prof_lik <- FALSE
fit_canada <- FALSE

if (fit_canada) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

source('src/functions/setup_global_likelilhood.R')

set.seed(12075)
sim_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- pars_top %>%
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

sims <- bind_rows(sim_list)

rm(sim_list, po_list, vir1, yr, resp_mod, pars_temp, sim_temp, i,
   dat_pomp, hk_dat, nrow_check, yr_index)

sims <- sims %>%
  pivot_longer(n_P1:n_P2, names_to = 'vir1', values_to = 'sim') %>%
  pivot_longer(obs1:obs2, names_to = 'vir2', values_to = 'obs') %>%
  mutate(vir1 = if_else(vir1 == 'n_P1', 'Influenza', 'RSV')) %>%
  mutate(vir2 = if_else(vir2 == 'obs1', 'Influenza', 'RSV')) %>%
  filter(vir1 == vir2) %>%
  mutate(vir = vir1) %>%
  select(time:.id, season, vir, sim, obs) %>%
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

fig8s <- ggplot(data = sims) +
  geom_line(aes(x = time, y = sim, group = paste(.id, vir), col = vir), lwd = 0.3) +
  geom_point(aes(x = time, y = obs, group = paste(.id, vir), col = vir), size = 0.75) +
  facet_wrap(~ season, scales = 'free', nrow = 2) +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Week #', y = '# of Cases', col = 'Virus')

# ggsave('results/plots/figures_for_manuscript/supp/FigureS8.svg', width = 14.5, height = 7.5, fig8s)

rm(fig8s, sims)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 9: Profile likelihoods of theta_lambda1 and theta_lambda2

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
    bind_cols(map_chr(str_split(res_files, '_'), 11),
              paste0('0.', map_chr(str_split(map_chr(str_split(res_files, '_'), 12), fixed('.')), 2))) %>%
    rename(run = '...14',
           profpar = '...15') %>%
    mutate(run = as.numeric(run),
           profpar = as.numeric(profpar)) %>%
    arrange(profpar, run)
  expect_true(nrow(res_temp) == length(res_files))
  expect_true(all(is.finite(res_temp$loglik)))
  
  res_temp <- res_temp %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # Return formatted results:
  return(res_temp)
  
}

get_maxloglik_and_ci_cutoff <- function(res) {
  
  # Get maximum log-likelihood value:
  maxloglik <- res %>%
    summarise(loglik = max(loglik)) %>%
    pull(loglik)
  
  # Get cutoff for 95% confidence interval:
  ci_cutoff <- maxloglik - 0.5 * qchisq(df = 1, p = 0.95)
  res <- res %>%
    mutate(ci = ci_cutoff)
  
  # Keep only best-fit run for each value of interaction parameter:
  res <- res %>%
    group_by(profpar) %>%
    filter(rank(-loglik) == 1) %>%
    ungroup
  
  # Return updated results:
  return(res)
  
}

res_proflik1_ZOOM <- load_and_format_proflik_results(filepath = 'results/prof_lik/prof_lik_thetalambda1ZOOM/',
                                                     prof_par = 'theta_lambda1',
                                                     shared_estpars = shared_estpars)
res_proflik1_ZOOM <- get_maxloglik_and_ci_cutoff(res_proflik1_ZOOM)

res_proflik2_ZOOM <- load_and_format_proflik_results(filepath = 'results/prof_lik/prof_lik_thetalambda2ZOOM/',
                                                     prof_par = 'theta_lambda2',
                                                     shared_estpars = shared_estpars)
res_proflik2_ZOOM <- get_maxloglik_and_ci_cutoff(res_proflik2_ZOOM)

# res_proflik1_ZOOM %>%
#   filter(loglik > (max(loglik) - 20))
# res_proflik2_ZOOM %>%
#   filter(loglik > (max(loglik) - 20))

res_proflik1_ZOOM <- res_proflik1_ZOOM %>%
  filter(profpar <= 0.011)
res_proflik2_ZOOM <- res_proflik2_ZOOM %>%
  filter(profpar <= 0.011)

p9a <- ggplot(res_proflik1_ZOOM, aes(x = profpar, y = loglik)) +
  geom_point() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.975)) +
  geom_smooth(method = 'loess', span = 0.7, color = 'black') +
  geom_hline(color = 'black', aes(yintercept = ci), linewidth = 1, lty = 2) +
  labs(x = bquote(theta[lambda*1]), y = 'Log-Likelihood', tag = 'A') +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(limits = c(-9685, -9644), n.breaks = 6)

p9b <- ggplot(res_proflik2_ZOOM, aes(x = profpar, y = loglik)) +
  geom_point() + theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.02, 0.975)) +
  geom_smooth(method = 'loess', span = 0.7, color = 'black') +
  geom_hline(color = 'black', aes(yintercept = ci), linewidth = 1, lty = 2) +
  labs(x = bquote(theta[lambda*2]), y = 'Log-Likelihood', tag = 'B') +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(limits = c(-9685, -9644), n.breaks = 6)

fig9s <- arrangeGrob(p9a, p9b, nrow = 1)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS9.svg', width = 14, height = 5, fig9s)

rm(fig9s, p9a, p9b, res_proflik1_ZOOM, res_proflik2_ZOOM)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 10: Simulated attack rates for flu and RSV at the MLE

mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')

if (fit_canada) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

source('src/functions/setup_global_likelilhood.R')
dat_temp <- hk_dat$h1_plus_b_rsv

traj_list = ar_list_seas = vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  yr <- seasons[i]
  resp_mod <- po_list[[i]]
  
  pars_temp <- mle[1, ] %>%
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
  
  traj_list[[i]] <- traj_temp
  
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
    mutate(season = yr)
  
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
    mutate(season = yr)
  
  ar_list_seas[[i]] <- inner_join(ar_tot, ar_obs, by = 'season') %>%
    select(season, H1:H2, obs1:obs2)
  
}
rm(i, yr, resp_mod, pars_temp, param_mat_temp, traj_temp, ar_tot, ar_obs, rho1, rho2, alpha, phi, rho1_w, rho2_w,
   dat_temp, hk_dat, dat_pomp, po_list, nrow_check, yr_index, obj_fun_list, mle)

ar_df <- bind_rows(ar_list_seas)
rm(ar_list_seas)

ar_df <- ar_df %>%
  select(-c(obs1:obs2)) %>%
  pivot_longer(H1:H2, names_to = 'virus', values_to = 'attack_rate') %>%
  mutate(virus = if_else(virus == 'H1', 'Influenza', 'RSV'),
         attack_rate = attack_rate * 100)

fig10s <- ggplot(data = ar_df, aes(x = virus, y = attack_rate, group = virus)) +
  geom_violin(fill = 'gray90') +
  geom_point(aes(col = season), size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_viridis(discrete = TRUE) +
  # guides(colour = guide_legend(nrow = 1)) +
  labs(x = 'Virus', y = 'Attack Rate (%)', col = 'Season')

# ggsave('results/plots/figures_for_manuscript/supp/FigureS10.svg', width = 6, height = 3.5, fig10s)

rm(fig10s, ar_df, unit_estpars, true_estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 11: Simulations at the MLE, with and without interaction effect

res <- bind_rows(traj_list) %>%
  as_tibble()
rm(traj_list)

res <- res %>%
  mutate(season = str_sub(season, 2))

p11a <- ggplot() +
  geom_line(data = res %>%
              filter(.id %in% c(1, 2)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H2, col = .id)) +
  geom_line(data = res %>%
              filter(.id == 1),
            aes(x = time, y = H1), lty = 2) +
  facet_wrap(~ season, scales = 'free_x', nrow = 1) +
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
  scale_color_manual(values = c('#d95f02', '#fdbf6f')) +
  labs(x = 'Week #', y = 'RSV Incidence', color = '', tag = 'A')
p11b <- ggplot() +
  geom_line(data = res %>%
              filter(.id %in% c(1, 3)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H1, col = .id)) +
  geom_line(data = res %>%
              filter(.id == 1),
            aes(x = time, y = H2), lty = 2) +
  facet_wrap(~ season, scales = 'free_x', nrow = 1) +
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
  scale_color_manual(values = c('#1b9e77', '#8dd3c7')) +
  labs(x = 'Week #', y = 'Influenza Incidence', color = '', tag = 'B')

fig11s <- arrangeGrob(p11a, p11b, ncol = 1)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS11.svg', width = 18, height = 7, fig11s)

rm(p11a, p11b, fig11s, res, pars_top, seasons, vir1, vir2, age_structured, sens,
   d2_max, debug_bool, prof_lik, Ri_max1, Ri_max2)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 12: Model schematic for vaccine simulation study
# Not generated in R

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 13: Impact of LAIV by season for all scenarios

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
  filter(vacc_cov <= 0.70) %>%
  mutate(vacc_cov = vacc_cov * 100)

res_metrics <- res_metrics %>%
  filter(!(climate == 'subtrop' & season == 's13-14') &
           !(climate == 'temp' & season == 's18-19'))

upper_bound_ar <- max(res_metrics$ar2_impact)

get_inset <- function(df){
  p <- ggplot(data = df %>% 
                group_by(season),# %>% 
              # slice(1),
              aes(x = time, y = val, col = Virus)) +
    geom_line() +
    # guides(fill=FALSE) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none',
          strip.text = element_blank()) +
    scale_y_continuous(limits = c(0, 4.1)) +
    scale_color_brewer(palette = 'Dark2') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)')
  return(p)
}
# Source: https://clarewest.github.io/blog-posts/ggplotInset.html

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}
# Source: https://clarewest.github.io/blog-posts/ggplotInset.html

p13a <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
               aes(group = season)) +
  geom_tile(aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  facet_wrap(~ season, nrow = 1, scales = 'fixed') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.position = 'none',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(title = expression(paste('Temperate (', theta[lambda[vacc]], ' = ', theta[lambda*1], ')')),
       x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'RR', tag = 'A')

p_legend1 <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'natural'),
                    aes(group = season)) +
  geom_tile(aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  facet_wrap(~ season, nrow = 1, scales = 'fixed') +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.width = unit(2.0, 'cm'),
        legend.key.height = unit(0.7, 'cm'),
        legend.position = 'bottom') +
  scale_fill_distiller(palette = 'RdBu',
                       values = c(0, 1 / upper_bound_ar, 1),
                       limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks = seq(10, 60, by = 10)) +
  labs(fill = 'RR')
p_legend1 <- ggplotGrob(p_legend1)$grobs[[which(sapply(ggplotGrob(p_legend1)$grobs, function(x) x$name) == 'guide-box')]]

res_simA <- res %>%
  filter(climate == 'temp',
         scenario == 'natural',
         season != 's18-19',
         .id == 1,
         vacc_time == '0',
         vacc_cov == '0.5') %>%
  # select(time:.id, season) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV'))

p_legend2 <- ggplot(data = res_simA, aes(x = time, y = val, col = Virus)) +
  geom_line() +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_brewer(palette = 'Dark2')
p_legend2 <- ggplotGrob(p_legend2)$grobs[[which(sapply(ggplotGrob(p_legend2)$grobs, function(x) x$name) == 'guide-box')]]

insetA <- res_simA %>%
  split(f = .$season) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)),
    data = data.frame(season = unique(.$season)),
    ymin = 2, ymax = 30, xmin = 25, xmax = 52)
  )

p13a <- p13a + insetA

p13b <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'natural'),
               aes(group = season)) +
  geom_tile(aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  facet_wrap(~ season, nrow = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
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

res_simB <- res %>%
  filter(climate == 'subtrop',
         scenario == 'natural',
         season != 's13-14',
         .id == 1,
         vacc_time == '0',
         vacc_cov == '0.5') %>%
  # select(time:.id, season) %>%
  pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
  mutate(val = val * 100) %>%
  mutate(Virus = if_else(Virus == 'H1', 'Influenza', 'RSV'))

insetB <- res_simB %>%
  split(f = .$season) %>%
  purrr::map(~annotation_custom2(
    grob = ggplotGrob(get_inset(.)),
    data = data.frame(season = unique(.$season)),
    ymin = 2, ymax = 30, xmin = 25, xmax = 52)
  )

p13b <- p13b + insetB

p13c <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'half'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
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

p13d <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'half'),
               aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
  geom_tile() +
  facet_wrap(~ season, nrow = 1) +
  theme_classic() +
  theme(title = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
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

# fig13s <- p13a / p13b / p13c / p13d + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
fig13s <- arrangeGrob(p13a, p13b, p13c, p13d,
                      arrangeGrob(p_legend1, p_legend2, layout_matrix = rbind(c(NA, 1, 2, NA)), widths = c(1.25, 2, 2, 1.25)),
                      nrow = 5, heights = c(10, 10, 10, 10, 2.5))
# ggsave('results/plots/figures_for_manuscript/supp/FigureS13.svg', fig13s, width = 13, height = 11)

rm(fig13s, p13a, p13b, p13c, p13d, insetA, insetB, res, res_metrics, res_simA, res_simB, upper_bound_ar)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 14: Sensitivity analyses for the vaccine simulation study

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
  filter(vacc_cov <= 0.70) %>%
  mutate(vacc_cov = vacc_cov * 100)

res_metrics <- res_metrics %>%
  filter((climate == 'subtrop' & season == 's13-14') |
           (climate == 'temp' & season == 's18-19'))

upper_bound_ar <- max(res_metrics$ar2_impact)

p14a <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'deltaShort'),
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

p14b <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'deltaShort'),
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

p14c <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'deltaLong'),
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

p14d <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'deltaLong'),
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

p14e <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'effLow'),
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

p14f <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'effLow'),
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

p14g <- ggplot(data = res_metrics %>% filter(climate == 'temp' & scenario == 'effHigh'),
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

p14h <- ggplot(data = res_metrics %>% filter(climate == 'subtrop' & scenario == 'effHigh'),
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

fig14s <- (p14a + p14b) / (p14c + p14d) / (p14e + p14f) / (p14g + p14h) + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
# ggsave('results/plots/figures_for_manuscript/supp/FigureS14.svg', fig14s, width = 10, height = 15)

rm(fig14s, p14a, p14b, p14c, p14d, p14e, p14f, p14g, p14h, res, res_metrics, upper_bound_ar)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 15: Comparison between observed data and synthetic data generated by the age-structured model

res_combined <- read_csv('data/age_structured_SA/synthetic_obs_combined_ORIG.csv')
seasons <- unique(res_combined$season)

hk_dat <- NULL
for (yr in seasons) {
  
  hk_dat_temp <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv %>%
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

fig15s <- ggplot(data = res_combined_long,
                 aes(x = time, color = virus)) +
  geom_point(aes(y = obs)) +
  geom_line(aes(y = synth)) +
  facet_wrap(~ season, scale = 'free', nrow = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = breaks_fxn,
                     labels = c(46, 51, seq(4, 45, by = 5))) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Week #', y = '# of Cases', color = 'Virus')

# ggsave('results/plots/figures_for_manuscript/supp/FigureS15.svg', width = 14.5, height = 7.5, fig15s)

rm(fig15s, res_combined, hk_dat, hk_dat_temp, yr, hk_dat_long, res_combined_long)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 16: Age-structured synthetic data

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
  mutate(obs1 = obs1 / n_samp * 100,
         obs2 = obs2 / n_samp * 100) %>%
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

fig16s <- ggplot(data = res_all_ages, aes(x = time, y = val, col = age)) +
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

# ggsave('results/plots/figures_for_manuscript/supp/FigureS16.svg', fig16s, width = 9.5, height = 10.5)

rm(fig16s, res_all_ages, covar_all_ages)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 17: Circulation of rhinovirus over time, vs. influenza (H1N1 + B) and RSV

dat_hk <- read_csv('data/formatted/dat_hk.csv')

dat_hk <- dat_hk %>%
  filter(Year < 2020 & !(Year == 2019 & Week > 45)) %>%
  select(Time:n_h1, n_b:n_rsv, n_rhino_est_rnd) %>%
  mutate(n_h1 = n_h1 / n_samp * 100,
         n_b = n_b / n_samp * 100,
         n_h1_b = n_h1_b / n_samp * 100,
         n_rsv = n_rsv / n_samp * 100,
         n_rhino = n_rhino_est_rnd / n_samp * 100) %>%
  select(-n_rhino_est_rnd)

dat_pos <- dat_hk %>%
  select(-n_samp) %>%
  pivot_longer(n_h1:n_rhino,
               names_to = 'virus',
               values_to = 'perc_pos') %>%
  mutate(virus = factor(virus, levels = c('n_h1', 'n_b', 'n_h1_b', 'n_rsv', 'n_rhino'))) %>%
  mutate(virus = recode(virus, n_h1 = 'Influenza A(H1N1)', n_b = 'Influenza (B)',
                        n_h1_b = 'Influenza (A(H1N1) + B)    ', n_rsv = 'RSV  ', n_rhino = 'Rhinovirus'))

x_lab_breaks <- dat_hk %>% filter(Week == 1) %>% pull(Time)
season_breaks <- dat_hk %>% filter(Week == 46) %>% pull(Time)

fig17s <- ggplot(data = dat_pos %>%
                   filter(virus %in% c('Influenza (A(H1N1) + B)    ', 'RSV  ', 'Rhinovirus')),
                 aes(x = Time, y = perc_pos, col = virus)) +
  geom_line() + theme_classic() +
  theme(legend.position = 'bottom',
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_x_continuous(breaks = x_lab_breaks, labels = 2014:2019) +
  scale_y_continuous(limits = c(0, 42)) +
  scale_color_brewer(palette = 'Dark2') +
  labs(x = 'Year', y = '\n% Positive', col = 'Virus') +
  geom_vline(xintercept = season_breaks, linetype = 'dashed')

fig17s <- reposition_legend(fig17s, position = 'top left', plot = FALSE)
# ggsave('results/plots/figures_for_manuscript/supp/FigureS17.svg', width = 9.5, height = 4, fig17s)

rm(fig17s, dat_hk, dat_pos, x_lab_breaks, season_breaks)

# ---------------------------------------------------------------------------------------------------------------------

# Supplementary Figure 18: MLE for shared parameters for a variety of sensitivity analyses
# For now only have 95% CIs for some of these - decide whether ranges need to be plotted, or just MLE values

res_dir_main <- 'results/round2_fit/round2_3_fluH1_plus_B/'
res_dir_noAH <- 'results/round2_fit/sens/no_ah/round2_3_fluH1_plus_B/'
res_dir_sinusoidal <- 'results/round2_fit/sens/sinusoidal_forcing/round2_3_fluH1_plus_B/'
res_dir_noint <- 'results/round2_fit/sens/no_int/round2_2_fluH1_plus_B/'
res_dir_noRSVimmune <- 'results/round2_fit/sens/no_rsv_immune/round2_5_fluH1_plus_B/'
res_dir_h3covar <- 'results/round2_fit/sens/h3_covar/round2_3_fluH1_plus_B/'
res_dir_lesscirch3 <- 'results/round2_fit/sens/less_circ_h3/round2_5_fluH1_plus_B/'
res_dir_rhino <- 'results/round2_fit/sens/rhino_covar/round2_3_fluH1_plus_B/'

res_main <- load_and_format_mega_results(res_dir_main, run_name = 'Main')
res_noah <- load_and_format_mega_results(res_dir_noAH, run_name = 'Temperature Only')
res_sinusoidal <- load_and_format_mega_results(res_dir_sinusoidal, run_name = 'Sinusoidal Forcing')
res_noint <- load_and_format_mega_results(res_dir_noint, run_name = 'No Interaction')
res_noRSVimmune <- load_and_format_mega_results(res_dir_noRSVimmune, run_name = 'R[20] + R[120] = 0')
res_h3covar <- load_and_format_mega_results(res_dir_h3covar, run_name = 'H3 as Covariate')
res_lesscirch3 <- load_and_format_mega_results(res_dir_lesscirch3, run_name = 'Low H3 Circulation')
res_rhino <- load_and_format_mega_results(res_dir_rhino, run_name = 'Rhino as Covariate')

# res_main <- read_csv('results/MLE_plus_95CI_from_boostrapping_HPDI.csv')
# res_h3covar <- read_csv('results/round2_fit/sens/h3_covar/MLE_plus_95CI_from_boostrapping_HPDI.csv')
# res_lesscirch3 <- read_csv('results/round2_fit/sens/less_circ_h3/MLE_plus_95CI_from_boostrapping_HPDI.csv')

res <- bind_rows(res_main,
                 res_noah,
                 res_sinusoidal,
                 res_noint,
                 res_noRSVimmune,
                 res_h3covar,
                 res_lesscirch3,
                 res_rhino)

res <- res %>%
  group_by(condition) %>%
  filter(loglik == max(loglik)) %>%
  ungroup() %>%
  select(all_of(shared_estpars), condition) %>%
  mutate(delta2 = d2 * delta1,
         delta1 = 7 / delta1,
         delta2 = 7 / delta2) %>%
  select(-d2) %>%
  pivot_longer(-condition,
               names_to = 'parameter')

res <- res %>%
  mutate(parameter = if_else(parameter == 'theta_lambda1', 'theta[lambda*1]', parameter),
         parameter = if_else(parameter == 'theta_lambda2', 'theta[lambda*2]', parameter),
         parameter = if_else(parameter == 'delta1', '7 / delta[1]', parameter),
         parameter = if_else(parameter == 'delta2', '7 / delta[2]', parameter),
         parameter = if_else(parameter == 'eta_temp1', 'eta[temp*1]', parameter),
         parameter = if_else(parameter == 'eta_temp2', 'eta[temp*2]', parameter),
         parameter = if_else(parameter == 'eta_ah1', 'eta[ah*1]', parameter),
         parameter = if_else(parameter == 'eta_ah2', 'eta[ah*2]', parameter),
         parameter = if_else(parameter == 'rho1', 'rho[1]', parameter),
         parameter = if_else(parameter == 'rho2', 'rho[2]', parameter))

res <- res %>%
  mutate(condition = factor(condition, levels = c('Main', 'H3 as Covariate', 'Low H3 Circulation', 'No Interaction', 'Temperature Only', 'Sinusoidal Forcing', 'R[20] + R[120] = 0', 'Rhino as Covariate')),
         parameter = factor(parameter, levels = c('theta[lambda*1]', 'theta[lambda*2]', '7 / delta[1]', '7 / delta[2]', 'eta[temp*1]', 'eta[temp*2]', 'eta[ah*1]', 'eta[ah*2]', 'rho[1]', 'rho[2]', 'alpha', 'phi')))

xlabels <- levels(res$condition)
xlabels[7] <- expression(R[20] + R[120] == 0)

fig18s <- ggplot(data = res, aes(x = condition, y = value)) +
  geom_point() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 40, vjust = 1, hjust = 1),
        strip.text = element_text(size = 14)) +
  facet_wrap(~ parameter, scales = 'free_y', ncol = 2, labeller = 'label_parsed') +
  scale_x_discrete(labels = xlabels) +
  labs(x = 'Analysis', y = 'Parameter Value')
# ggsave('results/plots/figures_for_manuscript/supp/FigureS18.svg', width = 7.5, height = 12, fig18s)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

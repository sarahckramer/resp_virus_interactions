# ---------------------------------------------------------------------------------------------------------------------
# Prepare and save figures for manuscript and supplement
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(lemon)
library(gridExtra)
library(patchwork)

# Figure 1: Model schematic
# Not generated in R

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

# Figure 3: Plot heatmaps from simulation study
file_list_hk <- list.files('results/vaccine_simulation_study/simulations/run_initial_hk/', pattern = 'sim', full.names = TRUE)
file_list_temp <- list.files('results/vaccine_simulation_study/simulations/run_initial_temperate/', pattern = 'sim', full.names = TRUE)

res_list <- vector('list', length(file_list_hk))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_hk[i]) %>% mutate(season = str_sub(file_list_hk[i], 72, 77))
}
res_hk <- bind_rows(res_list) %>%
  as_tibble() %>%
  mutate(climate = 'subtrop')

res_list <- vector('list', length(file_list_temp))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list_temp[i]) %>% mutate(season = str_sub(file_list_temp[i], 79, 84))
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

res_metrics_AVG <- res_metrics %>%
  group_by(climate, vacc_cov, vacc_time) %>%
  summarise(ar2_impact = median(ar2_impact))

upper_bound_ar <- max(res_metrics_AVG$ar2_impact)

# p3 <- ggplot(data= res_metrics_AVG,
#              aes(x = vacc_time, y = vacc_cov, fill = ar2_impact)) +
#   geom_tile() + facet_wrap(~ climate, ncol = 2) +
#   theme_classic() + theme(strip.text = element_blank()) +
#   scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1),
#                        breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
#   scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
#   labs(x = 'Week of Vaccination', y = 'Vaccine Coverage (%)', fill = 'Impact')
# p3

p3a <- ggplot(data= res_metrics_AVG %>% filter(climate == 'temp'),
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
p3b <- ggplot(data= res_metrics_AVG %>% filter(climate == 'subtrop'),
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

fig3 <- grid_arrange_shared_legend(p3a, p3b, ncol = 2, plot = FALSE)
ggsave('results/plots/figures_for_manuscript/Figure3.svg', fig3, width = 10, height = 5)

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

# ---------------------------------------------------------------------------------------------------------------------
# Explore fit parameter values from single-season fits
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(gridExtra)

# Set directory where results from round1 fits are stored:
res_dir <- 'results/round1_fitsharedFALSE/'

# Check that directory for storing plots exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/plots/')) {
  dir.create('results/plots')
}

# Specify parameters estimated:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')

# Get vectors of flu types/seasons:
flu_types <- c('flu_h1', 'flu_b')
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

# ---------------------------------------------------------------------------------------------------------------------

# Read in and format results

# Read in full results:
pars_list <- read_rds(paste0(res_dir, 'traj_match_round1_byvirseas_FULL.rds'))

# Get best-fit values and values within 95% CI:
pars_top <- vector('list', length = length(pars_list))
names(pars_top) <- names(pars_list)

for (i in 1:length(pars_list)) {
  no_best <- nrow(subset(pars_list[[i]], 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = length(estpars))))
  pars_top[[i]] <- pars_list[[i]][1:no_best, ]
}

# Combine results into data frame:
pars_df <- pars_top %>%
  bind_rows() %>%
  pivot_longer(all_of(estpars), names_to = 'param', values_to = 'val') %>%
  mutate(param = factor(param, levels = estpars))

# Get MLEs and ranges for each season/parameter:
mle_ranges_df <- pars_df %>%
  group_by(virus1, year, param) %>%
  summarise(mle = val[loglik == max(loglik)],
            min = min(val),
            max = max(val))

# ---------------------------------------------------------------------------------------------------------------------

# Plots

# Plot violin plots of all values within 95% CI:
p1 <- ggplot(data = pars_df, aes(x = year, y = val, group = paste(virus1, year), fill = virus1)) + geom_violin() +
  theme_classic() + facet_wrap(~ param, scales = 'free_y', ncol = 3) +
  scale_fill_brewer(palette = 'Set1') + labs(x = 'Season', y = 'Parameter Value')
print(p1)

# Plot parameter values and ranges:
p2 <- ggplot(data = mle_ranges_df, aes(x = year, y = mle, col = param)) + geom_point(size = 2.5) +
  theme_bw() + facet_grid(param ~ virus1, scales = 'free_y') +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.7))
p3 <- ggplot(data = mle_ranges_df, aes(x = year, y = mle, ymin = min, ymax = max, col = param)) + geom_pointrange() +
  theme_bw() + facet_grid(param ~ virus1, scales = 'free_y') +
  theme(legend.position = 'none', axis.text.x = element_text(angle = 45, vjust = 0.7))
grid.arrange(p2, p3, ncol = 2)

# Output plots to file:
pdf('results/plots/param_est_single_seasons_fitsharedFALSE.pdf', width = 9.5, height = 11)
# print(p1)
grid.arrange(p2, p3, ncol = 2)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

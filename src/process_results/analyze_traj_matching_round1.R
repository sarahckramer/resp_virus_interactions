# ---------------------------------------------------------------------------------------------------------------------
# Analyze results from trajectory mapping round 1
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(corrplot)
library(ppcor)
library(gridExtra)

# Specify parameters estimated:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# ---------------------------------------------------------------------------------------------------------------------

# Plot results by flu type/season

# Read in results:
pars_list <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')
slice_list <- read_rds('results/traj_match_round1_byvirseas_SLICE.rds')

# Get vectors of flu types/seasons:
flu_types <- c('flu_A', 'flu_B')
seasons <- 2006:2014

# Create lists to store results:
cor_list = pcor_list = vector('list', length(pars_list))
names(cor_list) = names(pcor_list) = names(pars_list)

# Loop through flu types and seasons:
for (vir1 in flu_types) {
  for (yr in seasons) {
    
    # Get list position:
    vir_seas <- paste(vir1, yr, sep = '_')
    print(vir_seas)
    
    # Check that results exist for this flu/season:
    if (vir_seas %in% names(pars_list)) {
      
      # Get results for specific vir1/yr:
      pars_temp <- pars_list[[vir_seas]]
      
      # Get correlation coefficients:
      cor_mat <- pars_temp %>% dplyr::select(Ri1:loglik) %>% cor(method = 'spearman')
      cor_list[[vir_seas]] <- cor_mat
      
      # Compute partial correlation coefficients:
      pcor_mat <- cor_mat
      for (ix in 1:nrow(pcor_mat)) {
        for (jx in 1:ncol(pcor_mat)) {
          if (ix != jx) {
            pcor_mat[ix, jx] <- pcor.test(x = pars_temp[, ix + 2],
                                          y = pars_temp[, jx + 2],
                                          z = pars_temp[, seq(1, nrow(pcor_mat))[-c(1:2, ix + 2, jx + 2)]],
                                          method = 'spearman')$estimate
          }
        }
      }
      rm(ix, jx)
      pcor_list[[vir_seas]] <- pcor_mat
      
      # Clean up:
      rm(pars_temp, cor_mat, pcor_mat)
    }
    
  }
}

# Clean up:
rm(vir_seas, vir1, yr)

# Output plots to file (for A and B separately):
pdf('results/plots/290721_trajectory_matching_round1_byVirSeas.pdf', width = 15, height = 10)

for (vir1 in flu_types) {
  pars_list_temp <- pars_list[str_detect(names(pars_list), vir1)]
  cor_list_temp <- cor_list[str_detect(names(cor_list), vir1)]
  pcor_list_temp <- pcor_list[str_detect(names(pcor_list), vir1)]
  slice_list_temp <- slice_list[str_detect(names(slice_list), vir1)]
  
  lapply(pars_list_temp, function(ix) {
    pairs(x = ix %>% dplyr::select(all_of(estpars), 'loglik'))
  })
  
  par(mfrow = c(4, ceiling(length(pars_list_temp) / 2)))
  lapply(cor_list_temp, function(ix) {
    corrplot(ix)
  })
  lapply(pcor_list_temp, function(ix) {
    try(corrplot(ix))
  })
  
  for (i in 1:length(slice_list_temp)) {
    par(mfrow = c(4, 2), bty = 'l')
    for (par in estpars) {
      slices_cur <- filter(slice_list_temp[[i]], slice == par)
      plot(slices_cur[[par]], slices_cur$ll, type = 'l',
           xlab = par, ylab = 'Log-Likelihood',
           main = par)
    }
    rm(par, slices_cur)
  }
  rm(i)
  
}

dev.off()

# Clean up:
rm(vir1, cor_list, pcor_list, slice_list,
   pars_list_temp, cor_list_temp, pcor_list_temp, slice_list_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Plot overall results

# Compile all parameter estimates:
pars_df <- bind_rows(pars_list)
rm(pars_list)

# Format results and get MLEs:
pars_df <- pars_df %>%
  pivot_longer(Ri1:R120, names_to = 'param', values_to = 'value') %>%
  group_by(virus1, year, param) %>%
  mutate(param = factor(param)) %>%
  mutate(mle = value[which.max(loglik)]) %>%
  mutate(par_mle = paste0(param, '=', signif(mle, 3))) %>%
  ungroup %>%
  mutate(year = as.numeric(year))
pars_df$param <- factor(pars_df$param, levels = estpars)

# Distinguish between RSV when fit alongside flu_A vs. flu_B:
pars_df <- pars_df %>%
  mutate(virus2 = 'rsv',
         virus_pair = paste(virus1, virus2, sep = '_')) %>%
  dplyr::select(virus1, virus2, virus_pair, year:par_mle)

# Plot best estimates by loglik and virus:
p1 <- ggplot(data = pars_df, aes(x = value, y = loglik, col = virus_pair)) + geom_point() +
  theme_classic() + facet_wrap(~ param, scales = 'free_x') +
  labs(x = 'Parameter Value', y = 'Log-Likelihood', col = 'Virus Pair') +
  scale_color_brewer(palette = 'Set2')

# Plot range of fits by year and virus:
p2 <- ggplot(data = pars_df, aes(x = year, y = value, group = paste(virus_pair, year), fill = virus_pair)) +
  geom_boxplot() + theme_classic() +
  facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Season', y = 'Parameter Value') +
  scale_x_continuous(breaks = 2006:2014,
                     labels = paste0('s', str_pad(c(5:13), width = 2, side = 'left', pad = '0'),
                                     '-', str_pad(c(6:14), width = 2, side = 'left', pad = '0'))) +
  scale_fill_brewer(palette = 'Set2')

# Plot overall fit range across all years, by virus:
p3 <- ggplot(data = pars_df, aes(x = virus_pair, y = value, group = virus_pair)) +
  geom_violin(fill = 'gray92') + geom_point(size = 0.75, position = position_jitter()) +
  theme_classic() + facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Virus', y = 'Parameter Value')

# Plot range of MLEs by virus:
pars_mle <- pars_df %>%
  group_by(virus_pair, year, param) %>%
  filter(loglik == max(loglik))
p4 <- ggplot(data = pars_mle, aes(x = virus_pair, y = mle, group = virus_pair)) +
  geom_violin(fill = 'gray92') + geom_point(size = 1.0, position = position_jitter()) +
  theme_classic() + facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Virus', y = 'Parameter Value', title = 'MLE Values')

# Save plots to file:
pdf('results/plots/290721_trajectory_matching_round1.pdf', width = 15, height = 8)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())
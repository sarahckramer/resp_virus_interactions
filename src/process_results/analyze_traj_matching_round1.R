# ---------------------------------------------------------------------------------------------------------------------
# Analyze results from trajectory mapping round 1
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(corrplot)
library(ppcor)
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

# Get/format date (for saving results):
date <- format(Sys.Date(), '%d%m%y')

# Specify parameters estimated:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')

# Sensitivity analysis?:
sens <- 'main'

# ---------------------------------------------------------------------------------------------------------------------

# Plot results by flu type/season

# Read in results:
pars_list <- read_rds(paste0(res_dir, 'traj_match_round1_byvirseas_TOP.rds'))
slice_list <- read_rds(paste0(res_dir, 'traj_match_round1_byvirseas_SLICE.rds'))

# Get vector of seasons:
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

# Create lists to store results:
cor_list = pcor_list = vector('list', length(pars_list))
names(cor_list) = names(pcor_list) = names(pars_list)

# Loop through flu types and seasons:
for (yr in seasons) {
  
  # Get list position:
  vir_seas <- paste('flu_h1_plus_b', yr, sep = '_')
  print(vir_seas)
  
  # Check that results exist for this flu/season:
  if (vir_seas %in% names(pars_list)) {
    
    # Get results for specific virus/yr:
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

# Remove empty list elements:
cor_list <- cor_list[lapply(cor_list, length) > 0]
pcor_list <- pcor_list[lapply(pcor_list, length) > 0]

# Clean up:
rm(vir_seas, yr)

# Output plots to file:
pdf(paste0('results/plots/', date, '_trajectory_matching_round1_byVirSeas_fitsharedFALSE.pdf'),
    width = 15, height = 10)

# Pairs plots:
lapply(pars_list, function(ix) {
  pairs(x = ix %>% dplyr::select(all_of(estpars), 'loglik'))
})

# Plot correlations:
par(mfrow = c(2, ceiling(length(pars_list) / 2)))
lapply(cor_list, function(ix) {
  corrplot(ix)
})
lapply(pcor_list, function(ix) {
  try(corrplot(ix))
})

# Plot slices:
for (i in 1:length(slice_list)) {
  
  par(mfrow = c(3, 3), bty = 'l')
  for (par in estpars) {
    slices_cur <- filter(slice_list[[i]], slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = par)
  }
  rm(par, slices_cur)
}
rm(i)

# Plot simulations:
try(detach('package:ppcor'))
try(detach('package:MASS'))

plot_list <- vector('list', length = length(pars_list))

for (i in 1:length(pars_list)) {
  yr <- pars_list[[i]] %>% pull(year) %>% unique()
  print(yr)
  vir1 <- 'flu_h1_plus_b'; vir2 <- 'rsv'; debug_bool <- FALSE; Ri_max1 <- 3.0; Ri_max2 <- 3.0; d2_max <- 10.0
  source('src/resp_interaction_model.R')
  
  x0s <- pars_list[[i]] %>% select(all_of(estpars))
  
  coef(resp_mod, estpars) <- x0s[1, ]
  
  # x0_trans <- coef(resp_mod, estpars, transform = TRUE)
  # obj_fun <- traj_objfun(data = resp_mod,
  #                        est = estpars,
  #                        partrans = resp_mod@partrans,
  #                        verbose = TRUE)
  # print(-1 * obj_fun(x0_trans))
  
  sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame') %>%
    select(time:.id, n_P1:n_P2) %>%
    arrange(.id) %>%
    cbind(t(resp_mod@data))
  names(sim_temp)[5:6] <- c('obs1', 'obs2')
  
  p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') + 
    geom_point(aes(x = time, y = obs1, group = .id)) + geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
    theme_classic() +
    labs(x = 'Time', y = '# Positive Tests', title = unique(pars_list[[i]]$year))
  plot_list[[i]] <- p_temp
}

do.call('grid.arrange', plot_list)

dev.off()

# Clean up:
rm(cor_list, pcor_list, slice_list)

# ---------------------------------------------------------------------------------------------------------------------

# Plot overall results

# Compile all parameter estimates:
pars_df <- bind_rows(pars_list)
rm(pars_list)

# Format results and get MLEs:
pars_df <- pars_df %>%
  pivot_longer(all_of(estpars), names_to = 'param', values_to = 'value') %>%
  group_by(year, param) %>%
  mutate(param = factor(param)) %>%
  mutate(mle = value[which.max(loglik)]) %>%
  mutate(par_mle = paste0(param, '=', signif(mle, 3))) %>%
  ungroup()
pars_df$param <- factor(pars_df$param, levels = estpars)

# Plot best estimates by loglik and virus:
p1 <- ggplot(data = pars_df, aes(x = value, y = loglik)) + geom_point() +
  theme_classic() + facet_wrap(~ param, scales = 'free_x') +
  labs(x = 'Parameter Value', y = 'Log-Likelihood')

# Plot range of fits by year and virus:
p2 <- ggplot(data = pars_df, aes(x = year, y = value, group = year)) +
  geom_boxplot(fill = 'gray90') + theme_classic() +
  facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Season', y = 'Parameter Value')

# Plot overall fit range across all years, by virus:
p3 <- ggplot(data = pars_df, aes(x = virus1, y = value)) +
  geom_violin(fill = 'gray92') + geom_point(size = 0.75, position = position_jitter()) +
  theme_classic() + facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Virus', y = 'Parameter Value')

# Plot range of MLEs by virus:
pars_mle <- pars_df %>%
  group_by(year, param) %>%
  filter(loglik == max(loglik))
p4 <- ggplot(data = pars_mle, aes(x = virus1, y = mle)) +
  geom_violin(fill = 'gray92') + geom_point(size = 1.0, position = position_jitter()) +
  theme_classic() + facet_wrap(~ param, scales = 'free_y') +
  labs(x = 'Virus', y = 'Parameter Value', title = 'MLE Values')

# Save plots to file:
pdf(paste0('results/plots/', date, '_trajectory_matching_round1_fitsharedFALSE.pdf'),
    width = 15, height = 8)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

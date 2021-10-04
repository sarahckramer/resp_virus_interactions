# ---------------------------------------------------------------------------------------------------------------------
# Code to generate synthetic "data" and explore the role of different parameters on impact of interactions
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
# library(tidyverse)
# library(pomp)
library(tgp)
# library(gridExtra)

# Set seed:
set.seed(1902574195)

# Set global parameters:
int_parms <- c('theta_lambda1', 'theta_lambda2') # eventually also do theta_rhos and delta
int_eff <- c(seq(0, 1.0, by = 0.1))#, seq(1.25, 2.0, by = 0.25))
n_lhs <- 300
n_sim <- 10
search_type <- 'round1CIs' # 'broad' or 'round1CIs'
realistic_only <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Prep pomp object

# Set necessary parameter values:
vir1 <- 'flu_A' # 'flu_A', 'flu_B'
vir2 <- 'rsv'
yr <- 2006

debug_bool <- FALSE

Ri_max1 <- 2.0
Ri_max2 <- 2.0
delta_min <- 7 / 60.0

# Load pomp object:
source('src/resp_interaction_model.R')

# ---------------------------------------------------------------------------------------------------------------------

# Get parameter values to try

# Set parameter names:
nm_pars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Set parameter ranges:
if (search_type == 'broad') {
  param_bound <- cbind(c(1.0, 1.0, 0, 0, 0, 0, 0),
                       c(Ri_max1, Ri_max2, 1e-3, 1e-3, 0.75, 0.75, 0.75))
} else if (search_type == 'round1CIs') {
  mle_list <- read_rds('results/traj_match_round1_byvirseas_MLE.rds')
  mle_list <- mle_list[-which(str_detect(names(mle_list), '2010'))]
  
  lower_bounds <- mle_list %>% bind_rows() %>% apply(., 2, min)
  upper_bounds <- mle_list %>% bind_rows() %>% apply(., 2, max)
  
  param_bound <- cbind(lower_bounds, upper_bounds)
  
} else {
  stop('Unrecognized search type.')
}

# Draw from ranges:
parms <- lhs(n_lhs, param_bound) %>%
  as.data.frame() %>%
  as_tibble()
names(parms) <- nm_pars

# Correct where R10+R20+R120 > 0.99:
parms <- parms %>%
  mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
  mutate(R10 = if_else(sum >= 1.0, R10 - ((sum - 0.9995) * (R10 / sum)), R10),
         R20 = if_else(sum >= 1.0, R20 - ((sum - 0.9995) * (R20 / sum)), R20),
         R120 = if_else(sum >= 1.0, R120 - ((sum - 0.9995) * (R120 / sum)), R120)) %>%
  # mutate(sum_new = I10 + I20 + R10 + R20 + R120) %>%
  select(-sum)

check_sum <- parms %>%
  mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
  pull(sum) %>%
  max()

while (check_sum > 1.0) {
  parms <- parms %>%
    mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
    mutate(R10 = if_else(sum >= 1.0, R10 - ((sum - 0.9995) * (R10 / sum)), R10),
           R20 = if_else(sum >= 1.0, R20 - ((sum - 0.9995) * (R20 / sum)), R20),
           R120 = if_else(sum >= 1.0, R120 - ((sum - 0.9995) * (R120 / sum)), R120)) %>%
    # mutate(sum_new = I10 + I20 + R10 + R20 + R120) %>%
    select(-sum)
  
  check_sum <- parms %>%
    mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
    pull(sum) %>%
    max()
}

# ggplot(data = a) + geom_line(aes(x = R10, y = R10)) + geom_point(aes(x = R10, y = R10_new, col = (sum > 1.0))) + theme_classic()
# ggplot(data = a) + geom_line(aes(x = R20, y = R20)) + geom_point(aes(x = R20, y = R20_new, col = (sum > 1.0))) + theme_classic()
# ggplot(data = a) + geom_line(aes(x = R120, y = R120)) + geom_point(aes(x = R120, y = R120_new, col = (sum > 1.0))) + theme_classic()

# Replicate:
parms <- map_dfr(seq_len(length(int_eff)), ~parms)
expect_equal(dim(parms)[1], n_lhs * length(int_eff))
# https://stackoverflow.com/questions/8753531/repeat-rows-of-a-data-frame-n-times

# ---------------------------------------------------------------------------------------------------------------------

# Generate synthetic data

# Set values in parameter input matrix:
parms_mat <- parmat(params = coef(resp_mod), nrep = n_lhs * length(int_eff))
for (ix in 1:length(nm_pars)) {
  parms_mat[nm_pars[ix], ] <- as.vector(pull(parms[, ix]))
}
rm(ix)
parms_mat_ORIG <- parms_mat

# Replicate int_eff:
int_eff_true <- rep(int_eff, each = n_lhs)
expect_equal(length(int_eff_true), dim(parms_mat)[2])

# Loop through interaction parameters we want to explore:
for (int_parm in int_parms) {
  
  # Remove packages as necessary:
  try(detach('package:ppcor'))
  try(detach('package:MASS'))
  
  # Set interaction parameters accordingly:
  parms_mat <- parms_mat_ORIG
  parms_mat[int_parm, ] <- int_eff_true
  
  # Run simulations:
  sim_dat <- simulate(object = resp_mod,
                      params = parms_mat,
                      nsim = n_sim,
                      format = 'data.frame')
  
  # Format results:
  sim_dat <- sim_dat %>%
    as_tibble() %>%
    select(time:.id, H1_tot:n_P2) %>%
    mutate(id.parm.orig = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
           id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
    mutate(int_parm = int_parm,
           int_parm_val = int_eff_true[id.parm.orig],
           id.parm = (id.parm.orig - 1) %% n_lhs + 1)
  
  # -------------------------------------------------------------------------------------------------------------------
  
  # Calculate outbreak summary statistics
  
  # Calculate attack rates, peak timings, and differences:
  sim_metrics <- sim_dat %>%
    group_by(.id) %>%
    summarise(ar1 = sum(n_P1), ar2 = sum(n_P2),
              pt1 = which.max(n_P1) + 39, pt2 = which.max(n_P2) + 39,
              pt_diff = pt1 - pt2)
  # pt_diff describes how many weeks flu peak occurs AFTER the RSV peak (negative values mean flu came first)
  
  # Join with identifying info:
  sim_metrics <- sim_metrics %>%
    left_join(sim_dat %>%
                select(.id, id.parm, id.sim, id.parm.orig, int_parm:int_parm_val) %>%
                unique(),
              by = '.id')
  
  # Remove where either or both viruses have no discernible outbreak in 50% or more of simulations:
  parm_sets_to_remove <- sim_metrics %>%
    filter(int_parm_val == 1.0) %>%
    mutate(no_outbreak1 = ar1 < 100,
           no_outbreak2 = ar2 < 100) %>%
    group_by(id.parm, no_outbreak1, no_outbreak2) %>%
    summarise(count_low_ar = n()) %>%
    filter(no_outbreak1 | no_outbreak2) %>%
    group_by(id.parm) %>%
    summarise(total_low = sum(count_low_ar)) %>%
    filter(total_low >= n_sim / 2) %>%
    pull(id.parm)
  
  sim_dat <- sim_dat %>% filter(!(id.parm %in% parm_sets_to_remove))
  sim_metrics <- sim_metrics %>% filter(!(id.parm %in% parm_sets_to_remove))
  
  # Identify parameter sets that can produce realistic outbreaks:
  remaining_ids <- unique(sim_metrics$id.parm)
  print(length(remaining_ids))
  
  if (realistic_only) {
    
    sim_metrics_red <- sim_metrics %>%
      filter(int_parm_val == 1.0) %>%
      filter(pt_diff >= 0 & pt_diff < 14)
    
    check_pt_diff <- unique(sim_metrics_red$id.parm)
    print(length(check_pt_diff) / length(remaining_ids))
    
    sim_metrics_red <- sim_metrics_red %>%
      filter(ar2 < 500 & ar1 < 1500 & ar1 > 200)
    
    check_ar <- unique(sim_metrics_red$id.parm)
    print(length(check_ar) / length(remaining_ids))
    
    sim_metrics_red <- sim_metrics_red %>%
      filter(pt1 >= 53 & pt1 <= 63 &
               pt2 >= 46 & pt2 <= 55)
    
    check_pt <- unique(sim_metrics_red$id.parm)
    print(length(check_pt) / length(remaining_ids))
    
    sim_dat <- sim_dat %>%
      filter(id.parm %in% check_pt)
    sim_metrics <- sim_metrics %>%
      filter(id.parm %in% check_pt)
    
  }
  
  # When starting with round1 CIs, AR metric seems to be most discerning.
  # When using broad ranges, pt_diff also plays a large role.
  
  # Save simulations/metrics:
  write_csv(sim_dat, file = paste0('results/synthetic_data/synth_dat_', search_type, '_realistic', realistic_only, '_', int_parm, '.csv'))
  write_csv(sim_metrics, file = paste0('results/synthetic_data/synth_metrics_', search_type, '_realistic', realistic_only, '_', int_parm, '.csv'))
  
  # Plot data to pdf:
  pdf(paste0('results/plots/synth_data_', search_type, '_realistic', realistic_only, '_', int_parm, '.pdf'),
      width = 15, height = 12)
  unique_ids <- sort(unique(sim_dat$id.parm))
  for (i in 1:ceiling(length(unique_ids) / 10)) {
    sim_dat_temp <- sim_dat %>% filter(id.parm %in% unique_ids[seq(i * 10 - 9, i * 10, by = 1)])

    p1 <- ggplot(data = sim_dat_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'steelblue2') +
      geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') +
      theme_classic() + facet_grid(int_parm_val ~ id.parm) +
      labs(x = 'Time', y = 'Lab-Confirmed Cases')
    print(p1)
  }
  dev.off()
  
  # -------------------------------------------------------------------------------------------------------------------
  
  # Plot parameter influences (univariate)
  
  # Join parameter values to data frames:
  parms$id.parm.orig <- 1:(dim(parms)[1])
  
  sim_dat <- sim_dat %>%
    left_join(parms, by = 'id.parm.orig')
  sim_metrics <- sim_metrics %>%
    left_join(parms, by = 'id.parm.orig')
  
  # Calculate correlations with int_eff value by parameter values:
  sim_corr <- sim_metrics %>%
    # group_by(id.parm, id.sim) %>%
    group_by(id.parm) %>%
    summarise(cor_ar1 = cor(int_parm_val, ar1, method = 'kendall'),
              cor_ar2 = cor(int_parm_val, ar2, method = 'kendall'),
              cor_pt1 = cor(int_parm_val, pt1, method = 'kendall'),
              cor_pt2 = cor(int_parm_val, pt2, method = 'kendall'),
              cor_pt_diff = cor(int_parm_val, pt_diff, method = 'kendall')) %>%
    # left_join(sim_metrics %>%
    #             dplyr::select(id.parm:id.sim, Ri1:R120) %>%#, ar1:pt_diff) %>%
    #             unique())
    left_join(sim_metrics %>%
                dplyr::select(id.parm, Ri1:R120) %>%#, ar1:pt_diff) %>%
                unique())
  # can do overall corr rather than by simulation to reduce # of points
  
  # # Option to remove results where we don't expect any effect:
  # if (int_parm %in% c('theta_lambda1', 'theta_rho1')) {
  #   sim_corr <- sim_corr %>%
  #     select(-c(cor_ar1, cor_pt1))
  # } else if (int_parm %in% c('theta_lambda2', 'theta_rho2')) {
  #   sim_corr <- sim_corr %>%
  #     select(-c(cor_ar2, cor_pt2))
  # }
  # # pairs(sim_corr[, 3:(dim(sim_corr)[2])], pch = 20)
  
  # Lengthen results for plotting:
  sim_corr_long <- sim_corr %>% mutate(R20plusR120 = R20 + R120) %>%
    pivot_longer(cor_ar1:cor_pt_diff, names_to = 'metric', values_to = 'corr') %>%
    pivot_longer(Ri1:R20plusR120, names_to = 'param', values_to = 'param_vals') %>%
    mutate(metric = factor(metric),
           param = factor(param))
  
  sim_corr_long$metric <- factor(sim_corr_long$metric, levels = levels(sim_corr_long$metric)[c(1:2, 4:5, 3)])
  sim_corr_long$param <- factor(sim_corr_long$param, levels = levels(sim_corr_long$param)[c(7:8, 1:3, 5, 4, 6)])
  
  # Plot correlations:
  pdf(paste0('results/plots/param_inf_', search_type, '_realistic', realistic_only, '_', int_parm, '.pdf'),
      width = 14, height = 9)
  p2 <- ggplot(data = sim_corr_long, aes(x = param_vals, y = corr, col = param)) +
    geom_point(size = 0.8) + facet_grid(metric ~ param, scales = 'free_x') +
    geom_hline(yintercept = 0, lty = 2) + theme_classic() + theme(legend.position = 'none') +
    labs(x = 'Parameter Values', y = 'Kendall Rank Corr') +
    scale_y_continuous(limits = c(-1, 1)) + scale_color_brewer(palette = 'Set2')
  print(p2)
  
  # ggplot(data = sim_metrics, aes(x = int_parm_val, y = ar1, col = R120, group = paste(id.parm, id.sim))) +
  #   geom_line() + theme_classic() + scale_color_viridis()
  # # too busy probably
  
  # What about if we limit analysis to where flu/RSV peaks at least one week before the other?
  flu_first <- sim_metrics %>% filter(int_parm_val == 1.0 & pt_diff < 0) %>% pull(id.parm) %>% unique()
  rsv_first <- sim_metrics %>% filter(int_parm_val == 1.0 & pt_diff > 0) %>% pull(id.parm) %>% unique()
  
  sim_corr_long_flu_first <- sim_corr_long %>% filter(id.parm %in% flu_first)
  sim_corr_long_rsv_first <- sim_corr_long %>% filter(id.parm %in% rsv_first)
  
  p3 <- ggplot(data = sim_corr_long_flu_first, aes(x = param_vals, y = corr, col = param)) +
    geom_point(size = 0.8) + facet_grid(metric ~ param, scales = 'free_x') +
    geom_hline(yintercept = 0, lty = 2) + theme_classic() + theme(legend.position = 'none') +
    labs(x = 'Parameter Values', y = 'Kendall Rank Corr', title = 'Flu First Only') +
    scale_y_continuous(limits = c(-1, 1)) + scale_color_brewer(palette = 'Set2')
  p4 <- ggplot(data = sim_corr_long_rsv_first, aes(x = param_vals, y = corr, col = param)) +
    geom_point(size = 0.8) + facet_grid(metric ~ param, scales = 'free_x') +
    geom_hline(yintercept = 0, lty = 2) + theme_classic() + theme(legend.position = 'none') +
    labs(x = 'Parameter Values', y = 'Kendall Rank Corr', title = 'RSV First Only') +
    scale_y_continuous(limits = c(-1, 1)) + scale_color_brewer(palette = 'Set2')
  print(p3)
  print(p4)
  
  # -------------------------------------------------------------------------------------------------------------------
  
  # Plot partial correlations:
  
  # Load libraries:
  library(corrplot)
  library(ppcor)
  
  # Loop through metrics and plot results:
  corr_names <- names(sim_corr)[str_detect(names(sim_corr), 'cor')]
  
  par(mfrow = c(length(corr_names), 2))
  for (metric in corr_names) {
    sim_corr_temp <- sim_corr %>%
      ungroup() %>%
      dplyr::select(Ri1:R120, metric) %>%
      drop_na()
    
    cor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      cor_mat_temp[1, ix] <- cor(x = sim_corr_temp[, metric], y = sim_corr_temp[, ix], method = 'kendall', use = 'pairwise.complete.obs')
    }
    
    rownames(cor_mat_temp) <- metric
    colnames(cor_mat_temp) <- nm_pars
    
    pcor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    p_vals <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      pcor_test_temp <- pcor.test(x = sim_corr_temp[, metric],
                                  y = sim_corr_temp[, ix],
                                  z = sim_corr_temp[, -c(ix, which(names(sim_corr_temp) == metric))],
                                  method = 'kendall')
      pcor_mat_temp[1, ix] <- pcor_test_temp$estimate
      p_vals[1, ix] <- pcor_test_temp$p.value
    }
    
    rownames(pcor_mat_temp) <- metric
    colnames(pcor_mat_temp) <- nm_pars
    
    rownames(p_vals) <- metric
    colnames(p_vals) <- nm_pars
    
    # par(mfrow = c(2, 1))
    corrplot(cor_mat_temp)
    corrplot(pcor_mat_temp, p.mat = p_vals, insig = 'p-val', sig.level = 0.005)
  }
  
  # Again, try separately for outbreaks with flu vs. RSV peaking first:
  par(mfrow = c(length(corr_names), 2))
  for (metric in corr_names) {
    sim_corr_temp <- sim_corr %>%
      filter(id.parm %in% flu_first) %>%
      ungroup() %>%
      dplyr::select(Ri1:R120, metric) %>%
      drop_na()
    
    cor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      cor_mat_temp[1, ix] <- cor(x = sim_corr_temp[, metric], y = sim_corr_temp[, ix], method = 'kendall', use = 'pairwise.complete.obs')
    }
    
    rownames(cor_mat_temp) <- metric
    colnames(cor_mat_temp) <- nm_pars
    
    pcor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    p_vals <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      pcor_test_temp <- pcor.test(x = sim_corr_temp[, metric],
                                  y = sim_corr_temp[, ix],
                                  z = sim_corr_temp[, -c(ix, which(names(sim_corr_temp) == metric))],
                                  method = 'kendall')
      pcor_mat_temp[1, ix] <- pcor_test_temp$estimate
      p_vals[1, ix] <- pcor_test_temp$p.value
    }
    
    rownames(pcor_mat_temp) <- metric
    colnames(pcor_mat_temp) <- nm_pars
    
    rownames(p_vals) <- metric
    colnames(p_vals) <- nm_pars
    
    # par(mfrow = c(2, 1))
    corrplot(cor_mat_temp)
    corrplot(pcor_mat_temp, p.mat = p_vals, insig = 'p-val', sig.level = 0.005)
  }
  
  par(mfrow = c(length(corr_names), 2))
  for (metric in corr_names) {
    sim_corr_temp <- sim_corr %>%
      filter(id.parm %in% rsv_first) %>%
      ungroup() %>%
      dplyr::select(Ri1:R120, metric) %>%
      drop_na()
    
    cor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      cor_mat_temp[1, ix] <- cor(x = sim_corr_temp[, metric], y = sim_corr_temp[, ix], method = 'kendall', use = 'pairwise.complete.obs')
    }
    
    rownames(cor_mat_temp) <- metric
    colnames(cor_mat_temp) <- nm_pars
    
    pcor_mat_temp <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    p_vals <- matrix(NA, nrow = 1, ncol = length(nm_pars))
    for (ix in 1:length(nm_pars)) {
      pcor_test_temp <- pcor.test(x = sim_corr_temp[, metric],
                                  y = sim_corr_temp[, ix],
                                  z = sim_corr_temp[, -c(ix, which(names(sim_corr_temp) == metric))],
                                  method = 'kendall')
      pcor_mat_temp[1, ix] <- pcor_test_temp$estimate
      p_vals[1, ix] <- pcor_test_temp$p.value
    }
    
    rownames(pcor_mat_temp) <- metric
    colnames(pcor_mat_temp) <- nm_pars
    
    rownames(p_vals) <- metric
    colnames(p_vals) <- nm_pars
    
    # par(mfrow = c(2, 1))
    corrplot(cor_mat_temp)
    corrplot(pcor_mat_temp, p.mat = p_vals, insig = 'p-val', sig.level = 0.005)
  }
  dev.off()
  
}

rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

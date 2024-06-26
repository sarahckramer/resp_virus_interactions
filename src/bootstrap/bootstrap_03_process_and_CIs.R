# ---------------------------------------------------------------------------------------------------------------------
# Process results of parametric bootstrapping and output 95% CIs and plots
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(rethinking)
library(gt)

# Get names of all results files:
file_list <- list.files(path = 'results/bootstrapping/flu_H1_plus_B/', full.names = TRUE)

# Are results from a sensitivity analysis?:
if (str_detect(file_list[[1]], 'sinusoidal')) {
  sens <- 'sinusoidal_forcing'
} else if (str_detect(file_list[[1]], 'h3_covar')) {
  sens <- 'h3_covar'
} else if (str_detect(file_list[[1]], 'less_circ_h3')) {
  sens <- 'less_circ_h3'
} else if (str_detect(file_list[[1]], 'no_rsv_immune')) {
  sens <- 'no_rsv_immune'
} else if (str_detect(file_list[[1]], 'no_ah')) {
  sens <- 'no_ah'
} else if (str_detect(file_list[[1]], 'no_int')) {
  sens <- 'no_int'
} else if (str_detect(file_list[[1]], 'rhino_covar')) {
  sens <- 'rhino_covar'
} else if (str_detect(file_list[[1]], 'susc_plus_sev')) {
  sens <- 'susc_plus_sev'
} else {
  sens <- 'main'
}

# Results from Canada?:
if (str_detect(file_list[[1]], 'canada')) {
  fit_canada <- TRUE
  sens <- 'sinusoidal_forcing'
} else {
  fit_canada <- FALSE
}

# Ensure no results missing:
if (str_detect(file_list[[1]], 'PARALLEL')) {
  
  expect_true(length(file_list) == 500)
  
  run_list <- str_split(file_list, '_') %>%
    purrr::map(~ .x[!is.na(as.numeric(.x))]) %>%
    unlist() %>%
    as.numeric()
  print(which(!(1:500 %in% run_list)))
  
} else {
  
  expect_true(length(file_list) == 500 * 10)
  
}

# Read in all results:
res_full = list()
for (i in seq_along(file_list)) {
  res_full[[i]] <- read_rds(file_list[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
if (str_detect(file_list[[1]], 'PARALLEL')) {
  
  res_df <- lapply(res_full, function(ix) {
    lapply(ix, getElement, 'estpars') %>%
      bind_rows() %>%
      bind_cols('loglik' = lapply(ix, getElement, 'll') %>%
                  unlist()) %>%
      # mutate(dataset = str_split(file_list, '_') %>% purrr::map(~ .x[which(!is.na(as.numeric(str_split(file_list, '_')[[1]])))]) %>% unlist()) %>%
      bind_cols('conv' = lapply(ix, getElement, 'message') %>%
                  unlist())
  })
  
  res_df <- lapply(1:length(res_df), function(ix) {
    res_df[[ix]]$dataset <- as.numeric(str_split(file_list[ix], '_')[[1]])[which(!is.na(as.numeric(str_split(file_list[ix], '_')[[1]])))]
    res_df[[ix]]
  })
  
  res_df <- res_df %>%
    bind_rows()
  
  expect_true(nrow(res_df) == length(file_list) * 10)
  
} else {
  
  res_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    mutate(dataset = str_split(file_list, '_') %>% purrr::map(~ .x[which(!is.na(as.numeric(str_split(file_list, '_')[[1]])))]) %>% unlist()) %>%
    bind_cols('conv' = lapply(res_full, getElement, 'message') %>%
                unlist())
  
  expect_true(nrow(res_df) == length(file_list))
  
}

expect_true(all(is.finite(res_df$loglik)))

# Remove estimates that did not converge:
table(res_df$conv)

res_df <- res_df %>%
  filter(str_detect(conv, 'XTOL_REACHED')) %>%
  select(-conv)

# Check whether at least one estimate remains for all synthetic datasets:
expect_true(length(unique(res_df$dataset)) == 500)

# Keep only top fit for each synthetic dataset:
res_df <- res_df %>%
  group_by(dataset) %>%
  filter(loglik == max(loglik)) %>%
  ungroup()
expect_true(nrow(res_df) == 500)

# Are all top fits within a reasonable range of log-likelihood values?
hist(res_df$loglik, breaks = 50)

# Get season-specific and shared parameter names:
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
shared_estpars <- res_df %>% select(!contains(unit_estpars) & !'loglik' & !'dataset') %>% names()
true_estpars <- c(shared_estpars, unit_estpars)

# # Check that no states go below zero for any of the top-fit parameter sets:
# prof_lik <- FALSE
# 
# if (fit_canada) {
#   vir1 <- 'flu'
# } else {
#   vir1 <- 'flu_h1_plus_b'
# }
# 
# source('src/functions/setup_global_likelilhood.R')
# 
# traj_list <- lapply(1:length(seasons), function(ix) {
#   pars_temp <- res_df %>%
#     select(all_of(shared_estpars), contains(seasons[ix]))
#   names(pars_temp) <- true_estpars
# 
#   p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(pars_temp))
#   for (param in names(pars_temp)) {
#     p_mat[param, ] <- pars_temp %>% pull(param)
#   }
# 
#   trajectory(object = po_list[[ix]],
#              params = p_mat,
#              format = 'data.frame') %>%
#     select(!(H1_tot:H2_tot)) %>%
#     pivot_longer(X_SS:H2,
#                  names_to = 'state')
# })
# 
# expect_false(any(lapply(traj_list, function(ix) {
#   any(ix[, 'value'] < 0)
# }) %>%
#   unlist()))

# Calculate composite parameters:
if (fit_canada) {
  
  res_df_unit <- res_df %>%
    select(contains('I') | contains('R') | dataset) %>%
    select(-c(phi, rho1, rho2, phi1, phi2))
  
} else {
  
  res_df_unit <- res_df %>%
    select(contains('I') | contains('R') | dataset) %>%
    select(-c(phi, rho1, rho2))
  
  if (sens == 'sinusoidal_forcing') {
    res_df_unit <- res_df_unit %>%
      select(-c(phi1, phi2))
  }
  
  if (sens == 'rhino_covar') {
    res_df_unit <- res_df_unit %>%
      select(-beta_rhino)
  }
  
  if (sens == 'susc_plus_sev') {
    res_df_unit <- res_df_unit %>%
      select(-c(theta_rho1, theta_rho2))
  }
  
}

res_df_unit <- res_df_unit %>%
  pivot_longer(-c(loglik, dataset)) %>%
  mutate(season = str_sub(name, 1, 6),
         name = str_sub(name, 8)) %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  mutate(`R10 + R120` = R10 + R120,
         `R20 + R120` = R20 + R120,
         R01 = Ri1 / (1 - (I10 + R10 + R120)),
         R02 = Ri2 / (1 - (I20 + R20 + R120))) %>%
  pivot_longer(Ri1:R02,
               names_to = 'parameter',
               values_to = 'value') %>%
  mutate(parameter = paste(season, parameter, sep = '_')) %>%
  select(parameter:value)

res_df_shared <- res_df %>%
  select(all_of(shared_estpars), loglik, dataset)

if (sens != 'no_int') {
  
  res_df_shared <- res_df_shared %>%
    mutate(delta2 = d2 * delta1)
  
}

res_df_shared <- res_df_shared %>%
  pivot_longer(-c(loglik, dataset),
               names_to = 'parameter',
               values_to = 'value') %>%
  select(parameter:value)

if (sens != 'no_int') {
  
  res_df <- bind_rows(res_df_shared, res_df_unit) %>%
    mutate(parameter = factor(parameter, levels = c(shared_estpars, 'delta2', unique(res_df_unit$parameter))))
  
} else {
  
  res_df <- bind_rows(res_df_shared, res_df_unit) %>%
    mutate(parameter = factor(parameter, levels = c(shared_estpars, unique(res_df_unit$parameter))))
  
}

# Calculate 95% confidence intervals for each parameter:
ci_res <- res_df %>%
  group_by(parameter) %>%
  summarise(lower = HPDI(value, p = 0.95)[1],
            upper = HPDI(value, p = 0.95)[2])
# ci_res <- res_df_long %>%
#   group_by(parameter) %>%
#   summarise(lower = quantile(value, p = 0.025),
#             upper = quantile(value, p = 0.975))

# Write results to file:
if (fit_canada) {
  write_csv(ci_res, file = 'results/round2_fit/sens/canada/95CI_from_bootstrapping_HPDI.csv')
} else {
  
  if (sens != 'main') {
    write_csv(ci_res, file = paste0('results/round2_fit/sens/', sens, '/95CI_from_bootstrapping_HPDI.csv'))
  } else {
    write_csv(ci_res, file = 'results/95CI_from_bootstrapping_HPDI.csv')
  }
  
}

# Read in MLEs and add to data frame:
if (fit_canada) {
  
  mle <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')[1, ]
  
  mle_unit <- mle %>%
    select(contains('I') | contains('R')) %>%
    select(-c(phi, rho1, rho2, phi1, phi2))
  
} else {
  
  if (sens != 'main') {
    mle <- read_rds(paste0('results/round2_fit/sens/', sens, '/MLEs_flu_h1_plus_b.rds'))[1, ]
  } else {
    mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')[1, ]
  }
  
  mle_unit <- mle %>%
    select(contains('I') | contains('R')) %>%
    select(-c(phi, rho1, rho2))
  
  if (sens == 'sinusoidal_forcing') {
    mle_unit <- mle_unit %>%
      select(-c(phi1, phi2))
  }
  
  if (sens == 'rhino_covar') {
    mle_unit <- mle_unit %>%
      select(-beta_rhino)
  }
  
  if (sens == 'susc_plus_sev') {
    mle_unit <- mle_unit %>%
      select(-c(theta_rho1, theta_rho2))
  }
  
}

mle_unit <- mle_unit %>%
  pivot_longer(everything()) %>%
  mutate(season = str_sub(name, 1, 6),
         name = str_sub(name, 8)) %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  mutate(`R10 + R120` = R10 + R120,
         `R20 + R120` = R20 + R120,
         R01 = Ri1 / (1 - (I10 + R10 + R120)),
         R02 = Ri2 / (1 - (I20 + R20 + R120))) %>%
  pivot_longer(Ri1:R02,
               names_to = 'parameter',
               values_to = 'mle') %>%
  mutate(parameter = paste(season, parameter, sep = '_')) %>%
  select(parameter:mle)

mle_shared <- mle %>%
  select(all_of(shared_estpars)) 

if (sens != 'no_int') {
  
  mle_shared <- mle_shared %>%
    mutate(delta2 = d2 * delta1)
  
}

mle_shared <- mle_shared %>%
  pivot_longer(everything(),
               names_to = 'parameter',
               values_to = 'mle')

mle <- bind_rows(mle_shared, mle_unit)

ci_res <- ci_res %>%
  left_join(mle, by = 'parameter') %>%
  select(parameter, mle, lower:upper)

# Write results to file:
if (fit_canada) {
  write_csv(ci_res, file = 'results/round2_fit/sens/canada/MLE_plus_95CI_from_bootstrapping_HPDI.csv')
} else {
  
  if (sens != 'main') {
    write_csv(ci_res, file = paste0('results/round2_fit/sens/', sens, '/MLE_plus_95CI_from_bootstrapping_HPDI.csv'))
  } else {
    write_csv(ci_res, file = 'results/MLE_plus_95CI_from_bootstrapping_HPDI.csv')
  }
  
}

# Generate tables of results:
res_table <- ci_res %>%
  gt() %>%
  tab_header(title = 'Best-Fit Parameter Values') %>%
  fmt_number(columns = c(mle, lower, upper), decimals = 3, suffixing = TRUE)
print(res_table)

if (fit_canada) {
  gtsave(res_table, filename = 'results/round2_fit/sens/canada/table_CIs.html')
} else {
  
  if (sens != 'main') {
    gtsave(res_table, filename = paste0('results/round2_fit/sens/', sens, '/table_CIs_h1_plus_b.html'))
  } else {
    gtsave(res_table, filename = 'results/plots/table_CIs_h1_plus_b.html')
  }
  
}

# Check whether MLEs fall within CIs:
ci_res %>% filter(mle <= lower)
ci_res %>% filter(mle >= upper)

# Plot range of fit values for each parameter:
res_df_unit <- res_df_unit %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_sub(parameter, 8)) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'R01', 'R02', 'I10', 'I20', 'R10', 'R20', 'R120', 'R10 + R120', 'R20 + R120')))

p1 <- ggplot(data = res_df %>% filter(parameter %in% c(shared_estpars, 'delta2')),
             aes(x = value, y = after_stat(count)/nrow(res_df))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2 <- ggplot(data = res_df_unit, aes(x = value, y = after_stat(count)/nrow(res_df), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1)
print(p2)

if (fit_canada) {
  pdf('results/plots/plot_params_plus_ci_CANADA.pdf', width = 15, height = 8)
} else {
  
  if (sens != 'main') {
    pdf(paste0('results/plots/plot_params_plus_ci_', sens, '.pdf'), width = 15, height = 8)
  } else {
    pdf('results/plots/plot_params_plus_ci.pdf', width = 15, height = 8)
  }
  
}
print(p1)
print(p2)
dev.off()

# Clean up:
rm(list = ls())

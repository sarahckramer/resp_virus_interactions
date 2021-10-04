# ---------------------------------------------------------------------------------------------------------------------
# Code to compile results from trajectory matching for 2009 pandemic, fixing Ri2/I20
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(corrplot)
library(gridExtra)

# Set parameter values:
vir1 <- 'flu_A'
vir2 <- 'rsv'
yr <- 2010
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

# # Output plots to file:
# pdf('results/plots/240921_trajectory_matching_2010only_fixRi2I20.pdf',
#     width = 15, height = 10)

# ---------------------------------------------------------------------------------------------------------------------

# Compile results (theta_lambda1 only)

# Set estimated parameter names:
estpars <- c('Ri1', 'I10', 'R10', 'R20', 'R120', 'theta_lambda1')

# Get list of results files:
res_files <- list.files(path = 'results/240921_fluA_thetalambda1_2010ONLY/', pattern = 'res_', full.names = TRUE)

# Read in all results:
res_full = list()

for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  mutate(virus1 = 'flu_A',
         year = 2010) %>%
  select(virus1:year, Ri1:loglik)

expect_true(nrow(pars_df) == length(res_files))
expect_true(ncol(pars_df) == (length(estpars) + 3))
expect_true(all(is.finite(pars_df$loglik)))

# ---------------------------------------------------------------------------------------------------------------------

# Compile results (theta_lambda1, delta, rho2)

# Set estimated parameter names:
estpars2 <- c('Ri1', 'I10', 'R10', 'R20', 'R120', 'theta_lambda1', 'delta', 'rho2')

# Get list of results files:
res_files <- list.files(path = 'results/240921_fluA_thetalambda1_2010ONLY_rho2_delta/', pattern = 'res_', full.names = TRUE)

# Read in all results:
res_full = list()

for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df2 <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  mutate(virus1 = 'flu_A',
         year = 2010) %>%
  select(virus1:year, Ri1:loglik)

expect_true(nrow(pars_df2) == length(res_files))
expect_true(ncol(pars_df2) == (length(estpars2) + 3))
expect_true(all(is.finite(pars_df2$loglik)))

# ---------------------------------------------------------------------------------------------------------------------

# Sort results and get best fits

# Sort results:
pars_df <- pars_df %>%
  arrange(desc(loglik))
pars_df2 <- pars_df2 %>%
  arrange(desc(loglik))

# Get only best estimates:
no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = length(estpars))))
no_best <- max(no_best, 50) # get top 50 if less than 50
print(no_best)

pars_df <- pars_df[1:no_best, ]

no_best <- nrow(subset(pars_df2, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = length(estpars2))))
no_best <- max(no_best, 50) # get top 50 if less than 50
print(no_best)

pars_df2 <- pars_df2[1:no_best, ]

# Get MLE:
mle <- setNames(object = as.numeric(pars_df[1, estpars]),
                nm = estpars)
mle2 <- setNames(object = as.numeric(pars_df2[1, estpars2]),
                 nm = estpars2)
# Interestingly, mle2 places RSV immunity near 0%

# Plot:
pairs(pars_df %>% select(all_of(estpars), 'loglik'))
pairs(pars_df2 %>% select(all_of(estpars2), 'loglik'))
# Both fit theta_lambda1 very low; fit quality looks decent in both cases?
# Although fits for common parameters look a little less clear in second case

# Clean up:
rm(res_files, res_full, no_best)

# ---------------------------------------------------------------------------------------------------------------

# Calculate correlations

# Get correlation coefficients:
cor_mat <- pars_df %>% dplyr::select(Ri1:loglik) %>% cor(method = 'spearman')
cor_mat2 <- pars_df2 %>% dplyr::select(Ri1:loglik) %>% cor(method = 'spearman')

# Get partial correlation coefficients:
library(ppcor)
pcor_mat <- cor_mat
pcor_mat2 <- cor_mat2

for (ix in 1:nrow(pcor_mat)) {
  for (jx in 1:ncol(pcor_mat)) {
    if (ix != jx) {
      pcor_mat[ix, jx] <- pcor.test(x = pars_df[, ix + 2],
                                    y = pars_df[, jx + 2],
                                    z = pars_df[, seq(1, nrow(pcor_mat))[-c(1:2, ix + 2, jx + 2)]],
                                    method = 'spearman')$estimate
    }
  }
}
rm(ix, jx)

for (ix in 1:nrow(pcor_mat2)) {
  for (jx in 1:ncol(pcor_mat2)) {
    if (ix != jx) {
      pcor_mat2[ix, jx] <- pcor.test(x = pars_df2[, ix + 2],
                                     y = pars_df2[, jx + 2],
                                     z = pars_df2[, seq(1, nrow(pcor_mat2))[-c(1:2, ix + 2, jx + 2)]],
                                     method = 'spearman')$estimate
    }
  }
}
rm(ix, jx)

try(detach('package:ppcor'))
try(detach('package:MASS'))

# Plot:
par(mfrow = c(2, 2))
corrplot(cor_mat)
corrplot(pcor_mat)
corrplot(cor_mat2)
corrplot(pcor_mat2)
# For the first scenario, only Ri1/I10 have strong correlations (-) after controlling for everything else (a
# weaker negative correlation exists between R20 and R120). In the second, delta and rho2 are highly negatively
# correlated, and delta/R20 are also negatively correlated. Loglik seems to positively rely on delta.

# Clean up:
rm(cor_mat, pcor_mat, cor_mat2, pcor_mat2)

# ---------------------------------------------------------------------------------------------------------------

# Calculate slice likelihoods (theta_lambda1 only)

# Load pomp object:
early_start_val <- TRUE
source('src/resp_interaction_model.R')

# Set Ri2/I20:
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

# # Check that log-likelihood calculations are correct:
# obj_fun <- traj_objfun(data = resp_mod,
# est = estpars,
# partrans = resp_mod@partrans,
# verbose = TRUE)
# 
# coef(resp_mod, estpars) <- unname(mle)
# x0_trans <- coef(resp_mod, estpars, transform = TRUE)
# -obj_fun(par = x0_trans)
# pars_df[1, ]
# 
# obj_fun <- traj_objfun(data = resp_mod,
#                        est = estpars2,
#                        partrans = resp_mod@partrans,
#                        verbose = TRUE)
# 
# coef(resp_mod, estpars2) <- unname(mle2)
# x0_trans <- coef(resp_mod, estpars2, transform = TRUE)
# -obj_fun(par = x0_trans)
# pars_df2[1, ]

# Create objective function for call to nloptr:
obj_fun <- traj_objfun(data = resp_mod,
                       est = estpars,
                       partrans = resp_mod@partrans,
                       verbose = TRUE)

# Take slices:
slices <- slice_design(center = mle,
                       Ri1 = seq(from = 0.9 * mle['Ri1'], to = 1.1 * mle['Ri1'], length.out = 20),
                       # Ri2 = seq(from = 0.9 * mle['Ri2'], to = 1.1 * mle['Ri2'], length.out = 20),
                       I10 = seq(from = 0.9 * mle['I10'], to = 1.1 * mle['I10'], length.out = 20),
                       # I20 = seq(from = 0.9 * mle['I20'], to = 1.1 * mle['I20'], length.out = 20),
                       R10 = seq(from = 0.9 * mle['R10'], to = 1.1 * mle['R10'], length.out = 20),
                       R20 = seq(from = 0.9 * mle['R20'], to = 1.1 * mle['R20'], length.out = 20),
                       R120 = seq(from = 0.9 * mle['R120'], to = 1.1 * mle['R120'], length.out = 20),
                       # theta_lambda1 = seq(from = 0.9 * mle['theta_lambda1'], to = 1.1 * mle['theta_lambda1'], length.out = 20),
                       theta_lambda1 = seq(from = 0, to = 1, length.out = 20)) %>%
  mutate(ll = NA)

# Calculate log likelihoods:
for (i in 1:nrow(slices)) {
  coef(resp_mod, estpars) <- unname(slices[i, estpars])
  x0_trans <- coef(resp_mod, estpars, transform = TRUE)
  slices$ll[i] <- -obj_fun(par = x0_trans)
}
rm(i, x0_trans)

# Check that any NAs are due to initial conditions >= 1.0:
init_sums <- slices %>%
  filter(is.na(ll)) %>%
  mutate(init_sum = I10 + 0.0000111 + R10 + R20 + R120) %>%
  pull(init_sum)
expect_true(all(init_sums > 1))

# Remove NAs:
slices <- slices %>%
  filter(!is.na(ll))

# Plot:
par(mfrow = c(3, 2), bty = 'l')
for (par in estpars) {
  slices_cur <- filter(slices, slice == par)
  plot(slices_cur[[par]], slices_cur$ll, type = 'l',
       xlab = par, ylab = 'Log-Likelihood',
       main = par)
}
# Fits mostly look okay

# Clean up:
rm(slices, slices_cur, dat_full, dat_pomp, init_sums, par)

# ---------------------------------------------------------------------------------------------------------------

# Calculate slice likelihoods (theta_lambda1, delta, rho2)

# Load pomp object:
early_start_val <- TRUE
source('src/resp_interaction_model.R')

# Set Ri2/I20:
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

# Create objective function for call to nloptr:
obj_fun <- traj_objfun(data = resp_mod,
                       est = estpars2,
                       partrans = resp_mod@partrans,
                       verbose = TRUE)

slices <- slice_design(center = mle2,
                       Ri1 = seq(from = 0.9 * mle2['Ri1'], to = 1.1 * mle2['Ri1'], length.out = 20),
                       I10 = seq(from = 0.9 * mle2['I10'], to = 1.1 * mle2['I10'], length.out = 20),
                       R10 = seq(from = 0.9 * mle2['R10'], to = 1.1 * mle2['R10'], length.out = 20),
                       R20 = seq(from = 0.9 * mle2['R20'], to = 1.1 * mle2['R20'], length.out = 20),
                       R120 = seq(from = 0.9 * mle2['R120'], to = 1.1 * mle2['R120'], length.out = 20),
                       # theta_lambda1 = seq(from = 0.9 * mle2['theta_lambda1'], to = 1.1 * mle2['theta_lambda1'], length.out = 20),
                       theta_lambda1 = seq(from = 0, to = 1, length.out = 20),
                       delta = seq(7 / 60, 7 / 1, length.out = 20),
                       # rho2 = seq(from = 0.9 * mle2['rho2'], to = 1.1 * mle2['rho2'], length.out = 20),
                       rho2 = seq(from = 0, to = 1, length.out = 20)) %>%
  mutate(ll = NA)

# Calculate log likelihoods:
for (i in 1:nrow(slices)) {
  coef(resp_mod, estpars2) <- unname(slices[i, estpars2])
  x0_trans <- coef(resp_mod, estpars2, transform = TRUE)
  slices$ll[i] <- -obj_fun(par = x0_trans)
}
rm(i, x0_trans)

# Check that any NAs are due to initial conditions >= 1.0:
init_sums <- slices %>%
  filter(is.na(ll)) %>%
  mutate(init_sum = I10 + 0.0000111 + R10 + R20 + R120) %>%
  pull(init_sum)
expect_true(all(init_sums > 1))

# Remove NAs:
slices <- slices %>%
  filter(!is.na(ll))

par(mfrow = c(4, 2), bty = 'l')
for (par in estpars2) {
  slices_cur <- filter(slices, slice == par)
  plot(slices_cur[[par]], slices_cur$ll, type = 'l',
       xlab = par, ylab = 'Log-Likelihood',
       main = par)
}
# In slices, it appears that theta_lambda1 values near 1.0 are actually best?
# Might need to rerun these? Loglikelihood calculations seem off

# Clean up:
rm(slices, slices_cur, dat_full, dat_pomp, init_sums, par)

# ---------------------------------------------------------------------------------------------------------------

# Check fit quality using log-likelihoods

# Load results when Ri2/I20 are fit:
pars_comp <- read_rds('results/240921_round1_earlystart/traj_match_round1_byvirseas_TOP.rds')
pars_comp <- pars_comp$flu_A_2010

# Load pomp object:
early_start_val <- TRUE
source('src/resp_interaction_model.R')

# Set Ri2/I20:
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

# Create objective function for call to nloptr:
estpars_union <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'theta_lambda1', 'delta', 'rho2')
obj_fun <- traj_objfun(data = resp_mod,
                       est = estpars_union,
                       partrans = resp_mod@partrans,
                       verbose = TRUE)

# Calculate loglik for all three runs using same function:
pars_comp$loglik_new <- NA
source('src/resp_interaction_model.R')
for (i in 1:nrow(pars_comp)) {
  coef(resp_mod, c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')) <- unname(pars_comp[i, 3:9])
  x0_trans <- coef(resp_mod, estpars_union, transform = TRUE)
  pars_comp$loglik_new[i] <- -obj_fun(par = x0_trans)
}
rm(i, x0_trans)

pars_df$loglik_new <- NA
source('src/resp_interaction_model.R')
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)
for (i in 1:nrow(pars_df)) {
  coef(resp_mod, estpars) <- unname(pars_df[i, 3:8])
  x0_trans <- coef(resp_mod, estpars_union, transform = TRUE)
  pars_df$loglik_new[i] <- -obj_fun(par = x0_trans)
}
rm(i, x0_trans)

pars_df2$loglik_new <- NA
source('src/resp_interaction_model.R')
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)
for (i in 1:nrow(pars_df2)) {
  coef(resp_mod, estpars2) <- unname(pars_df2[i, 3:10])
  x0_trans <- coef(resp_mod, estpars_union, transform = TRUE)
  pars_df2$loglik_new[i] <- -obj_fun(par = x0_trans)
}
rm(i, x0_trans)

# Join all results:
pars_comp_ll <- pars_comp %>%
  select(virus1, loglik_new) %>%
  mutate(fit = 'Fit Ri2/I20') %>%
  bind_rows(pars_df %>% select(virus1, loglik_new) %>% mutate(fit = 'Fix Ri2/I20')) %>%
  bind_rows(pars_df2 %>% select(virus1, loglik_new) %>% mutate(fit = 'Fit delta/rho2')) %>%
  mutate(fit = factor(fit))
pars_comp_ll$fit <- factor(pars_comp_ll$fit, levels = levels(pars_comp_ll$fit)[c(2:3, 1)])

p1 <- ggplot(data = pars_comp_ll, aes(x = fit, y = loglik_new)) + geom_violin(fill = 'steelblue2') +
  theme_classic() + labs(x = '', y = 'Log Likelihood')
print(p1)
# p2 <- ggplot(data = pars_comp_ll %>% filter(fit %in% levels(fit)[c(1, 3)]), aes(x = fit, y = loglik_new)) +
#   geom_violin(fill = 'steelblue2') +
#   theme_classic() + labs(x = '', y = 'Log Likelihood')
# print(p2)
# Fixing Ri2/I20 results in significantly lower log likelihood, even when flu blocks RSV completely; allowing
# delta and rho2 to also be fit fixes this to some extent, though the average ll is still ~8 points lower
# (ranges overlap, though)

# Compare values of other parameters:
pars_comp_params <- pars_comp %>%
  select(-c(loglik, year, Ri2, I20)) %>%
  pivot_longer(-virus1, names_to = 'Param', values_to = 'val') %>%
  mutate(fit = 'Fit Ri2/I20') %>%
  bind_rows(pars_df %>%
              select(-c(loglik, year, theta_lambda1)) %>%
              pivot_longer(-virus1, names_to = 'Param', values_to = 'val') %>%
              mutate(fit = 'Fix Ri2/I20')) %>%
  bind_rows(pars_df2 %>%
              select(-c(loglik, year, theta_lambda1, delta, rho2)) %>%
              pivot_longer(-virus1, names_to = 'Param', values_to = 'val') %>%
              mutate(fit = 'Fit delta/rho2')) %>%
  mutate(fit = factor(fit),
         fit = factor(fit, levels = levels(fit)[c(2:3, 1)]))

p3 <- ggplot(data = pars_comp_params, aes(x = fit, y = val, group = fit, fill = fit)) + geom_boxplot() +
  facet_wrap(~ Param, scale = 'free_y') + theme_classic() + labs(x = '', y = 'Param. Value')
print(p3)
# Fixing Ri2/I20 results in higher R120, and lower R10/R20; fitting delta/rho results in a much higher values for R10,
# and much lower values for R20/R120

# ---------------------------------------------------------------------------------------------------------------

# Check quality of simulations using fit values

# For situation where Ri2/I20 are fit:
pars_temp <- pars_comp %>%
  select(Ri1:R120)

source('src/resp_interaction_model.R')

par(mfrow = c(1, 1))
matplot(resp_mod@data %>% t(), pch = 20, type = 'b', lty = 1)
obs_df <- resp_mod@data %>%
  t() %>%
  as_tibble() %>%
  mutate(time = 1:ncol(resp_mod@data)) %>%
  rename('obs1' = 'n_P1',
         'obs2' = 'n_P2')

plot_list <- list()
for (i in 1:9) {
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p_temp <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') +
    theme_classic() + labs(x = 'Time', y = 'Cases')
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  # sim_temp <- sim_temp %>% inner_join(obs_df, by = 'time')
  # p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
  #   geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') +
  #   geom_point(aes(x = time, y = obs1, group = .id), col = 'black') +
  #   geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
  #   theme_classic() + labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p_temp
  # print(sum(traj_temp$H1))
}
do.call('grid.arrange', plot_list)

# For fixing Ri2/I20 and fitting theta_lambda1:
pars_temp <- pars_df %>%
  select(Ri1:theta_lambda1)

source('src/resp_interaction_model.R')
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

plot_list <- list()
for (i in 1:9) {
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p_temp <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') +
    theme_classic() + labs(x = 'Time', y = 'Cases')
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  # sim_temp <- sim_temp %>% inner_join(obs_df, by = 'time')
  # p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
  #   geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') +
  #   geom_point(aes(x = time, y = obs1, group = .id), col = 'black') +
  #   geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
  #   theme_classic() + labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p_temp
  # print(sum(traj_temp$H1))
}
do.call('grid.arrange', plot_list)

# For fitting theta_lambda1/delta/rho2:
pars_temp <- pars_df2 %>%
  select(Ri1:rho2)

source('src/resp_interaction_model.R')
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

plot_list <- list()
for (i in 1:9) {
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p_temp <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') +
    theme_classic() + labs(x = 'Time', y = 'Cases')
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  # sim_temp <- sim_temp %>% inner_join(obs_df, by = 'time')
  # p_temp <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
  #   geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') +
  #   geom_point(aes(x = time, y = obs1, group = .id), col = 'black') +
  #   geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
  #   theme_classic() + labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p_temp
  # print(sum(traj_temp$H1))
}
do.call('grid.arrange', plot_list)

# ---------------------------------------------------------------------------------------------------------------

# Clean up
# dev.off()
rm(list = ls())

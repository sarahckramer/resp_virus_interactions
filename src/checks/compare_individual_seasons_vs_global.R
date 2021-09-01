# ---------------------------------------------------------------------------------------------------------------------
# Compare estimates of season-specific parameters, as well as global parameters, under different fitting schemes
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Save as pdf:
pdf('results/plots/010921_trajectory_matching_round2.pdf',
    width = 18, height = 10)

# Load libraries:
library(tidyverse)
library(testthat)

# Set virus 1:
vir1 <- 'flu_B'

# ---------------------------------------------------------------------------------------------------------------------

# Get results from fitting individual seasons

# Read in results:
res_ind <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

# Compile to data frame:
res_ind <- bind_rows(res_ind)

# Keep flu_B only:
res_ind <- res_ind %>%
  filter(virus1 == 'flu_B')

# ---------------------------------------------------------------------------------------------------------------------

# Compile results using global likelihood

# Get list of results files:
res_files <- list.files(path = 'results/300721_CHECK_fluB_noGlobalParams/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# ---------------------------------------------------------------------------------------------------------------------

# Compare estimates

# Reformat for comparison:
res_glob <- pars_top %>%
  pivot_longer(-loglik, names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global')

# Combine data frames:
res_df <- res_ind %>%
  select(year:R120) %>%
  pivot_longer(-year, names_to = 'param') %>%
  mutate(method = 'Seasonal') %>%
  bind_rows(res_glob)

# Plot results:
p1 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p1)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/040821_fluB_1param2/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_broad')
# certainly better than with theta_lambda1, but still having some issues
# theta_lambda2 tends to be lower as Ri1 gets larger

# # Slice over theta_lambda2/delta for top 4 parameter sets:
# estpars <- c('theta_lambda2', 'delta', names(pars_top)[2:36])
# shared_estpars <- c('theta_lambda2', 'delta')
# true_estpars <- c(shared_estpars, unit_estpars)
# source('src/functions/setup_global_likelilhood.R')
# 
# for (i in 1:4) {
#   mle <- setNames(object = c(as.numeric(pars_top[i, 1]),
#                              7 / 5,
#                              as.numeric(pars_top[i, 2:36])),
#                   nm = estpars)
#   slices <- slice_design(center = mle, 
#                          theta_lambda2 = c(seq(from = 0, to = 1.0, by = 0.05),
#                                            seq(from = 1.0, to = 10, by = 0.5)),
#                          delta = 7 / seq(from = 30, to = 1, by = -1)) %>%
#     mutate(ll = NA)
#   
#   for (j in 1:nrow(slices)) {
#     x0 <- slices[j, 1:37]
#     x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
#     slices$ll[j] <- -1 * calculate_global_loglik(x0_trans)
#   }
#   rm(j, x0, x0_trans)
#   
#   par(mfrow = c(1, 2), bty = 'l')
#   for (par in shared_estpars) {
#     slices_cur <- filter(slices, slice == par)
#     plot(slices_cur[[par]], slices_cur$ll, type = 'l',
#          xlab = par, ylab = 'Log-Likelihood',
#          main = par)
#   }
#   rm(par, slices_cur)
# }
# # these are definitely the best estimates for these specific sets of season-specific conditions, but it seems like
# # small differences in these season-specific values can lead to large differences in the best-fitting theta_lambda2

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# lower theta_lambda2 when Ri1 is higher; no clear relationship with I10
# Ri2/I20 not super sensitive to changes; R10 consistently very low (almost 0)
# appears positively related to R20 and negatively to R120, but we know that these two have a tradeoff of their own
# basically, mostly inflexible to parameters describing RSV alone

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using round 1 CIs as starting values

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/050821_fluB_thetalambda2_round1ci/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_CI')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_CI')
# logliks are fairly similar for the full range of theta_lambda2, from 0 to ~2.5
# higher values tend to be fit for lower values of Ri1, but these are very small changes in Ri1
# in one season, also higher values as I0 gets higher - 2011, where outbreaks are simultaneous

# # Slice over theta_lambda2/delta for top 20 parameter sets:
# estpars <- c('theta_lambda2', 'delta', names(pars_top)[2:36])
# shared_estpars <- c('theta_lambda2', 'delta')
# true_estpars <- c(shared_estpars, unit_estpars)
# source('src/functions/setup_global_likelilhood.R')
# 
# for (i in 1:20) {
#   mle <- setNames(object = c(as.numeric(pars_top[i, 1]),
#                              7 / 5,
#                              as.numeric(pars_top[i, 2:36])),
#                   nm = estpars)
#   slices <- slice_design(center = mle,
#                          theta_lambda2 = c(seq(from = 0, to = 1.0, by = 0.05),
#                                            seq(from = 1.0, to = 2, by = 0.5)),
#                          delta = 7 / seq(from = 30, to = 1, by = -1)) %>%
#     mutate(ll = NA)
# 
#   for (j in 1:nrow(slices)) {
#     x0 <- slices[j, 1:37]
#     x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
#     slices$ll[j] <- -1 * calculate_global_loglik(x0_trans)
#   }
#   rm(j, x0, x0_trans)
# 
#   par(mfrow = c(1, 2), bty = 'l')
#   for (par in shared_estpars) {
#     slices_cur <- filter(slices, slice == par)
#     plot(slices_cur[[par]], slices_cur$ll, type = 'l',
#          xlab = par, ylab = 'Log-Likelihood',
#          main = par)
#   }
#   rm(par, slices_cur)
# }
# # again, these are the best estimates for these specific sets of season-specific conditions, but it seems like
# # small differences in these season-specific values can lead to large differences in the best-fitting theta_lambda2

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using broad starting values AND rho2 = 0.15

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/190821_fluB_thetalambda2_rho015_broad/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_broad_rho015')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_broad_rho=0.15')
# still a wide range of theta_lambda2 values; associated with Ri1; highest ll to very low values

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using round 1 CIs as starting values AND rho2 = 0.15

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/160821_fluB_thetalambda2_rho015/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_CI_rho015')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_CI_rho=0.15')
# theta_lambda2 values are narrower, but still quite broad, and associated with Ri1

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using broad starting values AND rho2 = 0.15 AND holding R20=R120=0

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/300821_fluB_thetalambda2_R20_R120_0_rho015_broad/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_broad_rho015_R20R120=0')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_broad_rho=0.15_R20=R120=0')
# still a wide range of theta_lambda2 values; associated with Ri1; highest ll to very low values; lower ll than when R20/R120 allowed to move

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using broad starting values AND rho2 = 0.15 AND holding R10/R20 to realistic values

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/310821_fluB_thetalambda2_R_realistic/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_broad_rho015_R_realistic')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_broad_rho=0.15_R_realistic')
# now theta_lambda2 tends to fit toward 1, but many seem to just get stuck at either 0 or 1; lower ll than when R allowed to be fit
# pretty highly correlated with all other parameters, but params also fit to a really narrow range - still suggests that a very small change in other parameters changes the best-fit theta_lambda2 value...

# How do season-specific values change with different theta_lambda2?:
pars_comp <- pars_top %>%
  select(theta_lambda2, contains('Ri1'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('Ri2'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('I10'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('I20'))
plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda1 using broad starting values AND rho2 = 0.15 AND holding R10/R20 to realistic values

# Set estpars:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/010921_fluB_thetalambda1_R_realistic/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_1_broad_rho015_R_realistic')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Plot results:
res_df$method <- factor(res_df$method)
res_df$method <- factor(res_df$method, levels = levels(res_df$method)[c(9, 1, 3:4, 7:8, 6, 5, 2)])
p2 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p2)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda1_broad_rho=0.15_R_realistic')
# theta_lambda1 tends to fit toward 0, but most seem to just get stuck at either 0; lower ll than when R allowed to be fit

# # How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('I20'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

# For flu_A

# Set virus 1:
vir1 <- 'flu_A'

# ---------------------------------------------------------------------------------------------------------------------

# Get results from fitting individual seasons

# Read in results:
res_ind <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

# Compile to data frame:
res_ind <- bind_rows(res_ind)

# Keep flu_A only:
res_ind <- res_ind %>%
  filter(virus1 == 'flu_A')

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using round 1 CIs as starting values

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/060821_fluA_thetalambda2_round1ci/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_CI')

# Get combined data frame:
res_df <- res_ind %>%
  select(year:R120) %>%
  pivot_longer(-year, names_to = 'param') %>%
  mutate(method = 'Seasonal') %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_CI')
# Almost all estimates are less than 1.0; best seem to concentrate around 0.2, although all within ~3 ll points
# Negative relationship with Ri1, as with flu B; positive with I10 (mostly in 2011/2014, but in general relationship is there)

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using rho2=0.15 and round 1 CIs as starting values

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/170821_fluA_thetalambda2_rho015/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_CI_rho15')

# Get combined data frame:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_CI_rho=0.15')
# Almost all estimates are less than 1.0; tend to cluster around 0.07-0.15, although ll's very similar
# Again, negative relationship with Ri1

# How do season-specific values change with different theta_lambda2?:
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri1'))
# plot(pars_comp, pch = 20)
# almost looks like there's a line that represents an upper bound on theta_lambda2 for different Ri1 values?

# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('Ri2'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('I20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda2, contains('R120'))
# plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda2 using rho2=0.15 and broad starting values and R realistic

# Set estpars:
shared_estpars <- c('theta_lambda2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/310821_fluA_thetalambda2_R_realistic/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_2_broad_rho15_R_realistic')

# Get combined data frame:
res_df <- res_df %>%
  bind_rows(res_glob)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda2, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda2_broad_rho=0.15_R_realistic')
# ll maximized at 1, but values fit often to 1.0 or 0.0, and ll is much lower than when R fit
# correlated with Ri1; with other params, too, but not as much

# How do season-specific values change with different theta_lambda2?:
pars_comp <- pars_top %>%
  select(theta_lambda2, contains('Ri1'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('Ri2'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('I10'))
plot(pars_comp, pch = 20)

pars_comp <- pars_top %>%
  select(theta_lambda2, contains('I20'))
plot(pars_comp, pch = 20)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda1 using rho2=0.15 and broad starting values and R realistic

# Set estpars:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/010921_fluA_thetalambda1_R_realistic/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = (dim(pars_df)[2] - 1))))
no_best <- max(no_best, 50)
print(no_best)

pars_top <- pars_df[1:no_best, ]

# Clean up:
rm(res_files, res_full, no_best)

# Compare:
res_glob <- pars_top %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(param:year) %>%
  mutate(method = 'Global_p1_1_broad_rho15_R_realistic')

# Get combined data frame:
res_df <- res_df %>%
  bind_rows(res_glob)

# Plot results:
res_df$method <- factor(res_df$method)
res_df$method <- factor(res_df$method, levels = levels(res_df$method)[c(5, 3:4, 2, 1)])
p3 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p3)

# Explore correlations:
pars_corr <- pars_top %>% #filter(loglik > -720) %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)
pairs(pars_corr, pch = 20, main = 'theta_lambda1_broad_rho=0.15_R_realistic')
# almost all fit to 0.0, and ll is much lower than when R fit

# ---------------------------------------------------------------------------------------------------------------------

# Close plots:
dev.off()

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Code to confirm that using global likelihood with no global parameters can achieve the same results as fitting
# each season separately (flu_B only)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

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
  # pivot_wider(names_from = param) %>%
  # select(year:R120)

# Combine data frames:
res_df <- res_ind %>%
  select(year:R120) %>%
  pivot_longer(-year, names_to = 'param') %>%
  mutate(method = 'Seasonal') %>%
  bind_rows(res_glob)
  # inner_join(res_glob, by = 'year')

# Plot results:
p1 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p1)

# ---------------------------------------------------------------------------------------------------------------------

# Compile estimates of theta_lambda1

# Set estpars:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars <- c(shared_estpars, unit_estpars)

# Get list of results files:
res_files <- list.files(path = 'results/030821_fluB_1param/', full.names = TRUE)

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
  mutate(method = 'Global_p1')

# Combine data frames:
res_df <- res_df %>%
  bind_rows(res_glob)

# Plot results:
p2 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p2)
# pretty much the same values for all of these parameters regardless of theta_lambda1 value

# Explore correlations:
pars_corr <- pars_top %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  select(-year)

pairs(pars_corr, pch = 20)

# Slice likelihood over theta_lambda1:
estpars <- names(pars_top)[1:36]
mle <- setNames(object = as.numeric(pars_top[1, 1:36]),
                nm = names(pars_top)[1:36])
slices <- slice_design(center = mle,
                       theta_lambda1 = c(seq(from = 0, to = 1.0, by = 0.05),
                                         seq(from = 1.0, to = 20, by = 1.0))) %>%
  mutate(ll = NA)

source('src/setup_global_likelilhood.R')

for (i in 1:nrow(slices)) {
  x0 <- slices[i, 1:36]
  x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
  slices$ll[i] <- -1 * calculate_global_loglik(x0_trans)
}
rm(i, x0, x0_trans)

plot(slices$theta_lambda1, slices$ll, type = 'l', xlab = 'theta_lambda1', ylab = 'Log-Likelihood')
# fits to wildly different values depending on very slight changes in the values for other parameters









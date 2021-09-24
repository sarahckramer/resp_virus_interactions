# ---------------------------------------------------------------------------------------------------------------------
# Code to compile results from trajectory matching for 2009 pandemic, fixing Ri2/I20
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(corrplot)
library(gridExtra)

# Set estimated parameter names:
estpars <- c('Ri1', 'I10', 'R10', 'R20', 'R120', 'theta_lambda1')

# Set parameter values:
vir1 <- 'flu_A'
vir2 <- 'rsv'
yr <- 2010
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

# Output plots to file:
pdf('results/plots/240921_trajectory_matching_2010only_fixRi2I20.pdf',
    width = 15, height = 10)

# ---------------------------------------------------------------------------------------------------------------------

# Compile results

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

# Sort results and get best fits

# Sort results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

# Get only best estimates:
no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.99, df = length(estpars))))
no_best <- max(no_best, 50) # get top 50 if less than 50
print(no_best)

pars_df <- pars_df[1:no_best, ]

# Get MLE:
mle <- setNames(object = as.numeric(pars_df[1, estpars]),
                nm = estpars)

# Plot:
pairs(pars_df %>% select(all_of(estpars), 'loglik'))

# Clean up:
rm(res_files, res_full, no_best)

# ---------------------------------------------------------------------------------------------------------------

# Calculate correlations

# Get correlation coefficients:
cor_mat <- pars_df %>% dplyr::select(Ri1:loglik) %>% cor(method = 'spearman')

# Get partial correlation coefficients:
library(ppcor)
pcor_mat <- cor_mat

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

try(detach('package:ppcor'))
try(detach('package:MASS'))

# Plot:
corrplot(cor_mat)
corrplot(pcor_mat)

# Clean up:
rm(cor_mat, pcor_mat)

# ---------------------------------------------------------------------------------------------------------------

# Calculate slice likelihoods

# Load pomp object:
source('src/resp_interaction_model.R')

# Set Ri2/I20:
coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111)

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

# Clean up:
rm(slices, slices_cur, dat_full, dat_pomp, init_sums, par)

# ---------------------------------------------------------------------------------------------------------------

# Check fit quality

# Compare to lls obtained when Ri2/I20 are fit:
pars_comp <- read_rds('results/240921_round1_earlystart/traj_match_round1_byvirseas_TOP.rds')
pars_comp <- pars_comp$flu_A_2010

pars_comp_ll <- pars_comp %>%
  select(virus1, loglik) %>%
  mutate(fit = 'Fit Ri2/I20') %>%
  bind_rows(pars_df %>% select(virus1, loglik) %>% mutate(fit = 'Fix Ri2/I20'))

p1 <- ggplot(data = pars_comp_ll, aes(x = fit, y = loglik)) + geom_boxplot(fill = 'steelblue2') +
  theme_classic() + labs(x = '', y = 'Log Likelihood')
print(p1)
# Fixing Ri2/I20 results in significantly lower log likelihood, even when flu blocks RSV completely

# Compare values of other parameters:
pars_comp_params <- pars_comp %>%
  select(-c(loglik, year, Ri2, I20)) %>%
  pivot_longer(-virus1, names_to = 'Param', values_to = 'val') %>%
  mutate(fit = 'Fit Ri2/I20') %>%
  bind_rows(pars_df %>%
              select(-c(loglik, year, theta_lambda1)) %>%
              pivot_longer(-virus1, names_to = 'Param', values_to = 'val') %>%
              mutate(fit = 'Fix Ri2/I20'))

p2 <- ggplot(data = pars_comp_params, aes(x = fit, y = val, group = fit, fill = fit)) + geom_boxplot() +
  facet_wrap(~ Param, scale = 'free_y') + theme_classic() + labs(x = '', y = 'Param. Value')
print(p2)
# Fixing Ri2/I20 results in higher R120, and lower R10/R20

# Check quality of simulations using fit values:
pars_temp <- pars_df %>%
  select(Ri1:theta_lambda1)

par(mfrow = c(1, 1))
matplot(resp_mod@data %>% t(), pch = 20, type = 'b', lty = 1)

plot_list <- list()
for (i in 1:56) {
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p3 <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') + theme_classic() +
    labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p3
  # print(sum(traj_temp$H1))
}
do.call('grid.arrange', plot_list)

# ---------------------------------------------------------------------------------------------------------------

# Clean up

dev.off()
rm(list = ls())

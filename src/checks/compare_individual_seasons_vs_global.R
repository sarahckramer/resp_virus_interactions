# ---------------------------------------------------------------------------------------------------------------------
# Compare estimates of season-specific parameters, as well as global parameters, under different fitting schemes
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# # Save as pdf:
# pdf('results/plots/170921_trajectory_matching_round2.pdf',
#     width = 18, height = 10)

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Set virus 1:
vir1 <- 'flu_A'

# ---------------------------------------------------------------------------------------------------------------------

# Function to read in and format results:
load_and_format_mega_results <- function(filepath, shared_estpars, unit_estpars, res_comp, run_name) {
  
  # Compile estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  
  # Read in results:
  res_full = list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
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
  
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = (dim(pars_df)[2] - 1))))
  no_best <- max(no_best, 50)
  
  pars_top <- pars_df[1:no_best, ]
  
  # Compare:
  res_glob <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 4),
           param = str_sub(param, 6, str_length(param))) %>%
    select(param:year) %>%
    mutate(method = run_name)
  
  # Get combined data frame:
  res_df <- res_comp %>%
    bind_rows(res_glob)
  
  # Explore correlations:
  pars_corr <- pars_top %>%
    pivot_longer(-c(all_of(shared_estpars), loglik), names_to = 'param') %>%
    mutate(year = str_sub(param, 1, 4),
           param = str_sub(param, 6, str_length(param))) %>%
    pivot_wider(names_from = param, values_from = value) %>%
    select(-year)
  
  # Return formatted results:
  return(list(pars_top, res_df, pars_corr))
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Get results from fitting individual seasons

# Read in results:
res_ind <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

# Compile to data frame:
res_ind <- bind_rows(res_ind)

# Keep flu_A only:
res_ind <- res_ind %>%
  filter(virus1 == 'flu_A') %>%
  select(year:R120) %>%
  pivot_longer(-year, names_to = 'param') %>%
  mutate(method = 'Seasonal')

# ---------------------------------------------------------------------------------------------------------------------

# Compile results of trajectory matching for theta_lambda1

# Set shared and unit parameters:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Read in and format results:
res_temp <- load_and_format_mega_results(filepath = 'results/160921_fluA_thetalambda1_round1ci_fitR_rho2=15/',
                                         shared_estpars = shared_estpars,
                                         unit_estpars = unit_estpars,
                                         res_comp = res_ind,
                                         run_name = 'Global_CI')
pars_top <- res_temp[[1]]
res_df <- res_temp[[2]]
pars_corr <- res_temp[[3]]

# Pairs plot:
pairs(pars_corr, pch = 20, main = 'theta_lambda1_CI')
# Best estimates are at 1.0, but estimates at 0.0 are within 2-3 ll points

# # Zoom in on theta_lambda1 and likelihood:
# ggplot(data = pars_top, aes(x = theta_lambda1, y = loglik)) + geom_point() + theme_classic()

# Slice over theta_lambda1 for top 5 parameter sets:
estpars <- names(pars_top)[1:64]
true_estpars <- c(shared_estpars, unit_estpars)
source('src/functions/setup_global_likelilhood.R')

par(mfrow = c(2, 3), bty = 'l')
for (i in 1:6) {
  mle <- setNames(object = as.numeric(pars_top[i, 1:64]),
                  nm = estpars)
  slices <- slice_design(center = mle,
                         theta_lambda1 = seq(from = 0, to = 1.0, by = 0.05)) %>%
    mutate(ll = NA)
  
  for (j in 1:nrow(slices)) {
    x0 <- slices[j, 1:66]
    x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
    slices$ll[j] <- -1 * calculate_global_loglik(x0_trans)
  }
  rm(j, x0, x0_trans)
  
  # par(mfrow = c(1, 2), bty = 'l')
  for (par in shared_estpars) {
    slices_cur <- filter(slices, slice == par)
    plot(slices_cur[[par]], slices_cur$ll, type = 'l',
         xlab = par, ylab = 'Log-Likelihood',
         main = par)
  }
  rm(par, slices_cur)
}
# In some cases, there is a difference of 6 ll points between 1.0 and 0.0, and in others it is only 2-3

# # How do season-specific values change with different theta_lambda1?:
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
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R120'))
# plot(pars_comp, pch = 20)
# # A bit of a relationship with Ri2/I20, but nothing terribly strong

# ---------------------------------------------------------------------------------------------------------------------

# Compile results of trajectory matching for theta_lambda1, fixing Ri2/I20 for 2010

# Set estpars:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Read in and format results:
res_temp <- load_and_format_mega_results(filepath = 'results/170921_fluA_thetalambda1_round1CI_fixRi2I20/',
                                         shared_estpars = shared_estpars,
                                         unit_estpars = unit_estpars,
                                         res_comp = res_df,
                                         run_name = 'Global_CI_fixRi2I20')
pars_top <- res_temp[[1]]
res_df <- res_temp[[2]]
pars_corr <- res_temp[[3]]

# Pairs plot:
pairs(pars_corr, pch = 20, main = 'theta_lambda1_CI')
# Here theta_lambda1=1.0 is the most likely value - but why?

# # Zoom in on theta_lambda1 and likelihood:
# ggplot(data = pars_top, aes(x = theta_lambda1, y = loglik)) + geom_point() + theme_classic()
# note that also highly correlated with I10 and R120 (higher I10 = lower theta_lambda1; higher R120 = higher theta_lambda1)
# strangely, this means that higher I10 associated with lower R120 - would expect this to lead to a giant outbreak

# Explore model fit - why fitting R120 high/I10 high/theta_lambda1=1.0?:
vir1 <- 'flu_A'
vir2 <- 'rsv'
yr <- 2010

debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

source('src/resp_interaction_model.R')

pars_temp <- pars_top %>%
  select(contains('2010_'),
         theta_lambda1) %>%
  rename('Ri1' = `2010_Ri1`,
         'I10' = `2010_I10`,
         'R10' = `2010_R10`,
         'R20' = `2010_R20`,
         'R120' = `2010_R120`)

# matplot(resp_mod@data %>% t(), pch = 20, type = 'b', lty = 1)
plot_list <- list()
for (i in 1:56) {
  coef(resp_mod, c('Ri2', 'I20')) <- c(1.59, 0.00137)
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p1 <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') + theme_classic() +
    labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p1
  print(sum(traj_temp$H1))
  # print(sum(traj_temp$H2))
  # print('')
}
do.call('grid.arrange', plot_list)

# slices <- slice_design(#center = pars_temp[1, ] %>% unlist(),
#   # center = tail(pars_temp, 1) %>% unlist(),
#   center = pars_temp[300, ] %>% unlist(),
#   Ri1 = seq(min(pars_temp$Ri1), max(pars_temp$Ri1), length.out = 10),
#   I10 = seq(min(pars_temp$I10), max(pars_temp$I10), length.out = 10),
#   R120 = seq(min(pars_temp$R120), 0.99, length.out = 10),
#   theta_lambda1 = seq(min(pars_temp$theta_lambda1), max(pars_temp$theta_lambda1), length.out = 10))
# for (s in unique(slices$slice)) {
#   print(s)
#   slices_temp <- slices %>% filter(slice == s) %>% select(-slice)
#   plot_list <- list()
#   
#   for (i in 1:nrow(slices_temp)) {
#     coef(resp_mod, c('Ri2', 'I20')) <- c(1.59, 0.00137)
#     coef(resp_mod, names(slices_temp)) <- slices_temp[i, ]
#     traj_temp <- trajectory(resp_mod, format = 'data.frame')
#     p1 <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
#       geom_line(aes(x = time, y = H2, group = .id), col = 'coral') + theme_classic() +
#       labs(x = 'Time', y = 'Cases', title = slices_temp[i, s]) + scale_y_sqrt()
#     if (i == 1 & s == 'theta_lambda1') {
#       plot(traj_temp$H2, type = 'l')
#     } else if (i != 1 & s == 'theta_lambda1') {
#       lines(traj_temp$H2, col = 'blue')
#     }
#     plot_list[[i]] <- p1
#   }
#   
#   do.call('grid.arrange', plot_list)
# }
# # Fitted range of Ri1 is actually really narrow (1.3767-1.3899)
# # Most I10 values lead to immediate drop off of cases
# # Higher R120 leads to lower and earlier peaks for both viruses
# # When Ri1/I10/R120 at highest likelihood values, little impact of theta_lambda1; when near "tail" of pars_top, stronger/lower theta_lambda1 leads to delayed peak

# Maybe suggesting that, with a high attack rate (20-25%), or just a really high start value, even a short overlap can delay outbreak?
# And needs a more powerful push downward when there are fewer people immune?

# I guess point is that it does seem more identifiable in a best-case scenario where we know the initial conditions,
# even if this particular scenario is not one that is terribly realistic

# # How do season-specific values change with different theta_lambda1?:
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
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R10'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R20'))
# plot(pars_comp, pch = 20)
# 
# pars_comp <- pars_top %>%
#   select(theta_lambda1, contains('R120'))
# plot(pars_comp, pch = 20)
# # A bit of a relationship with Ri2/I20, but nothing terribly strong

# ---------------------------------------------------------------------------------------------------------------------

# Compile results of trajectory matching for theta_lambda1, fixing Ri2/I20 for 2010, and starting runs 10 weeks earlier

# Set estpars:
shared_estpars <- c('theta_lambda1')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Read in and format results:
res_temp <- load_and_format_mega_results(filepath = 'results/240921_fluA_thetalambda1_round1CI_fixRi2I20_early/',
                                         shared_estpars = shared_estpars,
                                         unit_estpars = unit_estpars,
                                         res_comp = res_df,
                                         run_name = 'Global_CI_fixRi2I20_early')
pars_top <- res_temp[[1]]
res_df <- res_temp[[2]]
pars_corr <- res_temp[[3]]

# Pairs plot:
pairs(pars_corr, pch = 20, main = 'theta_lambda1_CI')
# Mostly fit to be 0.7 or above, with a couple of fits around 0.5-0.6ish
# But overall same pattern as when we started in week 40: lower theta_lambda1 can only explain pattern if I10 is unrealistically high

# Explore model fit:
vir1 <- 'flu_A'
vir2 <- 'rsv'
yr <- 2010

debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

early_start_val <- TRUE
source('src/resp_interaction_model.R')

pars_temp <- pars_top %>%
  select(contains('2010_'),
         theta_lambda1) %>%
  rename('Ri1' = `2010_Ri1`,
         'I10' = `2010_I10`,
         'R10' = `2010_R10`,
         'R20' = `2010_R20`,
         'R120' = `2010_R120`)

# matplot(resp_mod@data %>% t(), pch = 20, type = 'b', lty = 1)
plot_list <- list()
for (i in 1:56) {
  coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111) # c(1.59, 0.00137)
  coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
  # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
  traj_temp <- trajectory(resp_mod, format = 'data.frame')
  p1 <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
    geom_line(aes(x = time, y = H2, group = .id), col = 'coral') + theme_classic() +
    labs(x = 'Time', y = 'Cases')
  plot_list[[i]] <- p1
  print(sum(traj_temp$H1))
  # print(sum(traj_temp$H2))
  # print('')
}
do.call('grid.arrange', plot_list)

# ---------------------------------------------------------------------------------------------------------------------

# Compile results of trajectory matching for theta_lambda1/delta/rho2, fixing Ri2/I20 for 2010, and starting runs 10 weeks earlier

# Set estpars:
shared_estpars <- c('theta_lambda1', 'delta', 'rho2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Read in and format results:
res_temp <- load_and_format_mega_results(filepath = 'results/240921_fluA_thetalambda1_round1CI_fixRi2I20_early_rho2_delta/',
                                         shared_estpars = shared_estpars,
                                         unit_estpars = unit_estpars,
                                         res_comp = res_df,
                                         run_name = 'Global_CI_fixRi2I20_early_rho2delta')
pars_top <- res_temp[[1]]
res_df <- res_temp[[2]]
pars_corr <- res_temp[[3]]

# Plot results:
p1 <- ggplot(data = res_df, aes(x = year, y = value, fill = method)) + geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_brewer(palette = 'Set1')
print(p1)
# Maybe needs longer or more starting sets to run? Fits look really all over the place (or random effects format)

# Explore correlations:
pairs(pars_corr, pch = 20, main = 'theta_lambda1_CI')
# theta_lambda1 tends to fit low, but delta and rho2 highly uncertain

# # Explore model fit:
# vir1 <- 'flu_A'
# vir2 <- 'rsv'
# yr <- 2010
# 
# debug_bool <- FALSE
# 
# Ri_max1 <- 3.0
# Ri_max2 <- 3.0
# delta_min <- 7 / 60.0
# 
# early_start_val <- TRUE
# source('src/resp_interaction_model.R')
# 
# pars_temp <- pars_top %>%
#   select(contains('2010_'),
#          theta_lambda1,
#          delta,
#          rho2) %>%
#   rename('Ri1' = `2010_Ri1`,
#          'I10' = `2010_I10`,
#          'R10' = `2010_R10`,
#          'R20' = `2010_R20`,
#          'R120' = `2010_R120`)
# 
# matplot(resp_mod@data %>% t(), pch = 20, type = 'b', lty = 1)
# plot_list <- list()
# for (i in 1:56) {
#   coef(resp_mod, c('Ri2', 'I20')) <- c(1.61, 0.0000111) # c(1.59, 0.00137)
#   coef(resp_mod, names(pars_temp)) <- pars_temp[i, ]
#   # sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame')
#   traj_temp <- trajectory(resp_mod, format = 'data.frame')
#   p1 <- ggplot(data = traj_temp) + geom_line(aes(x = time, y = H1, group = .id), col = 'black') +
#     geom_line(aes(x = time, y = H2, group = .id), col = 'coral') + theme_classic() +
#     labs(x = 'Time', y = 'Cases')
#   plot_list[[i]] <- p1
#   print(sum(traj_temp$H1))
#   # print(sum(traj_temp$H2))
#   # print('')
# }
# do.call('grid.arrange', plot_list[!unlist(lapply(plot_list, is.null))])
# # Seems to have similar issues to other attempts to hold 2010 constant

# ---------------------------------------------------------------------------------------------------------------------

# Assess profile likelihood
# https://kingaa.github.io/sbied/pfilter/slides.pdf

# Get lists of results files for each theta_lambda1 value:
res_files_0 <- list.files(path = 'results/prof_0_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)
res_files_0.2 <- list.files(path = 'results/prof_2_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)
res_files_0.4 <- list.files(path = 'results/prof_4_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)
res_files_0.6 <- list.files(path = 'results/prof_6_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)
res_files_0.8 <- list.files(path = 'results/prof_8_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)
res_files_1 <- list.files(path = 'results/prof_10_160921_fluA_thetalambda1_round1ci_fitR_rho=15/', full.names = TRUE)

res_files <- list(res_files_0, res_files_0.2, res_files_0.4, res_files_0.6, res_files_0.8, res_files_1)

# Read in results:
pars_list <- list()
for (i in 1:length(res_files)) {
  
  # Load:
  res_full <- list()
  for (j in seq_along(res_files[[i]])) {
    res_full[[j]] <- read_rds(res_files[[i]][[j]])
  }
  
  # Compile:
  pars_temp <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist())
  expect_true(nrow(pars_temp) == length(res_files[[i]]))
  expect_true(all(is.finite(pars_temp$loglik)))
  
  # Organize:
  pars_temp <- pars_temp %>%
    arrange(desc(loglik))
  
  # Keep only top values:
  no_best <- nrow(subset(pars_temp, 2 * (max(loglik) - loglik) <= 1.92))
  pars_temp <- pars_temp[1:no_best, ]
  
  # Add theta_lambda1 value:
  pars_temp$theta_lambda1 <- seq(0, 1, by = 0.2)[i]
  
  # Add to list:
  pars_list[[i]] <- pars_temp
  
}

# Compile all results frames:
pars_top <- bind_rows(pars_list)

# Clean up:
rm(res_files, res_files_0, res_files_0.2, res_files_0.4, res_files_0.6, res_files_0.8, res_files_1,
   pars_list, pars_temp, res_full, i, j, no_best)

# Plot profile likelihood:
ggplot(data = pars_top %>% group_by(theta_lambda1) %>% summarise(loglik = max(loglik)),
       aes(x = theta_lambda1, y = loglik)) +
  geom_point() + geom_line() + theme_classic()

ggplot(data = pars_top, aes(x = theta_lambda1, y = loglik, group = theta_lambda1)) +
  geom_violin(fill = 'steelblue2') + theme_classic() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))

# Any differences in estimated parameter values?:
pars_top <- pars_top %>%
  pivot_longer(-c(theta_lambda1, loglik), names_to = 'param') %>%
  mutate(year = str_sub(param, 1, 4),
         param = str_sub(param, 6, str_length(param))) %>%
  select(year, param, value, theta_lambda1)

p1 <- ggplot(data = pars_top, aes(x = year, y = value, fill = theta_lambda1, group = paste(year, theta_lambda1))) +
  geom_boxplot() +
  facet_wrap(~param, scales = 'free_y') + theme_classic() + scale_fill_viridis()
print(p1)
# In general, as theta_lambda1 increases, fit values for R10/R20 decrease and for R120 increase, particularly for
# the pandemic season. Patterns for other parameters are either weak or vary from year to year.

# ---------------------------------------------------------------------------------------------------------------------

# # Close plots:
# dev.off()

# Clean up:
rm(list = ls())

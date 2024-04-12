# ---------------------------------------------------------------------------------------------------------------------
# Code to compare various solving algorithms for model fitting for both 1) time and 2) fit parameter values
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(gridExtra)

# Read in "main" results:
res_dir_sbplx <- 'results/round2_fit/round2_1_fluH1_plus_B/'

res_files <- list.files(path = res_dir_sbplx, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  mutate(num_fit = length(res_files)) %>%
  mutate(algorithm = 'sbplx')

expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

df_use <- pars_df %>% select(-c(loglik, message, niter, time, num_fit, algorithm)) %>% names() %>% length()
no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

pars_df <- pars_df %>%
  mutate(is_mle = no_best > 1,
         num_best = no_best)

pars_df <- pars_df %>%
  arrange(desc(loglik))
pars_top <- pars_df[1:25, ]

# Read in results using alternative algorithms:
res_dir_alt <- 'results/round2_fit/sens/check_algorithms/round2_1_fluH1_plus_B/'

res_files <- list.files(path = res_dir_alt, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df_alt <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  bind_cols('algorithm' = res_files %>%
              str_split('_') %>%
              purrr::map(~ .[[15]]) %>%
              str_split(fixed('.')) %>%
              purrr::map(~ .[[1]]) %>%
              unlist())

pars_df_alt <- pars_df_alt %>%
  mutate(algorithm = if_else(algorithm == 'sbplx', 'sbplx_check', algorithm))

pars_df_alt_LIST <- split(pars_df_alt, f = pars_df_alt$algorithm)

pars_df_alt_LIST <- lapply(pars_df_alt_LIST, function(ix) {
  
  no_best_ix <- nrow(subset(ix, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  
  ix %>%
    mutate(num_fit = nrow(ix),
           is_mle = no_best_ix > 1,
           num_best = no_best_ix)
  
})

pars_df_alt <- bind_rows(pars_df_alt_LIST)

pars_df_alt_LIST <- lapply(pars_df_alt_LIST, function(ix) {
  
  ix[1:25, ]
  
})
pars_top_alt <- bind_rows(pars_df_alt_LIST) %>%
  drop_na()
rm(pars_df_alt_LIST, df_use, no_best, res_dir_sbplx, res_dir_alt, res_full, res_files, i)

# Join all results:
pars_df <- pars_df %>%
  bind_rows(pars_df_alt)
pars_top <- pars_top %>%
  bind_rows(pars_top_alt)
rm(pars_df_alt, pars_top_alt)

# Compare:
pars_df %>%
  group_by(algorithm) %>%
  summarise(num_fit = unique(num_fit)) %>%
  print()

table(pars_df$algorithm, pars_df$message)

p_ll1 <- pars_df %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:unique(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_ll2 <- ggplot(data = pars_df, aes(x = loglik)) +
  geom_histogram(binwidth = 10000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Log-Likelihood', y = 'Count')

p_time1 <- pars_df %>%
  group_by(algorithm) %>%
  arrange(time) %>%
  mutate(x_val = 1:unique(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = time, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Time (Hours)', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_time2 <- ggplot(data = pars_df, aes(x = time)) +
  geom_histogram(binwidth = 0.1, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Time (Hours)', y = 'Count')

p_niter1 <- pars_df %>%
  group_by(algorithm) %>%
  arrange(niter) %>%
  mutate(x_val = 1:unique(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = niter, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Number of Iterations', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_niter2 <- ggplot(data = pars_df, aes(x = niter)) +
  geom_histogram(binwidth = 10000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Number of Iterations', y = 'Count')

p_niterpertime1 <- ggplot(data = pars_df, aes(x = time, y = niter, col = algorithm)) +
  geom_point() +
  theme_classic() +
  labs(x = 'Time (Hours)', y = 'Number of Iterations', col = 'Algorithm') +
  scale_color_brewer(palette = 'Set1')
p_niterpertime2 <- ggplot(data = pars_df, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

# print(p_ll1)
# print(p_ll2)
# 
# print(p_time1)
# print(p_time2)
# 
# print(p_niter1)
# print(p_niter2)
# 
# print(p_niterpertime1)
# print(p_niterpertime2)

p_parms <- pars_df %>%
  select(rho1:eta_ah2, algorithm) %>%
  pivot_longer(-algorithm,
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot() +
  geom_violin(aes(x = algorithm, y = value), fill = 'gray95') +
  facet_wrap(~ parameter, scales = 'free') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65)) +
  labs(x = 'Algorithm', y = 'Parameter Value')

# Compare only best-fit:
table(pars_top$algorithm, pars_top$message)

p_ll1_top_a <- pars_top %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_ll1_top_b <- pars_top %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(limits = c(max(pars_top$loglik) - 10000, max(pars_top$loglik) + 10)) +
  scale_color_brewer(palette = 'Set1')

p_time1_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(time) %>%
  mutate(x_val = 1:length(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = time, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Time (Hours)', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niter1_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(niter) %>%
  mutate(x_val = 1:length(num_fit)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = niter, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Number of Iterations', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niterpertime2_top <- ggplot(data = pars_top, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

p_parms_top <- pars_top %>%
  filter(str_detect(message, 'XTOL') | str_detect(message, 'success')) %>%
  select(rho1:eta_ah2, algorithm) %>%
  pivot_longer(-algorithm,
               names_to = 'parameter',
               values_to = 'value') %>%
  ggplot() +
  geom_violin(aes(x = algorithm, y = value), fill = 'gray95') +
  facet_wrap(~ parameter, scales = 'free') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.65)) +
  labs(x = 'Algorithm', y = 'Parameter Value')

plot(arrangeGrob(p_ll1, p_ll1_top_b), nrow = 1)
plot(arrangeGrob(arrangeGrob(p_time1, p_niter1, nrow = 1), p_niterpertime2, ncol = 1))
plot(arrangeGrob(arrangeGrob(p_time1_top, p_niter1_top, nrow = 1), p_niterpertime2_top, ncol = 1))
print(p_parms)
print(p_parms_top)

# Clean up:
rm(p_ll1, p_ll1_top_a, p_ll1_top_b, p_ll2, p_niter1, p_niter1_top, p_niter2, p_time1, p_time1_top, p_time2,
   p_niterpertime1, p_niterpertime2, p_niterpertime2_top, p_parms, p_parms_top)

# ---------------------------------------------------------------------------------------------------------------------

# Based on this, try fitting using "praxis" algorithm, and compare MLEs (Hong Kong)

# Read in results using "sbplx" algorithm:
res_dir_sbplx <- 'results/round2_fit/round2_3_fluH1_plus_B/'

res_files <- list.files(path = res_dir_sbplx, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  mutate(algorithm = 'sbplx')

expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

df_use <- pars_df %>% select(-c(loglik, message, niter, time, algorithm)) %>% names() %>% length()
no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

pars_df <- pars_df %>%
  arrange(desc(loglik))

pars_top <- pars_df[1:no_best, ]

pars_top <- pars_top %>%
  filter(!str_detect(message, 'maxtime'))

# Read in results using "praxis" algorithm:
res_dir_alt <- 'results/round2_fit/sens/check_algorithms/round2_6_fluH1_plus_B/'

res_files <- list.files(path = res_dir_alt, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df_alt <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  mutate(algorithm = 'praxis')

expect_true(nrow(pars_df_alt) == length(res_files))
expect_true(all(is.finite(pars_df_alt$loglik)))

df_use <- pars_df_alt %>% select(-c(loglik, message, niter, time, algorithm)) %>% names() %>% length()
no_best <- nrow(subset(pars_df_alt, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

pars_df_alt <- pars_df_alt %>%
  arrange(desc(loglik))

pars_top_alt <- pars_df_alt[1:no_best, ]

pars_top_alt <- pars_top_alt %>%
  filter(!str_detect(message, 'maxtime'))

rm(df_use, no_best, res_dir_sbplx, res_dir_alt, res_full, res_files, i)

# Join all results:
pars_df <- pars_df %>%
  bind_rows(pars_df_alt)
pars_top <- pars_top %>%
  bind_rows(pars_top_alt)
rm(pars_df_alt, pars_top_alt)

# Compare:
p_ll <- pars_df %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_ll_top1 <- pars_top %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_ll_top2 <- pars_top %>%
  ggplot(aes(x = algorithm, y = loglik)) +
  geom_violin(fill = 'gray95') +
  theme_classic() +
  labs(x = 'Algorithm', y = 'Log-Likelihood')

p_niterpertime <- ggplot(data = pars_df, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

p_time_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(time) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = time, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Time (Hours)', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niter_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(niter) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = niter, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Number of Iterations', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niterpertime_top <- ggplot(data = pars_top, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

pars_top <- pars_top %>%
  mutate(is_mle = c(TRUE, rep(FALSE, 6), TRUE, rep(FALSE, 39))) %>%
  filter(str_detect(message, 'XTOL') | str_detect(message, 'success')) %>%
  select(rho1:eta_ah2, algorithm, is_mle) %>%
  pivot_longer(-c(algorithm, is_mle),
               names_to = 'parameter',
               values_to = 'value')

cis <- read_csv('results/MLE_plus_95CI_from_bootstrapping_HPDI.csv') %>%
  filter(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2'))

p_parms_top <- ggplot() +
  geom_violin(data = pars_top, aes(x = algorithm, y = value), fill = 'gray95') +
  geom_point(data = pars_top %>% filter(is_mle), aes(x = algorithm, y = value), shape = '\u2605', size = 4) +
  geom_linerange(data = cis, x = 'sbplx', aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ parameter, scales = 'free') +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.65)) +
  labs(x = 'Algorithm', y = 'Parameter Value')

plot(arrangeGrob(p_ll_top2, arrangeGrob(p_ll, p_ll_top1, nrow = 2), ncol = 2))
print(p_parms_top)
plot(arrangeGrob(arrangeGrob(p_time_top, p_niter_top, nrow = 1), p_niterpertime_top, p_niterpertime, ncol = 1))

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------

# Also check for Canada

# Read in results using "sbplx" algorithm:
res_dir_sbplx <- 'results/round2_fit/sens/canada/round2_3_flu/'

res_files <- list.files(path = res_dir_sbplx, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  mutate(algorithm = 'sbplx')

expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

df_use <- pars_df %>% select(-c(loglik, message, niter, time, algorithm)) %>% names() %>% length()
no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

pars_df <- pars_df %>%
  arrange(desc(loglik))

pars_top <- pars_df[1:no_best, ]

pars_top <- pars_top %>%
  filter(!str_detect(message, 'maxtime'))

# Read in results using "praxis" algorithm:
res_dir_alt <- 'results/round2_fit/sens/canada/check_algorithms/round2_2_flu/'

res_files <- list.files(path = res_dir_alt, full.names = TRUE)

res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df_alt <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  bind_cols('message' = lapply(res_full, getElement, 'message') %>%
              unlist()) %>%
  bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
              unlist()) %>%
  bind_cols('time' = lapply(res_full, getElement, 'etime') %>%
              unlist()) %>%
  mutate(algorithm = 'praxis')

expect_true(nrow(pars_df_alt) == length(res_files))
expect_true(all(is.finite(pars_df_alt$loglik)))

df_use <- pars_df_alt %>% select(-c(loglik, message, niter, time, algorithm)) %>% names() %>% length()
no_best <- nrow(subset(pars_df_alt, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

pars_df_alt <- pars_df_alt %>%
  arrange(desc(loglik))

pars_top_alt <- pars_df_alt[1:no_best, ]

pars_top_alt <- pars_top_alt %>%
  filter(!str_detect(message, 'maxtime'))

rm(df_use, no_best, res_dir_sbplx, res_dir_alt, res_full, res_files, i)

# Join all results:
pars_df <- pars_df %>%
  bind_rows(pars_df_alt)
pars_top <- pars_top %>%
  bind_rows(pars_top_alt)
rm(pars_df_alt, pars_top_alt)

# Compare:
p_ll <- pars_df %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_ll_top1 <- pars_top %>%
  group_by(algorithm) %>%
  arrange(desc(loglik)) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = loglik, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Log-Likelihood', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')
p_ll_top2 <- pars_top %>%
  ggplot(aes(x = algorithm, y = loglik)) +
  geom_violin(fill = 'gray95') +
  theme_classic() +
  labs(x = 'Algorithm', y = 'Log-Likelihood')

p_niterpertime <- ggplot(data = pars_df, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

p_time_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(time) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = time, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Time (Hours)', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niter_top <- pars_top %>%
  group_by(algorithm) %>%
  arrange(niter) %>%
  mutate(x_val = 1:length(rho1)) %>%
  ungroup() %>%
  ggplot(aes(x = x_val, y = niter, group = algorithm, col = algorithm)) +
  geom_line() +
  geom_point() +
  theme_classic() +
  labs(x = '', y = 'Number of Iterations', col = 'Algorithm') +
  scale_x_continuous(breaks = NULL) +
  scale_color_brewer(palette = 'Set1')

p_niterpertime_top <- ggplot(data = pars_top, aes(x = niter / time)) +
  geom_histogram(binwidth = 5000, col = 'black', fill = 'gray90') +
  facet_wrap(~ algorithm, scales = 'free_y', nrow = 1) +
  theme_classic() +
  labs(x = 'Iterations per Hour', y = 'Count')

pars_top <- pars_top %>%
  mutate(is_mle = c(TRUE, rep(FALSE, 183), TRUE, rep(FALSE, 108))) %>%
  filter(str_detect(message, 'XTOL') | str_detect(message, 'success')) %>%
  select(rho1:phi2, algorithm, is_mle) %>%
  pivot_longer(-c(algorithm, is_mle),
               names_to = 'parameter',
               values_to = 'value')

cis <- read_csv('results/round2_fit/sens/canada/MLE_plus_95CI_from_bootstrapping_HPDI.csv') %>%
  filter(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2'))

p_parms_top <- ggplot() +
  geom_violin(data = pars_top, aes(x = algorithm, y = value), fill = 'gray95') +
  geom_point(data = pars_top %>% filter(is_mle), aes(x = algorithm, y = value), shape = '\u2605', size = 4) +
  geom_linerange(data = cis, x = 'sbplx', aes(ymin = lower, ymax = upper)) +
  facet_wrap(~ parameter, scales = 'free') +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.65)) +
  labs(x = 'Algorithm', y = 'Parameter Value')

plot(arrangeGrob(p_ll_top2, arrangeGrob(p_ll, p_ll_top1, nrow = 2), ncol = 2))
print(p_parms_top)
plot(arrangeGrob(arrangeGrob(p_time_top, p_niter_top, nrow = 1), p_niterpertime_top, p_niterpertime, ncol = 1))

# Clean up:
rm(list = ls())

# ---------------------------------------------------------------------------------------------------------------------
# Process results of parametric bootstrapping and output 95% CIs and plots
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Get names of all results files:
file_list_h1 <- list.files(path = 'results/bootstrapping/flu_H1/', full.names = TRUE)
file_list_b <- list.files(path = 'results/bootstrapping/flu_b/', full.names = TRUE)

# Ensure no results missing:
expect_true(length(file_list_h1) == 500 * 10)
expect_true(length(file_list_b) == 500 * 10)

# Read in all results:
res_full_h1 = list()
for (i in seq_along(file_list_h1)) {
  res_full_h1[[i]] <- read_rds(file_list_h1[[i]])
}
res_full_b = list()
for (i in seq_along(file_list_b)) {
  res_full_b[[i]] <- read_rds(file_list_b[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
res_df_h1 <- lapply(res_full_h1, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full_h1, getElement, 'll') %>%
              unlist()) %>%
  mutate(dataset = str_split(file_list_h1, '_') %>% map(~ .x[6]) %>% unlist()) %>%
  bind_cols('conv' = lapply(res_full_h1, getElement, 'message') %>%
              unlist())
expect_true(nrow(res_df_h1) == length(file_list_h1))
expect_true(all(is.finite(res_df_h1$loglik)))

res_df_b <- lapply(res_full_b, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full_b, getElement, 'll') %>%
              unlist()) %>%
  mutate(dataset = str_split(file_list_b, '_') %>% map(~ .x[6]) %>% unlist()) %>%
  bind_cols('conv' = lapply(res_full_b, getElement, 'message') %>%
              unlist())
expect_true(nrow(res_df_b) == length(file_list_b))
expect_true(all(is.finite(res_df_b$loglik)))

# Remove estimates that did not converge:
table(res_df_h1$conv)
table(res_df_b$conv)

res_df_h1 <- res_df_h1 %>%
  filter(str_detect(conv, 'XTOL_REACHED')) %>%
  select(-conv)
res_df_b <- res_df_b %>%
  filter(str_detect(conv, 'XTOL_REACHED')) %>%
  select(-conv)

# Check whether at least one estimate remains for all synthetic datasets:
expect_true(length(unique(res_df_h1$dataset)) == 500)
expect_true(length(unique(res_df_b$dataset)) == 500)

# Keep only top fit for each synthetic dataset:
res_df_h1 <- res_df_h1 %>%
  group_by(dataset) %>%
  filter(loglik == max(loglik)) %>%
  ungroup()
res_df_b <- res_df_b %>%
  group_by(dataset) %>%
  filter(loglik == max(loglik)) %>%
  ungroup()
expect_true(nrow(res_df_h1) == 500)
expect_true(nrow(res_df_b) == 500)

# Are all top fits within a reasonable range of log-likelihood values?
hist(res_df_h1$loglik, breaks = 50)
hist(res_df_b$loglik, breaks = 50)

# Calculate 99% confidence intervals for each parameter:
res_df_h1_long <- res_df_h1 %>%
  select(-c(loglik, dataset)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = factor(parameter, levels = names(res_df_h1)[1:47]))

ci_res <- res_df_h1_long %>%
  group_by(parameter) %>%
  summarise(lower = quantile(value, p = 0.005),
            upper = quantile(value, p = 0.995)) %>%
  mutate(vir1 = 'flu_h1')

res_df_b_long <- res_df_b %>%
  select(-c(loglik, dataset)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = factor(parameter, levels = names(res_df_b)[1:47]))

ci_res <- ci_res %>%
  bind_rows(res_df_b_long %>%
              group_by(parameter) %>%
              summarise(lower = quantile(value, p = 0.005),
                        upper = quantile(value, p = 0.995)) %>%
              mutate(vir1 = 'flu_b'))

# Write results to file:
write_csv(ci_res, file = 'results/99CI_from_boostrapping.csv')

# Plot range of fit values for each parameter:
res_df_h1_NOSHARED <- res_df_h1_long %>%
  filter(!(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2'))) %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_remove(parameter, paste0(season, '_'))) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')))

p1_h1 <- ggplot(data = res_df_h1_long %>% filter(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')),
                aes(x = value, y = ..count../nrow(res_df_h1))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2_h1 <- ggplot(data = res_df_h1_NOSHARED, aes(x = value, y = ..count../nrow(res_df_h1), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1_h1)
print(p2_h1)

res_df_b_NOSHARED <- res_df_b_long %>%
  filter(!(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2'))) %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_remove(parameter, paste0(season, '_'))) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')))

p1_b <- ggplot(data = res_df_b_long %>% filter(parameter %in% c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')),
               aes(x = value, y = ..count../nrow(res_df_h1))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2_b <- ggplot(data = res_df_b_NOSHARED, aes(x = value, y = ..count../nrow(res_df_h1), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1_b)
print(p2_b)

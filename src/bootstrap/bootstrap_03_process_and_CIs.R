# ---------------------------------------------------------------------------------------------------------------------
# Process results of parametric bootstrapping and output 95% CIs and plots
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)
library(rethinking)
library(gt)

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
  mutate(dataset = str_split(file_list_h1, '_') %>% purrr::map(~ .x[6]) %>% unlist()) %>%
  bind_cols('conv' = lapply(res_full_h1, getElement, 'message') %>%
              unlist())
expect_true(nrow(res_df_h1) == length(file_list_h1))
expect_true(all(is.finite(res_df_h1$loglik)))

res_df_b <- lapply(res_full_b, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full_b, getElement, 'll') %>%
              unlist()) %>%
  mutate(dataset = str_split(file_list_b, '_') %>% purrr::map(~ .x[6]) %>% unlist()) %>%
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

# # Check that no states go below zero for any of the top-fit parameter sets:
# unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
# shared_estpars <- res_df_h1 %>% select(!contains(unit_estpars) & !'loglik' & !'dataset') %>% names()
# true_estpars <- c(shared_estpars, unit_estpars)
# 
# prof_lik <- FALSE
# lag_val <- 0
# 
# vir1 <- 'flu_h1'
# source('src/functions/setup_global_likelilhood.R')
# 
# traj_list <- lapply(1:length(seasons), function(ix) {
#   pars_temp <- res_df_h1 %>%
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
# 
# vir1 <- 'flu_b'
# source('src/functions/setup_global_likelilhood.R')
# 
# traj_list <- lapply(1:length(seasons), function(ix) {
#   pars_temp <- res_df_h1 %>%
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
res_df_h1 <- res_df_h1 %>%
  mutate(delta2 = d2 * delta1,
         `s13-14_R10 + R120` = `s13-14_R10` + `s13-14_R120`,
         `s14-15_R10 + R120` = `s14-15_R10` + `s14-15_R120`,
         `s15-16_R10 + R120` = `s15-16_R10` + `s15-16_R120`,
         `s16-17_R10 + R120` = `s16-17_R10` + `s16-17_R120`,
         `s17-18_R10 + R120` = `s17-18_R10` + `s17-18_R120`,
         `s18-19_R10 + R120` = `s18-19_R10` + `s18-19_R120`,
         `s13-14_R20 + R120` = `s13-14_R20` + `s13-14_R120`,
         `s14-15_R20 + R120` = `s14-15_R20` + `s14-15_R120`,
         `s15-16_R20 + R120` = `s15-16_R20` + `s15-16_R120`,
         `s16-17_R20 + R120` = `s16-17_R20` + `s16-17_R120`,
         `s17-18_R20 + R120` = `s17-18_R20` + `s17-18_R120`,
         `s18-19_R20 + R120` = `s18-19_R20` + `s18-19_R120`) %>%
  select(rho1:`s18-19_R120`, delta2:`s18-19_R20 + R120`, loglik:dataset)
res_df_b <- res_df_b %>%
  mutate(delta2 = d2 * delta1,
         `s13-14_R10 + R120` = `s13-14_R10` + `s13-14_R120`,
         `s14-15_R10 + R120` = `s14-15_R10` + `s14-15_R120`,
         `s15-16_R10 + R120` = `s15-16_R10` + `s15-16_R120`,
         `s16-17_R10 + R120` = `s16-17_R10` + `s16-17_R120`,
         `s17-18_R10 + R120` = `s17-18_R10` + `s17-18_R120`,
         `s18-19_R10 + R120` = `s18-19_R10` + `s18-19_R120`,
         `s13-14_R20 + R120` = `s13-14_R20` + `s13-14_R120`,
         `s14-15_R20 + R120` = `s14-15_R20` + `s14-15_R120`,
         `s15-16_R20 + R120` = `s15-16_R20` + `s15-16_R120`,
         `s16-17_R20 + R120` = `s16-17_R20` + `s16-17_R120`,
         `s17-18_R20 + R120` = `s17-18_R20` + `s17-18_R120`,
         `s18-19_R20 + R120` = `s18-19_R20` + `s18-19_R120`) %>%
  select(rho1:`s18-19_R120`, delta2:`s18-19_R20 + R120`, loglik:dataset)

# Calculate 95% confidence intervals for each parameter:
res_df_h1_long <- res_df_h1 %>%
  select(-c(loglik, dataset)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = factor(parameter, levels = names(res_df_h1)[1:(length(names(res_df_h1)) - 2)]))

ci_res <- res_df_h1_long %>%
  group_by(parameter) %>%
  summarise(lower = HPDI(value, p = 0.95)[1],
            upper = HPDI(value, p = 0.95)[2]) %>%
  mutate(vir1 = 'flu_h1')
# ci_res <- res_df_h1_long %>%
#   group_by(parameter) %>%
#   summarise(lower = quantile(value, p = 0.025),
#             upper = quantile(value, p = 0.975)) %>%
#   mutate(vir1 = 'flu_h1')

res_df_b_long <- res_df_b %>%
  select(-c(loglik, dataset)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = factor(parameter, levels = names(res_df_b)[1:(length(names(res_df_b)) - 2)]))

ci_res <- ci_res %>%
  bind_rows(res_df_b_long %>%
              group_by(parameter) %>%
              summarise(lower = HPDI(value, p = 0.95)[1],
                        upper = HPDI(value, p = 0.95)[2]) %>%
              mutate(vir1 = 'flu_b'))
# ci_res <- ci_res %>%
#   bind_rows(res_df_b_long %>%
#               group_by(parameter) %>%
#               summarise(lower = quantile(value, p = 0.025),
#                         upper = quantile(value, p = 0.975)) %>%
#               mutate(vir1 = 'flu_b'))

# Write results to file:
write_csv(ci_res, file = 'results/95CI_from_boostrapping_HPDI.csv')

# Read in MLEs and add to data frame:
mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
mle_b <- read_rds('results/MLEs_flu_b.rds')

mle_h1 <- mle_h1 %>%
  mutate(delta2 = d2 * delta1,
         `s13-14_R10 + R120` = `s13-14_R10` + `s13-14_R120`,
         `s14-15_R10 + R120` = `s14-15_R10` + `s14-15_R120`,
         `s15-16_R10 + R120` = `s15-16_R10` + `s15-16_R120`,
         `s16-17_R10 + R120` = `s16-17_R10` + `s16-17_R120`,
         `s17-18_R10 + R120` = `s17-18_R10` + `s17-18_R120`,
         `s18-19_R10 + R120` = `s18-19_R10` + `s18-19_R120`,
         `s13-14_R20 + R120` = `s13-14_R20` + `s13-14_R120`,
         `s14-15_R20 + R120` = `s14-15_R20` + `s14-15_R120`,
         `s15-16_R20 + R120` = `s15-16_R20` + `s15-16_R120`,
         `s16-17_R20 + R120` = `s16-17_R20` + `s16-17_R120`,
         `s17-18_R20 + R120` = `s17-18_R20` + `s17-18_R120`,
         `s18-19_R20 + R120` = `s18-19_R20` + `s18-19_R120`)
mle_b <- mle_b %>%
  mutate(delta2 = d2 * delta1,
         `s13-14_R10 + R120` = `s13-14_R10` + `s13-14_R120`,
         `s14-15_R10 + R120` = `s14-15_R10` + `s14-15_R120`,
         `s15-16_R10 + R120` = `s15-16_R10` + `s15-16_R120`,
         `s16-17_R10 + R120` = `s16-17_R10` + `s16-17_R120`,
         `s17-18_R10 + R120` = `s17-18_R10` + `s17-18_R120`,
         `s18-19_R10 + R120` = `s18-19_R10` + `s18-19_R120`,
         `s13-14_R20 + R120` = `s13-14_R20` + `s13-14_R120`,
         `s14-15_R20 + R120` = `s14-15_R20` + `s14-15_R120`,
         `s15-16_R20 + R120` = `s15-16_R20` + `s15-16_R120`,
         `s16-17_R20 + R120` = `s16-17_R20` + `s16-17_R120`,
         `s17-18_R20 + R120` = `s17-18_R20` + `s17-18_R120`,
         `s18-19_R20 + R120` = `s18-19_R20` + `s18-19_R120`)

mle <- bind_rows(mle_h1[1, ], mle_b[1, ]) %>%
  mutate(vir1 = c('flu_h1', 'flu_b')) %>%
  pivot_longer(cols = !vir1,
               names_to = 'parameter',
               values_to = 'mle')

ci_res <- ci_res %>%
  left_join(mle, by = c('vir1', 'parameter')) %>%
  select(parameter, mle, lower:vir1)

# Write results to file:
write_csv(ci_res, file = 'results/MLE_plus_95CI_from_boostrapping_HPDI.csv')

# Generate tables of results:
table_h1 <- ci_res %>%
  filter(vir1 == 'flu_h1') %>%
  select(-vir1) %>%
  gt() %>%
  tab_header(title = 'Best-Fit Parameter Values (Flu H1)') %>%
  fmt_number(columns = c(mle, lower, upper), decimals = 3, suffixing = TRUE)
table_b <- ci_res %>%
  filter(vir1 == 'flu_b') %>%
  select(-vir1) %>%
  gt() %>%
  tab_header(title = 'Best-Fit Parameter Values (Flu B)') %>%
  fmt_number(columns = c(mle, lower, upper), decimals = 3, suffixing = TRUE)

print(table_h1)
print(table_b)

gtsave(table_h1, filename = 'results/plots/table_CIs_h1.html')
gtsave(table_b, filename = 'results/plots/table_CIs_b.html')

# Check whether MLEs fall within CIs:
ci_res %>% filter(mle <= lower) # several season-specific parameters for B (although they are close to the boundaries)
ci_res %>% filter(mle >= upper)

# Plot range of fit values for each parameter:
composite_params <- c('delta2', res_df_h1_long %>%
                        filter(grepl(' ', parameter)) %>%
                        pull(parameter) %>%
                        as.character() %>%
                        unique())

res_df_h1_NOSHARED <- res_df_h1_long %>%
  filter(!(parameter %in% composite_params)) %>%
  filter(!(parameter %in% shared_estpars)) %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_remove(parameter, paste0(season, '_'))) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')))

p1_h1 <- ggplot(data = res_df_h1_long %>% filter(parameter %in% shared_estpars),
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
  filter(!(parameter %in% composite_params)) %>%
  filter(!(parameter %in% shared_estpars)) %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_remove(parameter, paste0(season, '_'))) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')))

p1_b <- ggplot(data = res_df_b_long %>% filter(parameter %in% shared_estpars),
               aes(x = value, y = ..count../nrow(res_df_h1))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2_b <- ggplot(data = res_df_b_NOSHARED, aes(x = value, y = ..count../nrow(res_df_h1), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1_b)
print(p2_b)

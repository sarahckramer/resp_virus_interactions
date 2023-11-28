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

# Ensure no results missing:
expect_true(length(file_list) == 500 * 10)

# Read in all results:
res_full = list()
for (i in seq_along(file_list)) {
  res_full[[i]] <- read_rds(file_list[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
res_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist()) %>%
  mutate(dataset = str_split(file_list, '_') %>% purrr::map(~ .x[which(!is.na(as.numeric(str_split(file_list, '_')[[1]])))]) %>% unlist()) %>%
  bind_cols('conv' = lapply(res_full, getElement, 'message') %>%
              unlist())
expect_true(nrow(res_df) == length(file_list))
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
# lag_val <- 0
# 
# vir1 <- 'flu_h1_plus_b'
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
res_df <- res_df %>%
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
res_df_long <- res_df %>%
  select(-c(loglik, dataset)) %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'value') %>%
  mutate(parameter = factor(parameter, levels = names(res_df)[1:(length(names(res_df)) - 2)]))

ci_res <- res_df_long %>%
  group_by(parameter) %>%
  summarise(lower = HPDI(value, p = 0.95)[1],
            upper = HPDI(value, p = 0.95)[2]) %>%
  mutate(vir1 = 'flu_h1_plus_b')
# ci_res <- res_df_long %>%
#   group_by(parameter) %>%
#   summarise(lower = quantile(value, p = 0.025),
#             upper = quantile(value, p = 0.975)) %>%
#   mutate(vir1 = 'flu_h1_plus_b')

# Write results to file:
write_csv(ci_res, file = 'results/95CI_from_boostrapping_HPDI.csv')

# Read in MLEs and add to data frame:
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')

mle <- mle %>%
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

mle <- mle[1, ] %>%
  pivot_longer(cols = everything(),
               names_to = 'parameter',
               values_to = 'mle')

ci_res <- ci_res %>%
  left_join(mle, by = 'parameter') %>%
  select(parameter, mle, lower:vir1)

# Write results to file:
write_csv(ci_res, file = 'results/MLE_plus_95CI_from_boostrapping_HPDI.csv')

# Generate tables of results:
res_table <- ci_res %>%
  select(-vir1) %>%
  gt() %>%
  tab_header(title = 'Best-Fit Parameter Values') %>%
  fmt_number(columns = c(mle, lower, upper), decimals = 3, suffixing = TRUE)
print(res_table)

gtsave(res_table, filename = 'results/plots/table_CIs_h1_plus_b.html')

# Check whether MLEs fall within CIs:
ci_res %>% filter(mle <= lower)
ci_res %>% filter(mle >= upper)

# Plot range of fit values for each parameter:
composite_params <- c('delta2', res_df_long %>%
                        filter(grepl(' ', parameter)) %>%
                        pull(parameter) %>%
                        as.character() %>%
                        unique())

res_df_NOSHARED <- res_df_long %>%
  filter(!(parameter %in% composite_params)) %>%
  filter(!(parameter %in% shared_estpars)) %>%
  mutate(season = str_sub(parameter, 1, 6),
         parameter = str_remove(parameter, paste0(season, '_'))) %>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')))

p1 <- ggplot(data = res_df_long %>% filter(parameter %in% shared_estpars),
             aes(x = value, y = after_stat(count)/nrow(res_df))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2 <- ggplot(data = res_df_NOSHARED, aes(x = value, y = after_stat(count)/nrow(res_df), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1)
print(p2)

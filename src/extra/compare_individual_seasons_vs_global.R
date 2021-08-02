# ---------------------------------------------------------------------------------------------------------------------
# Code to confirm that using global likelihood with no global parameters can achieve the same results as fitting
# each season separately (flu_B only)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

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

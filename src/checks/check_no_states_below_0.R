# ---------------------------------------------------------------------------------------------------------------------
# Code to check that best-fit parameter values do not lead trajectories to drop below 0
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Specify location of results to check:
res_dir <- 'results/round2_fit/round2_3_fluH1/'

# Check for missing results files:
res_files <- list.files(path = res_dir, full.names = FALSE)
res_exist <- res_files %>%
  str_split('_') %>%
  map(~ .x[length(.x)]) %>%
  str_split(fixed('.')) %>%
  map(~ .x[1]) %>%
  unlist() %>%
  as.numeric() %>%
  sort()
print(which(!(1:500 %in% res_exist)))

# Read in and compile all results:
res_files <- list.files(path = res_dir, full.names = TRUE)
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}

pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Limit to top results:
pars_df <- pars_df %>%
  arrange(desc(loglik))

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
shared_estpars <- pars_df %>% select(!contains(unit_estpars) & !'loglik') %>% names()
expect_equal(length(shared_estpars), 12)
true_estpars <- c(shared_estpars, unit_estpars)

df_use <- pars_df %>% select(-loglik) %>% names() %>% length()
# expect_equal(df_use, 54)

no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
no_best <- max(no_best, 50)

pars_top <- pars_df[1:no_best, ]

# Get/set relevant model parameters:
if (str_detect(res_dir, 'H1_plus_B')) {
  vir1 <- 'flu_h1_plus_b'
} else if (str_detect(res_dir, 'H1')) {
  vir1 <- 'flu_h1'
} else if (str_detect(res_dir, 'B')) {
  vir1 <- 'flu_b'
} else {
  print('Invalid flu subtype!')
}

prof_lik <- FALSE
lag_val <- 0

# Load pomp objects:
if (!any(str_detect(names(pars_df), 's13-14'))) {
  sens <- 'less_circ_h3'
}
source('src/functions/setup_global_likelilhood.R')

# Run trajectories for all seasons:
traj_list <- lapply(1:length(seasons), function(ix) {
  pars_temp <- pars_top %>%
    select(all_of(shared_estpars), contains(seasons[ix]))
  names(pars_temp) <- true_estpars
  
  p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(pars_temp))
  for (param in names(pars_temp)) {
    p_mat[param, ] <- pars_temp %>% pull(param)
  }
  
  trajectory(object = po_list[[ix]],
             params = p_mat,
             format = 'data.frame') %>%
    select(!(H1_tot:H2_tot)) %>%
    pivot_longer(X_SS:H2,
                 names_to = 'state')
})

# Check for 0s in state variables:
expect_false(any(lapply(traj_list, function(ix) {
  any(ix[, 'value'] < 0)
}) %>%
  unlist()))

# Clean up:
rm(list = ls())

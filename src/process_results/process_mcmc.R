# ---------------------------------------------------------------------------------------------------------------------
# Code to process results obtained using MCMC (using BayesianTools)
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(BayesianTools)
library(MCMCvis)
library(GGally)
library(rethinking)

# Set MCMC parameters:
n_chains <- 5
n_iter <- 1e6

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'
sens <- 'main'
fit_canada <- FALSE

if (fit_canada) {
  vir1 <- 'flu'
} else {
  vir1 <- 'flu_h1_plus_b'
}

seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
if (fit_canada) {
  seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
}

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

# Get parameters estimated:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

unit_sp_estpars <- c()
for (i in 1:length(seasons)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons[i], unit_estpars, sep = '_'))
}
rm(i)

true_estpars <- c(shared_estpars, unit_estpars)
estpars <- c(shared_estpars, unit_sp_estpars)

# ---------------------------------------------------------------------------------------------------------------------

# Functions

create_obj_fxn <- function(po, estpars) {
  # Creates the objective function for a given season
  # param po: A pomp model object for a specific season
  # param estpars: A vector listing the parameters to be fit
  # returns: The negative log likelihood
  
  ofun <- traj_objfun(data = po, 
                      est = estpars, 
                      partrans = po@partrans,
                      verbose=TRUE)
  
  return(ofun)
}

calculate_global_loglik <- function(trans_vals) {
  # Calculates the log-likelihood for each season, and combines to yield a global log-likelihood
  # param trans_vals: Unnamed vector of transformed parameters; fxn only works if this is the only input?
  # returns: The global, negative log-likelihood
  
  # Add names to vector:
  if(is.null(names(trans_vals))) names(trans_vals) <- x0_trans_names
  
  # Split params into shared and unit params:
  unit_in <- as.data.frame(sapply(paste0('^', seasons, '_'), grep, x = names(trans_vals)))
  if (ncol(unit_in) == 1) {
    unit_in <- as.data.frame(matrix(data = c(unlist(unit_in)), nrow = 1))
  }
  names(unit_in) <- seasons
  shared_params <- trans_vals[-unique(unlist(unit_in))]
  
  unit_params <- list()
  for (i in 1:length(seasons)) {
    unit <- trans_vals[unit_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 2]
      unit_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season:
  units_ll <- rep(NA, length(seasons))
  
  for (i in 1:length(seasons)) {
    if (length(unit_params) > 0) {
      params_temp <- c(shared_params, unit_params[[i]])
    } else {
      params_temp <- c(shared_params)
    }
    
    units_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Calculate global -ll:
  glob_ll <- sum(units_ll)
  return(-glob_ll)
}

transform_params <- function(orig_vals, po, seas, params_all, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_shared <- orig_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared) <- orig_vals_shared
  trans_vals_shared <- coef(po, params_shared, transform = TRUE)
  trans_vals <- trans_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    orig_vals_unit <- orig_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(orig_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(orig_vals_unit))
    
    coef(po, names(orig_vals_unit)) <- orig_vals_unit
    trans_vals_unit <- coef(po, names(orig_vals_unit), transform = TRUE)
    names(trans_vals_unit) <- paste0(seas[i], '_', names(orig_vals_unit))
    trans_vals <- c(trans_vals, trans_vals_unit)
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po, seas, params_all, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_shared <- trans_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared, transform = TRUE) <- trans_vals_shared
  orig_vals_shared <- coef(po, params_shared)
  orig_vals <- orig_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    trans_vals_unit <- trans_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(trans_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(trans_vals_unit))
    
    coef(po, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
    orig_vals_unit <- coef(po, names(trans_vals_unit))
    
    names(orig_vals_unit) <- paste0(seas[i], '_', names(trans_vals_unit))
    orig_vals <- c(orig_vals, orig_vals_unit)
  }
  
  return(orig_vals)
}

# ---------------------------------------------------------------------------------------------------------------------

# Create pomp model objects

# Loop through years and construct pomp models:
po_list <- vector('list', length(seasons))
for (yr_index in 1:length(seasons)) {
  yr <- seasons[yr_index]
  print(yr)
  
  # Load data and create pomp object:
  source('src/resp_interaction_model.R')
  
  # Check whether any data for given season:
  if (exists('resp_mod')) {
    
    # Add pomp object to list:
    po_list[[yr_index]] <- resp_mod
    
  }
  
  # Remove pomp object before repeating loop:
  rm(resp_mod)
  
}

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# Clean up:
rm(hk_dat, can_dat, us_dat, dat_pomp, age_structured, sens, yr_index, yr,
   nrow_check, debug_bool, unit_sp_estpars, vir1, vir2)

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# ---------------------------------------------------------------------------------------------------------------------

# Get results from maximum likelihood estimation

# Read in best-fit parameter values:
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')

# Get likelihoods:
ll <- apply(mle, 1, function(ix) {
  transform_params(as.numeric(ix), po_list[[1]], seasons, estpars, shared_estpars) %>%
    calculate_global_loglik()
})
mle <- mle %>%
  mutate(Llikelihood = ll)
rm(ll)

# Get MLEs and confidence intervals:
mle_ci <- read_csv('results/MLE_plus_95CI_from_bootstrapping_HPDI.csv')

# ---------------------------------------------------------------------------------------------------------------------

# Get results from MCMC

# Read in full fits for each set of starting parameter values and extract key info:
file_list <- list.files('results/mcmc/trans_scale/', full.names = TRUE)

res_list = map_list = vector('list', length(file_list))
for (i in 1:length(file_list)) {
  print(i)
  
  res_temp <- read_rds(file_list[i])
  
  res_list[[i]] <- getSample(sampler = res_temp, coda = TRUE, start = (n_iter / 10) * 9 + 1, thin = 10)
  map_list[[i]] <- MAP(res_temp)
  
  rm(res_temp)
}
rm(file_list, i)

# Get parameter values and likelihoods:
res <- lapply(map_list, function(ix) {
  back_transform_params(ix$parametersMAP %>% as.numeric(),
                        po_list[[1]],
                        seasons,
                        estpars,
                        shared_estpars)
}) %>%
  bind_rows() %>%
  mutate(jobid = 1:length(map_list), .before = rho1) %>%
  bind_cols(lapply(map_list, getElement, 'valuesMAP') %>%
              bind_rows())
rm(map_list)

# Get 95% confidence intervals (for each run):
res_ci_byRun <- lapply(res_list, function(ix) {
  lapply(ix, as_tibble) %>%
    bind_rows() %>%
    pivot_longer(everything()) %>%
    group_by(name) %>%
    summarise(mean = mean(value),
              median = median(value),
              lower = HPDI(value, p = 0.95)[1],
              upper = HPDI(value, p = 0.95)[2]) %>%
    ungroup()
}) %>%
  bind_rows(.id = 'run')

# Get 95% confidence intervals (overall):
res_ci <- lapply(res_list, function(ix) {
  lapply(ix, as_tibble) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  as.matrix() %>%
  apply(1, function(ix) {
    back_transform_params(ix, po_list[[1]], seasons, estpars, shared_estpars)
  }) %>%
  t() %>%
  as_tibble()

res_ci_UNIT <- res_ci %>%
  select(contains('I') | contains('R')) %>%
  select(-c(phi, rho1, rho2)) %>%
  mutate(iter = 1:nrow(res_ci)) %>%
  pivot_longer(-iter) %>%
  mutate(season = str_sub(name, 1, 6),
         name = str_sub(name, 8)) %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  mutate(`R10 + R120` = R10 + R120,
         `R20 + R120` = R20 + R120,
         R01 = Ri1 / (1 - (I10 + R10 + R120)),
         R02 = Ri2 / (1 - (I20 + R20 + R120))) %>%
  pivot_longer(Ri1:R02,
               names_to = 'parameter',
               values_to = 'value') %>%
  mutate(parameter = paste(season, parameter, sep = '_')) %>%
  select(parameter:value)
res_ci_SHARED <- res_ci %>%
  select(all_of(shared_estpars)) %>%
  mutate(delta2 = d2 * delta1) %>%
  pivot_longer(everything(),
               names_to = 'parameter',
               values_to = 'value') %>%
  select(parameter:value)

res_ci <- bind_rows(res_ci_SHARED, res_ci_UNIT) %>%
  group_by(parameter) %>%
  summarise(mean = mean(value),
            median = median(value),
            lower = HPDI(value, p = 0.95)[1],
            upper = HPDI(value, p = 0.95)[2]) %>%
  ungroup()
rm(res_ci_UNIT, res_ci_SHARED)

# # Check convergence:
# lapply(res_list, MCMCsummary, round = 4)

# Check a few trace plots:
set.seed(8942705)
to_check <- sample(1:length(res_list), size = 3, replace = FALSE)
for (i in to_check) {
  MCMCtrace(res_list[[i]], pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = n_iter / 10, params = shared_estpars)
}
rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Compare MCMC vs. MLE results

# Compare parameter estimates:
p1 <- res %>%
  mutate(analysis = 'mcmc') %>%
  bind_rows(mle %>%
              mutate(analysis = 'mle')) %>%
  pivot_longer(rho1:`s18-19_R120`) %>%
  mutate(name = factor(name, levels = estpars)) %>%
  filter(name %in% shared_estpars) %>%
  ggplot() +
  geom_violin(aes(x = analysis, y = value, group = analysis, fill = analysis), alpha = 0.6) +
  facet_wrap(~ name, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = NULL, y = NULL, title = 'Best-Fit Parameter Values (Shared)')
p2 <- res %>%
  mutate(analysis = 'mcmc') %>%
  bind_rows(mle %>%
              mutate(analysis = 'mle')) %>%
  pivot_longer(rho1:`s18-19_R120`) %>%
  mutate(name = factor(name, levels = estpars)) %>%
  filter(!(name %in% shared_estpars)) %>%
  ggplot() +
  geom_violin(aes(x = analysis, y = value, group = analysis, fill = analysis), alpha = 0.6) +
  facet_wrap(~ name, scales = 'free') +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = NULL, y = NULL, title = 'Best-Fit Parameter Values (Season-Specific)')

# Check colinearity between shared parameters:
p3 <- ggpairs(res %>% select(all_of(shared_estpars)),
              upper = list(continuous = wrap(ggally_cor, size = 3, method = 'kendall', digits = 2, display_grid = FALSE)),
              lower = list(continuous = wrap('points', size = 1.1)),
              labeller = 'label_parsed') +
  theme_classic() +
  theme(axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 55, vjust = 0.6),
        strip.text = element_text(size = 11),
        panel.border = element_rect(linewidth = 0.5, fill = NA))

# Check consistency of estimates across runs:
res_stdev <- res %>%
  mutate(analysis = 'mcmc') %>%
  bind_rows(
    mle %>%
      mutate(analysis = 'mle')
  ) %>%
  select(analysis, all_of(estpars)) %>%
  pivot_longer(-analysis) %>%
  group_by(analysis, name) %>%
  summarise(stdev = sd(value) / abs(mean(value))) %>%
  ungroup()
p4 <- ggplot(res_stdev %>% filter(name %in% shared_estpars),
             aes(x = analysis, y = stdev, color = analysis)) +
  geom_point(size = 3) +
  facet_wrap(~ name, scales = 'free') +
  # scale_y_continuous(limits = c(0, 0.05)) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = NULL, y = 'St. Dev. (Relative to Mean)')
p5 <- ggplot(res_stdev %>% filter(!(name %in% shared_estpars)),
             aes(x = analysis, y = stdev, color = analysis)) +
  geom_point(size = 3) +
  facet_wrap(~ name, scales = 'free', ncol = 7) +
  # scale_y_continuous(limits = c(0, 0.05)) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  theme(legend.position = 'none') +
  labs(x = NULL, y = 'St. Dev. (Relative to Mean)')

p6 <- ggplot(res_ci_byRun %>% filter(name %in% shared_estpars)) +
  geom_pointrange(aes(x = median, xmin = lower, xmax = upper, y = as.numeric(run))) +
  geom_point(aes(x = mean, y = as.numeric(run)), color = 'red2') +
  facet_wrap(~ name, scales = 'free_x', ncol = 3) +
  theme_classic()
p7 <- ggplot(res_ci_byRun %>% filter(!(name %in% shared_estpars))) +
  geom_pointrange(aes(x = median, xmin = lower, xmax = upper, y = as.numeric(run))) +
  geom_point(aes(x = mean, y = as.numeric(run)), color = 'red2') +
  facet_wrap(~ name, scales = 'free_x', ncol = 7) +
  theme_classic()

# Compare estimates and 95% CIs:
res_ci <- res_ci %>%
  select(-mean) %>%
  rename('mle' = 'median') %>%
  mutate(analysis = 'mcmc') %>%
  bind_rows(mle_ci %>%
              mutate(analysis = 'mle'))
res_ci_PLOT <- res_ci %>%
  filter(parameter %in% c(shared_estpars) | str_detect(parameter, 'Ri') | str_detect(parameter, 'I10') | str_detect(parameter, 'I20') |
           str_detect(parameter, 'R10 ') | str_detect(parameter, 'R20 ') | parameter == 'delta2') %>%
  filter(parameter != 'd2') %>%
  mutate(parameter = factor(parameter))

p8 <- ggplot(res_ci_PLOT %>% filter(parameter %in% c(shared_estpars, 'delta2'))) +
  geom_pointrange(aes(x = mle, xmin = lower, xmax = upper, y = analysis, color = analysis)) +
  facet_wrap(~ parameter, scales = 'free_x') +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_color_brewer(palette = 'Set1')
p9 <- ggplot(res_ci_PLOT %>% filter(!(parameter %in% shared_estpars) & parameter != 'delta2')) +
  geom_pointrange(aes(x = mle, xmin = lower, xmax = upper, y = analysis, color = analysis)) +
  facet_wrap(~ parameter, scales = 'free_x', ncol = 6) +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_color_brewer(palette = 'Set1')

# Compare log-likelihoods:
p10 <- res %>%
  mutate(analysis = 'mcmc') %>%
  bind_rows(mle %>%
              mutate(analysis = 'mle')) %>%
  select(analysis, Llikelihood, Lposterior) %>%
  pivot_longer(-analysis) %>%
  drop_na() %>%
  mutate(name = if_else(name == 'Llikelihood', 'LL', 'LPost'),
         name = paste(analysis, name, sep = '_'),
         name = factor(name, levels = c('mle_LL', 'mcmc_LL', 'mcmc_LPost'))) %>%
  ggplot() +
  geom_hline(yintercept = max(mle$Llikelihood) + qchisq(p = 0.95, df = length(estpars)), lty = 2) +
  geom_violin(aes(x = name, y = value, group = name, fill = name), alpha = 0.6) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_brewer(palette = 'Set2') +
  labs(x = NULL, y = 'Log-Likelihood')

# Output plots to file:
pdf('results/plots/check_fit_MCMC_trans.pdf', width = 14, height = 10)
for (i in to_check) {
  MCMCtrace(res_list[[i]], pdf = FALSE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, iter = n_iter / 10, params = shared_estpars)
}; rm(i)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)
print(p8)
print(p9)
print(p10)
dev.off()

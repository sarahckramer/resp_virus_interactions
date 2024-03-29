#######################################################################################################
#  Run simulations of the age-structured model of flu-RSV interaction
#######################################################################################################

rm(list = ls())

library(tidyverse)
library(testthat)
library(pomp)
library(socialmixr)
library(rootSolve)

source("src/age_structured_SA/f-CreateInteractionMod.R")

debug_bool <- FALSE
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))

# Functions

# Calculate reproduction number using next-generation matrix approach
Ri_fun <- function(beta, gamma, C_mat, tau_vec, r0_vec, N_vec) {
  NGM_mat <- matrix(data = 0, nrow = length(r0_vec), ncol = length(r0_vec))
  
  for(i in 1:nrow(NGM_mat)) {
    for(j in 1:ncol(NGM_mat)) {
      NGM_mat[i, j] <- beta * C_mat[i, j] / gamma * tau_vec[i] * (1 - r0_vec[i]) * N_vec[i] / N_vec[j]
    }
  }
  
  out <- NGM_mat %>% 
    eigen() %>% 
    getElement("values") %>% 
    max()
  return(out)
}

# Define function to find beta for a target Ri
find_beta <- function(beta, gamma, C_mat, tau_vec, r0_vec, N_vec, Ri) {
  Ri_fun(beta, gamma, C_mat, tau_vec, r0_vec, N_vec) - Ri
}

# Create contact matrix -------------------------------------------------------
# Contact data from Hong-Kong: https://zenodo.org/record/3874808
# 2021 mid-year population estimates from Hong-Kong: https://www.censtatd.gov.hk/en/web_table.html?id=1A
# Could not find pop <1 yr, assume it represents one fifth of the population aged 0-4 yr 
age_limits <- c(0, 1, 5, 16, 65) # Age limits for age categories
data(polymod)
CM_all <- contact_matrix(survey = get_survey(survey = "https://doi.org/10.5281/zenodo.3874808"), 
                         age.limits = age_limits,
                         survey.pop = data.frame(lower.age.limit = age_limits, population = c(45.8, 183.2, 578.8, 5153.8, 1451.5) * 1e3),
                         symmetric = T)

# Plot matrix
CM <- CM_all$matrix
rownames(CM) <- colnames(CM)
n_ages <- nrow(CM)
N <- CM_all$demography$population

# age1 represents the reporter age, age2 the contact age
CM_long <- CM %>%
  melt() %>%
  rename('age1' = 1,
         'age2' = 2,
         'contacts' = 3)

# Plot matrix; x-axis: contact age, y-axis: reporter age
# CAUTION: rates are per day!
pl <- ggplot(data = CM_long, 
             mapping = aes(x = age2, y = age1, fill = contacts)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  theme_minimal() + 
  #theme(legend.position = "top") + 
  labs(x = "Contact age", y = "Reported age", fill = "Daily contacts")
if (debug_bool) print(pl)

# Loop through seasons to generate synthetic data -----------------------------
seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
fit_params_names <- c('delta1', 'd2', 'theta_lambda1', 'theta_lambda2', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2', 'I10', 'I20', 'R10', 'R20', 'R120')

res_all_ages = res_combined = vector('list', length = length(seasons))

for (yr_index in 1:length(seasons)) {
  
  # Get season:
  yr <- seasons[yr_index]
  print(yr)
  
  # Get season-specific parameters --------------------------------------------
  mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')
  fit_params <- mle[1, ] %>%
    select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  # Flu and RSV parameters
  gamma_vec <- c(1 / 5, 1 / 10) # Recovery rates (per day)
  tau_l <- list(rep(1, n_ages), c(1, 0.75, 0.65, 0.65, 0.65)) # Susceptibility
  r0_l <- list(rep(unname(fit_params['R10'] + fit_params['R120']), n_ages),
               rep(unname(fit_params['R20'] + fit_params['R120']), n_ages)) # Initial fraction recovered
  R_vec <- c(c(unname(fit_params['Ri1']), unname(fit_params['Ri2']))) # Initial reproduction number
  beta_vec <- c(0, 0)
  
  for(i in 1:2) {
    
    # Find root of find_beta function
    tmp <- uniroot(f = find_beta, 
                   lower = 0, 
                   upper = 1, 
                   gamma = gamma_vec[i], 
                   C_mat = CM, 
                   tau_vec = tau_l[[i]],
                   r0_vec = r0_l[[i]], 
                   N_vec = N, 
                   Ri = R_vec[i])
    
    beta_vec[i] <- tmp$root
    print(sprintf("beta[%d] = %.2f", i, beta_vec[i]))
    
  }
  
  # Load and format necessary data --------------------------------------------
  hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv %>%
    filter(season == yr)
  
  dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
  
  hk_dat <- hk_dat %>%
    inner_join(dat_clim,
               by = c('Year' = 'year',
                      'Week' = 'week')) %>%
    select(time, n_T, GOPC, temp, ah) %>%
    rename('i_ILI' = 'GOPC')
  
  rm(dat_clim)
  
  # Create pomp object --------------------------------------------------------
  interMod <- CreateInteractionMod(dat = hk_dat, 
                                   nA = n_ages, 
                                   debug_bool = F)
  
  # Set parameters
  # CAUTION: rates of contact matrix are per day
  coef(interMod, paste0("N_", 1:n_ages)) <- N
  coef(interMod, paste0("CM_", 1:(n_ages * n_ages))) <- as.numeric(t(7 * CM))
  coef(interMod, c("b1", "b2")) <- beta_vec
  coef(interMod, fit_params_names) <- fit_params[fit_params_names]
  
  base_pars <- coef(interMod)
  init_pars <- base_pars[c("I10", "I20", "R10", "R20", "R120")]
  
  # Check initial conditions
  x0 <- rinit(interMod)
  x0 <- data.frame(agecat = rownames(x0), 
                   N0 = as.numeric(x0))
  
  if (debug_bool) {
    
    print(sprintf("Model tot pop: %d, data: %d", sum(x0$N0), sum(N)))
    
    init_vars <- c("X_IS", "X_SI", "X_RS", "X_SR", "X_RR")
    for(i in seq_along(init_vars)) {
      print(sprintf("Target initial conditions: %.6f", init_pars[i]))
      print("Realized initial conditions: ")
      print(x0$N0[x0$agecat %in% paste0(init_vars[i], "_", 1:n_ages)] / N)
    }
    
  }
  
  # Run simulation ------------------------------------------------------------
  coef(interMod, names(base_pars)) <- unname(base_pars)
  #coef(interMod, c("theta_lambda1", "theta_lambda2")) <- 1
  tj <- trajectory(object = interMod, 
                   format = "data.frame", 
                   ode_control = list(method = "lsoda"))
  
  # Add prevalence variables
  for(i in 1:n_ages) {
    tj[[paste0("p1_", i)]] <- rowSums(tj[, paste0(c("X_IS_", "X_II_", "X_IT_", "X_IR_"), i)])
    tj[[paste0("p2_", i)]] <- rowSums(tj[, paste0(c("X_SI_", "X_II_", "X_TI_", "X_RI_"), i)])
  }
  
  # Pivot to long format
  tj_long <- tj %>% 
    pivot_longer(cols = -c("time", ".id"), names_to = "var_nm_full") %>% 
    mutate(var_type = ifelse(str_detect(var_nm_full, "X"), "state_var", ifelse(str_detect(var_nm_full, "H"), "obs_var", "other_var")),
           var_nm = str_remove(string = var_nm_full, pattern = "_[0-9]"),
           age_no = str_extract(string = var_nm_full, pattern = "_[0-9]")) %>%
    mutate(age_no = str_remove(age_no, "_") %>% as.factor()) %>% 
    group_by(time, .id, age_no) %>% 
    mutate(N = sum(value[var_type == "state_var"])) %>% 
    ungroup()
  
  # Plot population size
  pl <- ggplot(data = tj_long %>% select(time, age_no, N) %>% unique(), 
               mapping = aes(x = time, y = N / 1e6, color = age_no)) + 
    geom_line() + theme_classic() +
    labs(x = "Time (weeks)", y = "Population size (millions)", color = "Age")
  if (debug_bool) print(pl)
  
  # Plot prevalence of infection
  var_to_plot <- c("H1", "H2")
  pl <- ggplot(data = tj_long %>% filter(var_nm %in% var_to_plot), 
               mapping = aes(x = time, y = value / N, color = age_no)) + 
    geom_line() + 
    facet_wrap(~ var_nm) +
    theme_classic() +
    labs(x = "Time (weeks)", y = "Fraction", color = "Age", title = var_to_plot)
  print(pl)
  
  # Generate synthetic observation data ---------------------------------------
  
  # Load and format age-structured covariate "data":
  hk_dat_covar <- read_csv('data/age_structured_SA/synthetic_covariate_data.csv') %>%
    filter(season == yr)
  
  hk_dat_covar <- hk_dat_covar %>%
    pivot_longer(n_T1:n_T5, names_to = 'age', values_to = 'n_T_age') %>%
    mutate(age = str_sub(age, 4, 4)) %>%
    arrange(time, age)
  
  # Set season-specific parameter values:
  rho1 <- fit_params['rho1']
  rho2 <- fit_params['rho2']
  alpha <- fit_params['alpha']
  phi <- fit_params['phi']
  
  # Loop through ages and generate synthetic data:
  res_by_age_list <- vector('list', length = n_ages)
  omega <- 2 * pi / 52.25
  for (i in 1:n_ages) {
    
    # Get age-specific covariate "data":
    dat_covar_temp <- hk_dat_covar %>%
      filter(age == i) %>%
      arrange(time)
    
    # Get age-specific simulated total cases:
    tj_temp <- tj %>%
      select(time,
             paste0('H1_', i),
             paste0('H2_', i)) %>%
      rename('H1' = paste0('H1_', i),
             'H2' = paste0('H2_', i))
    expect_equal(nrow(dat_covar_temp), nrow(tj_temp))
    
    # Calculate weekly probability of viral detection (assume same ILI attack rate in all age groups):
    rho1_w <- rho1 * (1.0 + alpha * cos(omega * (dat_covar_temp$time - phi))) * (tj_temp$H1 / dat_covar_temp$i_ILI)
    rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
    
    rho2_w <- rho2 * (1.0 + alpha * cos(omega * (dat_covar_temp$time - phi))) * (tj_temp$H2 / dat_covar_temp$i_ILI)
    rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
    
    # Generate synthetic observations (assume all age groups report):
    set.seed(1890435)
    
    obs1 <- rbinom(n = length(dat_covar_temp$n_T_age), size = dat_covar_temp$n_T_age, prob = rho1_w)
    obs2 <- rbinom(n = length(dat_covar_temp$n_T_age), size = dat_covar_temp$n_T_age, prob = rho2_w)
    
    # Compile synthetic observations into tibble:
    res_temp <- dat_covar_temp %>%
      select(time:age) %>%
      bind_cols(obs1, obs2) %>%
      rename('obs1' = '...8',
             'obs2' = '...9') %>%
      mutate(pop = N[i])
    
    # Store observations:
    res_by_age_list[[i]] <- res_temp
    
  }
  
  # Compile observations from all age groups and store:
  res_temp_all_ages <- bind_rows(res_by_age_list)
  res_all_ages[[yr_index]] <- res_temp_all_ages
  
  # Calculate total observations for all scenarios and store:
  res_temp_combined <- res_temp_all_ages %>%
    group_by(time, Year, Week, season) %>%
    summarise(across(obs1:obs2, sum),
              i_ILI = unique(i_ILI),
              n_T = unique(n_T),
              pop = sum(pop)) %>%
    ungroup() %>%
    select(time:season, i_ILI:pop, obs1:obs2)
  res_combined[[yr_index]] <- res_temp_combined
  
}

# Combine synthetic observations from all seasons:
res_all_ages <- bind_rows(res_all_ages)
res_combined <- bind_rows(res_combined)

# Reorder columns:
res_combined <- res_combined %>%
  select(time:Year, n_T, i_ILI, Week:season, pop, obs1:obs2)

# Visualize "observations" based on different sets of assumptions:
res_combined_long <- res_combined %>%
  pivot_longer(obs1:obs2,
               names_to = 'virus',
               values_to = 'val') %>%
  mutate(virus = if_else(virus == 'obs1', 'Influenza', 'RSV'))

pl <- ggplot(data = res_combined_long,
             aes(x = time, y = val, color = virus)) +
  geom_line() +
  facet_wrap(~ season, scale = 'free', ncol = 1) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(x = 'Time (Weeks)', y = '# of Cases', color = 'Virus')
print(pl)

# Compare to observed data, synthetic data from homogeneous mixing model ------
res_homogeneous <- NULL
for (yr_index in 1:length(seasons)) {
  
  # Get season:
  yr <- seasons[yr_index]
  print(yr)
  
  # Get season-specific parameters --------------------------------------------
  mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')
  fit_params <- mle[1, ] %>%
    select(rho1:eta_ah2, contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  
  # Create pomp object --------------------------------------------------------
  vir1 <- 'flu_h1'; vir2 <- 'rsv'
  Ri_max1 <- 3.0; Ri_max2 <- 2.0; d2_max <- 10.0
  sens <- 'main'
  fit_canada <- FALSE
  source('src/resp_interaction_model.R')
  
  # Set parameters
  coef(resp_mod, names(fit_params)) <- fit_params
  
  # Run simulation ------------------------------------------------------------
  sim <- simulate(object = resp_mod,
                  format = 'data.frame') %>%
    select(time, n_P1:n_P2)
  
  sim_long <- sim %>%
    pivot_longer(n_P1:n_P2,
                 names_to = 'virus',
                 values_to = 'synth_homogeneous') %>%
    mutate(season = yr,
           virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'))
  
  res_homogeneous <- bind_rows(res_homogeneous, sim_long)
  
}

hk_dat <- NULL
for (yr in seasons) {
  
  hk_dat_temp <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$h1_plus_b_rsv %>%
    filter(season == yr) %>%
    select(time, season, n_P1:n_P2)
  hk_dat <- bind_rows(hk_dat, hk_dat_temp)
  
}

hk_dat_long <- hk_dat %>%
  pivot_longer(n_P1:n_P2,
               names_to = 'virus',
               values_to = 'obs') %>%
  mutate(virus = if_else(virus == 'n_P1', 'Influenza', 'RSV'))

res_combined_long <- res_combined_long %>%
  rename('synth' = 'val')

res_combined_long <- res_combined_long %>%
  inner_join(hk_dat_long,
             by = c('time', 'season', 'virus')) %>%
  inner_join(res_homogeneous,
             by = c('time', 'season', 'virus'))

pl <- ggplot(data = res_combined_long,
             aes(x = time, color = virus)) +
  geom_point(aes(y = obs)) +
  geom_line(aes(y = synth)) +
  geom_line(aes(y = synth_homogeneous), lty = 2) +
  facet_wrap(~ season, scale = 'free') +
  theme_classic() +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Time (Weeks)', y = '# of Cases', color = 'Virus')
print(pl)

# Write to file:
write_csv(res_all_ages, 'data/age_structured_SA/synthetic_obs_by_age.csv')
write_csv(res_combined, 'data/age_structured_SA/synthetic_obs_combined.csv')

#######################################################################################################

# Clean up:
rm(list = ls())

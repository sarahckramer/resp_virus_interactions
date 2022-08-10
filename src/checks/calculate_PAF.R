# ---------------------------------------------------------------------------------------------------------------------
# Calculate population attributable fraction (PAF) of the interaction effect at the MLE
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)

# Get names of fitted parameters:
shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                    'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

# Set parameter values necessary for loading models:
prof_lik <- FALSE
lag_val <- 0

# ---------------------------------------------------------------------------------------------------------------------

# Generate synthetic data at the MLE, with and without interaction

# Read in MLEs:
mle_h1 <- read_rds('results/MLEs_flu_h1.rds')
mle_b <- read_rds('results/MLEs_flu_b.rds')

# Get synthetic data (H1N1/RSV):
vir1 <- 'flu_h1'
source('src/functions/setup_global_likelilhood.R')

traj_list_H1N1 <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Get season:
  yr <- seasons[i]
  
  # Get pomp object:
  resp_mod <- po_list[[i]]
  
  # Get best-fit model parameters:
  pars_temp <- mle_h1[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  # Generate matrix of model parameters for fitting:
  param_mat_temp <- parmat(coef(resp_mod), nrep = 5)
  
  # Remove impact of flu from some parameter sets:
  param_mat_temp['theta_lambda1', 2] <- 1.0
  param_mat_temp['I10', 3] <- 0
  
  # Remove impact of RSV from some parameter sets:
  param_mat_temp['theta_lambda2', 4] <- 1.0
  param_mat_temp['I20', 5] <- 0
  
  # Simulate using deterministic model:
  traj_temp <- trajectory(resp_mod, params = param_mat_temp, format = 'data.frame') %>%
    mutate(season = yr) %>%
    select(time:season, H1:H2)
  
  # Check that removal of interaction works as expected:
  expect_true(all.equal(traj_temp %>% filter(.id == 2) %>% pull(H2), traj_temp %>% filter(.id == 3) %>% pull(H2)))
  expect_true(all.equal(traj_temp %>% filter(.id == 4) %>% pull(H1), traj_temp %>% filter(.id == 5) %>% pull(H1)))
  
  # Remove unneeded simulations:
  traj_temp <- traj_temp %>%
    filter(.id %in% c(1, 3, 5))
  
  # Store:
  traj_list_H1N1[[i]] <- traj_temp
  
}

# Get synthetic data (B/RSV):
vir1 <- 'flu_b'
source('src/functions/setup_global_likelilhood.R')

traj_list_B <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Get season:
  yr <- seasons[i]
  
  # Get pomp object:
  resp_mod <- po_list[[i]]
  
  # Get best-fit model parameters:
  pars_temp <- mle_b[1, ] %>%
    select(all_of(shared_estpars),
           contains(yr))
  names(pars_temp)[(length(names(pars_temp)) - 6):length(names(pars_temp))] <- unit_estpars
  
  coef(resp_mod, true_estpars) <- pars_temp
  
  # Generate matrix of model parameters for fitting:
  param_mat_temp <- parmat(coef(resp_mod), nrep = 5)
  
  # Remove impact of flu from some parameter sets:
  param_mat_temp['theta_lambda1', 2] <- 1.0
  param_mat_temp['I10', 3] <- 0
  
  # Remove impact of RSV from some parameter sets:
  param_mat_temp['theta_lambda2', 4] <- 1.0
  param_mat_temp['I20', 5] <- 0
  
  # Simulate using deterministic model:
  traj_temp <- trajectory(resp_mod, params = param_mat_temp, format = 'data.frame') %>%
    mutate(season = yr) %>%
    select(time:season, H1:H2)
  
  # Check that removal of interaction works as expected:
  expect_true(all.equal(traj_temp %>% filter(.id == 2) %>% pull(H2), traj_temp %>% filter(.id == 3) %>% pull(H2)))
  expect_true(all.equal(traj_temp %>% filter(.id == 4) %>% pull(H1), traj_temp %>% filter(.id == 5) %>% pull(H1)))
  
  # Remove unneeded simulations:
  traj_temp <- traj_temp %>%
    filter(.id %in% c(1, 3, 5))
  
  # Store:
  traj_list_B[[i]] <- traj_temp
  
}

# Compile:
res_H1N1 <- bind_rows(traj_list_H1N1) %>%
  mutate(virus_pair = 'Influenza H1N1 / RSV') %>%
  as_tibble()

res_B <- bind_rows(traj_list_B) %>%
  mutate(virus_pair = 'Influenza B / RSV') %>%
  as_tibble()

res <- bind_rows(res_H1N1, res_B)
rm(res_H1N1, res_B)

# ---------------------------------------------------------------------------------------------------------------------

# Calculate statistics and plot

# Calculate PAF (or negative of PAF, i.e., % increase in AR expected if interaction were not present):
res_ars <- res %>%
  group_by(virus_pair, season, .id) %>%
  summarise(attack_rate_H1 = sum(H1),
            attack_rate_H2 = sum(H2))

paf_h1_fluOnRSV = paf_h1_RSVonFlu = paf_b_fluOnRSV = c()
for (yr in unique(res_ars$season)) {
  
  res_temp_H1N1 <- res_ars %>%
    filter(season == yr & virus_pair == 'Influenza H1N1 / RSV')
  res_temp_B <- res_ars %>%
    filter(season == yr & virus_pair == 'Influenza B / RSV')
  
  if (nrow(res_temp_H1N1) > 0) {
    
    # Impact of flu on RSV:
    ar_pop <- res_temp_H1N1$attack_rate_H2[res_temp_H1N1$.id == 1]
    ar_unexposed <- res_temp_H1N1$attack_rate_H2[res_temp_H1N1$.id == 3]
    
    # paf_h1_fluOnRSV <- c(paf_h1_fluOnRSV, (ar_pop - ar_unexposed) / ar_pop) # negative
    paf_h1_fluOnRSV <- c(paf_h1_fluOnRSV, (ar_unexposed - ar_pop) / ar_pop)
    
    # Impact of RSV on flu:
    ar_pop <- res_temp_H1N1$attack_rate_H1[res_temp_H1N1$.id == 1]
    ar_unexposed <- res_temp_H1N1$attack_rate_H1[res_temp_H1N1$.id == 5]
    
    paf_h1_RSVonFlu <- c(paf_h1_RSVonFlu, (ar_unexposed - ar_pop) / ar_pop)
    
  }
  
  if (nrow(res_temp_B) > 0) {
    
    # Impact of flu on RSV:
    ar_pop <- res_temp_B$attack_rate_H2[res_temp_H1N1$.id == 1]
    ar_unexposed <- res_temp_B$attack_rate_H2[res_temp_H1N1$.id == 3]
    
    paf_b_fluOnRSV <- c(paf_b_fluOnRSV, (ar_unexposed - ar_pop) / ar_pop)
    
  }
  
}

print(summary(paf_h1_fluOnRSV))
print(summary(paf_b_fluOnRSV))

print(summary(paf_h1_RSVonFlu))

# Plot how interaction influences outbreak dynamics:
res <- res %>%
  mutate(virus_pair = factor(virus_pair, levels = c('Influenza H1N1 / RSV', 'Influenza B / RSV')))
# res <- res %>%
#   mutate(Week = time + 45) %>%
#   mutate(Week = if_else(Week > 52 & season != 's16-17', Week - 52, Week),
#          Week = if_else(Week > 53 & season == 's16-17', Week - 53, Week))

p1 <- ggplot() +
  geom_line(data = res %>%
              filter(.id %in% c(1, 3)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H2, col = .id)) +
  geom_line(data = res %>%
              filter(.id == 1),
            aes(x = time, y = H1), lty = 2) +
  facet_grid(virus_pair ~ season, scales = 'free_x') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_continuous(breaks = seq(0, 53, by = 10)) +
  scale_color_manual(values = c('#3182bd', '#9ecae1')) +
  labs(x = 'Time (Weeks)', y = 'RSV Incidence', color = '', tag = 'A')

p2 <- ggplot() +
  geom_line(data = res %>%
              filter(virus_pair == 'Influenza H1N1 / RSV') %>%
              filter(.id %in% c(1, 5)) %>%
              mutate(.id = if_else(.id == 1, 'Interaction', 'No Interaction')),
            aes(x = time, y = H1, col = .id)) +
  geom_line(data = res %>%
              filter(virus_pair == 'Influenza H1N1 / RSV') %>%
              filter(.id == 1),
            aes(x = time, y = H2), lty = 2) +
  facet_grid(virus_pair ~ season, scales = 'free_x') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom',
        plot.tag = element_text(size = 22),
        plot.tag.position = c(0.01, 0.98)) +
  scale_x_continuous(breaks = seq(0, 53, by = 10)) +
  scale_color_manual(values = c('#de2d26', '#fc9272')) +
  labs(x = 'Time (Weeks)', y = 'Influenza Incidence', color = '', tag = 'B')

fig <- arrangeGrob(p1, p2, ncol = 1, heights = c(1.6, 1))
plot(fig)

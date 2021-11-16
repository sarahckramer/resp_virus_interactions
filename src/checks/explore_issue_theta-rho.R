# ---------------------------------------------------------------------------------------------------------------------
# Determine why changing values of theta_rho1 appear to have no influence on simulations
# ---------------------------------------------------------------------------------------------------------------------

# Set seed:
set.seed(1098302745)

# Load libraries:
library(tidyverse)
library(gridExtra)

# ---------------------------------------------------------------------------------------------------------------------

# Set range of interaction parameter values:
int_vals <- seq(0, 1, by = 0.1)

# Set number of simulations:
n_sim <- 5

# Set model parameters:
vir1 <- 'flu_h1' # 'flu_A', 'flu_B'
vir2 <- 'rsv'
yr <- 's15-16' #c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

debug_bool <- FALSE

Ri_max1 <- 2.0
Ri_max2 <- 3.0
delta_min <- 7 / 30.0

# Load pomp object:
source('src/resp_interaction_model.R')

# ---------------------------------------------------------------------------------------------------------------------

# Select ~5-10 sets of "realistic" parameter values
pars_df <- read_csv('results/round1_interaction/res_traj_match_round1.csv')
pars_df <- pars_df %>%
  filter(virus1 == vir1,
         year == yr) %>%
  arrange(desc(loglik)) %>%
  filter(rho1 < 0.99,
         rho2 < 0.99)
which_params <- c(1, sample(2:100, 4)) # 1 96 28 22 68
pars_df <- pars_df[which_params, ] %>%
  select(Ri1:rho2)

# Add ~5 sets of additional parameter values, to explore
pars_check <- cbind(c(1.3, 1.4, 1.5, 1.3, 1.3),
                    c(1.3, 1.3, 1.3, 1.1, 1.4),
                    rep(3.1e-5, 5),
                    rep(4.5e-4, 5),
                    rep(0.1, 5),
                    c(0.1, 0.1, 0.1, 0.1, 0.1),
                    c(0.4, 0.4, 0.4, 0.4, 0.4),
                    c(0.1, 0.1, 0.1, 0.1, 0.1),
                    c(0.15, 0.15, 0.15, 0.15, 0.15)) %>%
  as.data.frame()
names(pars_check) <- names(pars_df)
pars_df <- rbind(pars_df, pars_check)

# Set up parameter matrix:
p_mat <- parmat(params = coef(resp_mod), nrep = nrow(pars_df))
p_mat[names(pars_df), ] <- t(pars_df)

# ---------------------------------------------------------------------------------------------------------------------

# How does theta_lambda behave with these parameter sets?:
res_list_det = res_list_stoch = vector('list', length = length(int_vals))
for (i in 1:length(int_vals)) {
  p_mat_temp <- p_mat
  p_mat_temp['theta_lambda1', ] <- int_vals[i]
  
  sim_determ <- trajectory(object = resp_mod,
                           params = p_mat_temp,
                           format = 'data.frame') %>%
    select(time:.id, H1_tot:H2) %>%
    mutate(theta_lambda1 = int_vals[i]) %>%
    as_tibble()
  sim_stoch <- simulate(object = resp_mod,
                        params = p_mat_temp,
                        nsim = n_sim,
                        format = 'data.frame') %>%
    select(time:.id, H1_tot:n_P2) %>%
    # select(time:.id, n_P1:n_P2) %>%
    mutate(theta_lambda1 = int_vals[i],
           id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
           id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
    as_tibble()
  
  res_list_det[[i]] <- sim_determ
  res_list_stoch[[i]] <- sim_stoch
}
rm(i, p_mat_temp, sim_determ, sim_stoch)

sim_determ <- bind_rows(res_list_det) %>%
  pivot_longer(H1_tot:H2, names_to = 'state', values_to = 'val') %>%
  mutate(state_type = if_else(str_detect(state, 'tot'), 'tot', 'red'))
sim_stoch <- bind_rows(res_list_stoch) %>%
  pivot_longer(H1_tot:n_P2, names_to = 'state', values_to = 'val')

p1 <- ggplot(data = sim_determ,
             aes(x = time, y = val, col = state, lty = state_type)) +
  geom_line() + theme_classic() +
  facet_grid(theta_lambda1 ~ .id) +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'State Var.') +
  scale_linetype_discrete(guide = 'none')
# print(p1) # check that H1_tot = H1, H2_tot = H2

check_det <- bind_rows(res_list_det)
expect_true(all.equal(check_det$H1, check_det$H1_tot))
expect_true(all.equal(check_det$H2, check_det$H2_tot))

p2 <- ggplot(data = sim_determ %>% filter(state_type == 'red'),
             aes(x = time, y = val, col = theta_lambda1, group = theta_lambda1)) +
  geom_line() + theme_classic() +
  # facet_grid(.id ~ state, scales = 'free') +
  facet_grid(state ~ .id, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_lambda1') +
  scale_color_viridis()
# print(p2) # check influence of theta_lambda1 on total cases eligible for reporting

# p2 <- ggplot(data = sim_determ %>% filter(state_type != 'red'),
#              aes(x = time, y = val, col = theta_lambda1, group = theta_lambda1)) +
#   geom_line() + theme_classic() +
#   # facet_grid(.id ~ state, scales = 'free') +
#   facet_grid(state ~ .id, scales = 'free') +
#   labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_lambda1') +
#   scale_color_viridis()
# print(p2) # should be same as above

p3 <- ggplot(data = sim_stoch %>% filter(state %in% c('H1', 'H2')),
             aes(x = time, y = val, col = theta_lambda1, group = paste(theta_lambda1, id.sim))) +
  geom_line() + theme_classic() +
  facet_grid(state ~ id.parm, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_lambda1') +
  scale_color_viridis()
# print(p3) # ensure that stochastic simulations show same pattern
grid.arrange(p2, p3, ncol = 1)

check_stoch <- bind_rows(res_list_stoch)
expect_true(all.equal(check_stoch$H1, check_stoch$H1_tot))
expect_true(all.equal(check_stoch$H2, check_stoch$H2_tot))

p4 <- ggplot(data = sim_stoch %>% filter(state %in% c('n_P1', 'n_P2')),
             aes(x = time, y = val, col = theta_lambda1, group = paste(theta_lambda1, id.sim))) +
  geom_line() + theme_classic() +
  facet_grid(state ~ id.parm, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Pos. Tests', col = 'theta_lambda1') +
  scale_color_viridis()
print(p4) # check results for simulated cases (positive tests)

# Store plots:
plots_theta_lambda1 <- list(p1, p2, p3, p4)

# ---------------------------------------------------------------------------------------------------------------------

# Now check behavior of theta_rho:
res_list_det = res_list_stoch = vector('list', length = length(int_vals))
for (i in 1:length(int_vals)) {
  p_mat_temp <- p_mat
  p_mat_temp['theta_rho1', ] <- int_vals[i]
  
  sim_determ <- trajectory(object = resp_mod,
                           params = p_mat_temp,
                           format = 'data.frame') %>%
    select(time:.id, H1_tot:H2) %>%
    mutate(theta_rho1 = int_vals[i]) %>%
    as_tibble()
  sim_stoch <- simulate(object = resp_mod,
                        params = p_mat_temp,
                        nsim = n_sim,
                        format = 'data.frame') %>%
    select(time:.id, H1_tot:n_P2) %>%
    # select(time:.id, n_P1:n_P2) %>%
    mutate(theta_rho1 = int_vals[i],
           id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
           id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
    as_tibble()
  
  res_list_det[[i]] <- sim_determ
  res_list_stoch[[i]] <- sim_stoch
}
rm(i, p_mat_temp, sim_determ, sim_stoch)

sim_determ <- bind_rows(res_list_det) %>%
  pivot_longer(H1_tot:H2, names_to = 'state', values_to = 'val') %>%
  mutate(state_type = if_else(str_detect(state, 'tot'), 'tot', 'red'))
sim_stoch <- bind_rows(res_list_stoch) %>%
  pivot_longer(H1_tot:n_P2, names_to = 'state', values_to = 'val')

p1 <- ggplot(data = sim_determ,
             aes(x = time, y = val, col = state, lty = state_type)) +
  geom_line() + theme_classic() +
  facet_grid(theta_rho1 ~ .id) +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'State Var.') +
  scale_linetype_discrete(guide = 'none')
# print(p1) # check that H1_tot = H1, H2_tot != H2

check_det <- bind_rows(res_list_det)
expect_true(all.equal(check_det$H1, check_det$H1_tot))
expect_false(all.equal(check_det$H2, check_det$H2_tot) == TRUE)

p2 <- ggplot(data = sim_determ %>% filter(state_type == 'red'),
             aes(x = time, y = val, col = theta_rho1, group = theta_rho1)) +
  geom_line() + theme_classic() +
  # facet_grid(.id ~ state, scales = 'free') +
  facet_grid(state ~ .id, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_rho1') +
  scale_color_viridis()
# print(p2) # check influence of theta_rho1 on total cases eligible for reporting

# p2 <- ggplot(data = sim_determ %>% filter(state_type != 'red'),
#              aes(x = time, y = val, col = theta_rho1, group = theta_rho1)) +
#   geom_line() + theme_classic() +
#   # facet_grid(.id ~ state, scales = 'free') +
#   facet_grid(state ~ .id, scales = 'free') +
#   labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_rho1') +
#   scale_color_viridis()
# print(p2) # should show no difference

p3 <- ggplot(data = sim_stoch %>% filter(state %in% c('H1', 'H2')),
             aes(x = time, y = val, col = theta_rho1, group = paste(theta_rho1, id.sim))) +
  geom_line() + theme_classic() +
  facet_grid(state ~ id.parm, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Cases', col = 'theta_rho1') +
  scale_color_viridis()
# print(p3) # ensure that stochastic simulations show same pattern
grid.arrange(p2, p3, ncol = 1)

check_stoch <- bind_rows(res_list_stoch)
expect_true(all.equal(check_stoch$H1, check_stoch$H1_tot))
expect_false(all.equal(check_stoch$H2, check_stoch$H2_tot) == TRUE)

p4 <- ggplot(data = sim_stoch %>% filter(state %in% c('n_P1', 'n_P2')),
             aes(x = time, y = val, col = theta_rho1, group = paste(theta_rho1, id.sim))) +
  geom_line() + theme_classic() +
  facet_grid(state ~ id.parm, scales = 'free') +
  labs(x = 'Time (Weeks)', y = 'Pos. Tests', col = 'theta_rho1') +
  scale_color_viridis()
print(p4) # check results for simulated cases (positive tests)

# Store plots:
plots_theta_rho1 <- list(p1, p2, p3, p4)

# ---------------------------------------------------------------------------------------------------------------------

# Save plots to file:
pdf('results/plots/check_thetarho.pdf', width = 15, height = 8.5)
print(plots_theta_lambda1)
print(plots_theta_rho1)
dev.off()

# Clean up:
rm(list = ls())

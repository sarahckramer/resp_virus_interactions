# ---------------------------------------------------------------------------------------------------------------------
# Code to explore the influence of individual parameters on model dynamics
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
# library(tidyverse)
# library(pomp)
# library(tgp)
# library(gridExtra)

# Set global parameters:
n_sim <- 5

# ---------------------------------------------------------------------------------------------------------------------

# Prep pomp object

# Set necessary parameter values:
vir1 <- 'flu_A' # 'flu_A', 'flu_B'
vir2 <- 'rsv'
yr <- 2006

debug_bool <- FALSE

Ri_max1 <- 2.0
Ri_max2 <- 2.0
delta_min <- 7 / 60.0

# Load pomp object:
source('src/resp_interaction_model.R')

# ---------------------------------------------------------------------------------------------------------------------

# First let's briefly explore the impact of all parameters:

coef(resp_mod, c('Ri1', 'Ri2', 'I10', 'I20', 'R20', 'R120')) <- c(1.3, 1.8, 0.00005, 0.002, 0.45, 0.45)

slices_a <- slice_design(center = coef(resp_mod),
                         Ri1 = seq(1.0, 1.5, by = 0.1),
                         Ri2 = seq(1.0, 2.0, by = 0.2),
                         I10 = seq(0, 0.0001, by = 0.00002),
                         I20 = seq(0, 0.004, by = 0.0005),
                         rho1 = c(0.2, 0.5, 0.75),
                         rho2 = c(0.15, 0.35, 0.5),
                         # R10 = seq(0, 1, by = 0.1),
                         R20 = seq(0, 0.5, by = 0.1),
                         R120 = seq(0, 0.5, by = 0.1))
slices_df <- slices_a %>% as_tibble() %>% mutate(id.parm = 1:nrow(slices_a))

parms_mat <- parmat(params = coef(resp_mod), nrep = nrow(slices_a))
for (i in 1:nrow(slices_a)) {
  parms_mat[, i] <- unlist(slices_a[i, 1:(ncol(slices_a) - 1)])
}
rm(i)

parms_df <- parms_mat %>% t() %>% as_tibble() %>% mutate(id.parm = 1:ncol(parms_mat))

sim_dat <- simulate(object = resp_mod,
                    params = parms_mat,
                    nsim = n_sim,
                    format = 'data.frame')

sim_dat <- sim_dat %>%
  as_tibble() %>%
  select(time:.id, H1_tot:n_P2) %>%
  mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
         id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
  inner_join(slices_df) %>%
  select(-c(gamma1:theta_lambda2, theta_rho1:N))

sim_dat <- sim_dat %>%
  pivot_longer(Ri1:R120, names_to = 'param', values_to = 'param_val') %>%
  filter(slice == param)

p.ri1 <- ggplot(data = sim_dat %>% filter(slice == 'Ri1')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'Ri1')
p.ri2 <- ggplot(data = sim_dat %>% filter(slice == 'Ri2')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'Ri2')
p.I10 <- ggplot(data = sim_dat %>% filter(slice == 'I10')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'I10')
p.I20 <- ggplot(data = sim_dat %>% filter(slice == 'I20')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'I20')
# p.R10 <- ggplot(data = sim_dat %>% filter(slice == 'R10')) +
#   geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
#   geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
#   facet_wrap(~ param_val, ncol = 1) +
#   theme_classic() +
#   labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'R10')
p.R20 <- ggplot(data = sim_dat %>% filter(slice == 'R20')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'R20')
p.R120 <- ggplot(data = sim_dat %>% filter(slice == 'R120')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'R120')
p.rho1 <- ggplot(data = sim_dat %>% filter(slice == 'rho1')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'rho1')
p.rho2 <- ggplot(data = sim_dat %>% filter(slice == 'rho2')) +
  geom_line(aes(x = time, y = n_P1, group = id.sim), col = 'coral') +
  geom_line(aes(x = time, y = n_P2, group = id.sim), col = 'steelblue2') +
  facet_wrap(~ param_val, ncol = 1) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'Reported Cases', title = 'rho2')

p_list <- list(p.ri1, p.ri2, p.I10, p.I20, p.R20, p.R120, p.rho1, p.rho2)
print(p_list)
# Ri1: 1.3 seems to produce most realistic pattern; clear difference with 0.1 change in value
# Ri2: Outbreaks small regardless of value, but this could change with reporting rate; using 1.6, 1.8, 2.0
# I10: Larger values lead to earlier outbreaks, but only slighly; can use 0.00005 maybe
# I20: Like with Ri2, higher values lead to larger and earlier outbreaks, but all are small; use around 0.002
# R10: Need to check this one later; with high R20/R120, can't go above 0.1
# R20: Higher values lead to suppressed and shorter outbreaks
# R120: Higher values suppress and shorten outbreaks of both viruses; peaks also tend to occur earlier
# rho1: Higher values lead to more cases being reported; so outbreak shape is the same, but all higher
# rho2: Same as above

# ---------------------------------------------------------------------------------------------------------------------

# Now let's see how changing R10/R20/R120 while holding other parameters to various realistic values impacts dynamics:
# Or else for several sets of the above parameter sets

R_vals <- seq(0, 1, by = 0.1)
R_df <- NULL
for (ix in R_vals) {
  for (jx in R_vals) {
    for (kx in R_vals) {
      
      if (ix + jx + kx < 1.0) {
        R_df <- rbind(R_df, c(ix, jx, kx))
      }
      
    }
  }
}
rm(ix, jx, kx)
R_df <- R_df %>%
  as_tibble() %>%
  rename('R10' = 'V1',
         'R20' = 'V2',
         'R120' = 'V3') %>%
  mutate(id.parm = 1:nrow(R_df))

slices_b <- slice_design(center = coef(resp_mod),
                         Ri2 = c(1.8, 1.6, 2.0))#,
                         # rho1 = 0.25,
                         # rho2 = c(0.15, 0.35))

sim_list <- vector('list', nrow(slices_b))
for (i in 1:nrow(slices_b)) {
  
  parms_mat <- parmat(params = unlist(slices_b[i, 1:(ncol(slices_b) - 1)]), nrep = nrow(R_df))
  for (j in 1:nrow(R_df)) {
    parms_mat[c('R10', 'R20', 'R120'), j] <- unlist(R_df[j, c('R10', 'R20', 'R120')])
  }
  rm(j)
  
  # parms_mat['theta_lambda1', ] <- 0
  # # parms_mat['delta', ] <- 7 / 30
  # parms_mat['rho2', ] <- 0.15
  
  sim_dat <- simulate(object = resp_mod,
                      params = parms_mat,
                      nsim = n_sim,
                      format = 'data.frame')
  
  sim_dat <- sim_dat %>%
    as_tibble() %>%
    select(time:.id, H1_tot:n_P2) %>%
    mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
           id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
    inner_join(R_df)
  
  # sim_dat <- sim_dat %>%
  #   pivot_longer(Ri1:R120, names_to = 'param', values_to = 'param_val') %>%
  #   filter(slice == param)
  
  sim_list[[i]] <- sim_dat
}
rm(i)

sim_names <- c('Default', 'Low Ri2', 'High Ri2')#, 'Rho1=0.25', 'Rho2=0.15', 'Rho2=0.35')
for (i in 1:length(sim_list)) {
  sim_dat <- sim_list[[i]]
  
  p1 <- ggplot(data = sim_dat, aes(group = paste(id.sim, R120, R20))) +
    geom_line(aes(x = time, y = n_P1, col = R120)) +
    facet_wrap(~ R10) + theme_classic() + scale_color_viridis() +
    labs(x = 'Time', y = 'Flu Cases', title = sim_names[i])
  p2 <- ggplot(data = sim_dat, aes(group = paste(id.sim, R120, R10))) +
    geom_line(aes(x = time, y = n_P2, col = R120)) +
    facet_wrap(~ R20) + theme_classic() + scale_color_viridis() +
    labs(x = 'Time', y = 'RSV Cases', title = sim_names[i])
  grid.arrange(p1, p2)
}
rm(i)

# ---------------------------------------------------------------------------------------------------------------------

# Start with high Ri2/R20/R120 and lower all proportionally - how does this change dynamics?

parms_mat <- parmat(params = coef(resp_mod), nrep = 5)
for (i in 2:ncol(parms_mat)) {
  parms_mat[c('Ri2', 'R20', 'R120'), i] <- 0.9 * parms_mat[c('Ri2', 'R20', 'R120'), i - 1]
}
parms_df <- parms_mat %>% t() %>% as_tibble() %>% mutate(id.parm = 1:5)

sim_dat <- simulate(object = resp_mod,
                    params = parms_mat,
                    nsim = n_sim,
                    format = 'data.frame')

sim_dat <- sim_dat %>%
  as_tibble() %>%
  select(time:.id, H1_tot:n_P2) %>%
  mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
         id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
  inner_join(parms_df) %>%
  select(time:id.sim, Ri2, R20, R120)

ggplot(data = sim_dat, aes(group = id.sim)) +
  geom_line(aes(x = time, y = n_P1), col = 'coral') +
  geom_line(aes(x = time, y = n_P2), col = 'steelblue2') +
  facet_wrap(~ id.parm, ncol = 1) +
  theme_classic() +
  labs(x = 'Time', y = 'Cases')
# Reducing both Ri2 and R20/R120 seems to lead to a later start and longer duration

# ---------------------------------------------------------------------------------------------------------------------

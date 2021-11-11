# ---------------------------------------------------------------------------------------------------------------------
# Explore causes of poor fit to H1 outbreaks in 2017-18 and 2018-19
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(tidyverse)



# Set necessary parameters:
vir2 <- 'rsv'
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

# ---------------------------------------------------------------------------------------------------------------------

# Load and format results

# Set parameters estimated:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2', 'theta_lambda1', 'delta')

# Read in results:
pars_df <- read_csv('results/round1_interaction/res_traj_match_round1.csv')
pars_df_1314 <- read_csv('results/round1_interaction_leadNAs/res_traj_match_round1.csv')

# Compile:
pars_df <- pars_df_1314 %>%
  bind_rows(pars_df %>%
              filter(year != 's13-14'))
rm(pars_df_1314)

# ---------------------------------------------------------------------------------------------------------------------

# Simulate and plot

# Save to pdf:
pdf('results/plots/issue_with_H1_simulations.pdf',
    width = 14, height = 9)

# Loop through viruses/seasons:
for (vir1 in c('flu_h1', 'flu_b')) {
  pars_vir <- pars_df %>%
    filter(virus1 == vir1)
  
  plot_list <- vector('list', length = 4 * length(unique(pars_vir$year)))
  
  for (yr in unique(pars_vir$year)) {
    print(yr)
    
    # Subset and format data:
    pars_temp <- pars_vir %>%
      filter(year == yr) %>%
      select(-c(year, virus1))
    
    # Order by ll:
    pars_temp <- pars_temp %>%
      arrange(desc(loglik))
    
    # # Look at LL and simulations for all fits:
    # summary(pars_temp$loglik)
    # plot(pars_temp$loglik, pch = 20)
    
    # Load pomp object:
    source('src/resp_interaction_model.R')
    
    # Simulate with top parameter set:
    coef(resp_mod, estpars) <- pars_temp[1, 1:length(estpars)]
    rho1_const <- coef(resp_mod, 'rho1')
    rho2_const <- coef(resp_mod, 'rho2')
    
    sim_temp <- simulate(resp_mod, nsim = 5, format = 'data.frame') %>%
      select(time:.id, H1:n_P2) %>%
      arrange(.id) %>%
      inner_join(dat_pomp %>%
                   select(time, n_T:i_ARI) %>%
                   rename('obs1' = 'n_P1',
                          'obs2' = 'n_P2'),
                 by = 'time') %>%
      mutate(rho1_w = rho1_const * (H1 / i_ARI),
             rho2_w = rho2_const * (H2 / i_ARI)) %>%
      mutate(ILI_test_rat = i_ARI / n_T,
             ILI_test_rat_norm = ILI_test_rat / min(ILI_test_rat, na.rm = TRUE)) %>%
      as_tibble()
    
    # Create and store plots:
    p1 <- ggplot(data = sim_temp) + #geom_line(aes(x = time, y = i_ARI, group = .id)) +
      geom_line(aes(x = time, y = H1, group = .id), col = 'steelblue2') +
      geom_line(aes(x = time, y = H2, group = .id), col = 'coral') +
      theme_classic() + labs(x = 'Time', y = 'ILI Incidence', title = yr)
    
    p2 <- ggplot(data = sim_temp) + #geom_line(aes(x = time, y = i_ARI, group = .id)) +
      geom_line(aes(x = time, y = rho1_w, group = .id), col = 'steelblue2') +
      geom_line(aes(x = time, y = rho2_w, group = .id), col = 'coral') +
      theme_classic() + labs(x = 'Time', y = 'Prob. of Detection', title = yr)
    
    # p3 <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = ILI_test_rat_norm, group = .id)) +
    #   theme_classic() + labs(x = 'Time', y = 'Ratio of ILI:n_T', title = yr)# + scale_y_continuous(limits = c(1.0, 5.0))
    # p3 <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = i_ARI / n_T * 100000, group = .id)) +
    #   geom_line(aes(x = time, y = rho1_w, group = .id), col = 'steelblue2') +
    #   theme_classic() + labs(x = 'Time', y = 'Ratio of i_ILI:n_T / Prob. of Detection', title = yr)
    # p3 <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = rho1_w / i_ARI, group = .id)) +
    #   theme_classic() + labs(x = 'Time', y = 'Ratio of P(Detect):i_ILI', title = yr)
    p3 <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = (rho1_w)/n_T, group = .id)) +
      theme_classic() + labs(x = 'Time', y = 'Ratio of P(Detect):n_T', title = yr)
    
    p4 <- ggplot(data = sim_temp) + geom_line(aes(x = time, y = n_P1, group = .id), col = 'black') +
      geom_line(aes(x = time, y = n_P2, group = .id), col = 'coral') + 
      geom_point(aes(x = time, y = obs1, group = .id)) + geom_point(aes(x = time, y = obs2, group = .id), col = 'coral') +
      theme_classic() +
      labs(x = 'Time', y = '# Positive Tests', title = yr)
    
    plot_list[[4 * which(unique(pars_vir$year) == yr) - 3]] <- p1
    plot_list[[4 * which(unique(pars_vir$year) == yr) - 2]] <- p2
    plot_list[[4 * which(unique(pars_vir$year) == yr) - 1]] <- p3
    plot_list[[4 * which(unique(pars_vir$year) == yr)]] <- p4
    
  }
  
  # Print plots:
  do.call('grid.arrange', c(plot_list, ncol = 4))
}

# Close plots:
dev.off()

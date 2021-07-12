# ---------------------------------------------------------------------------------------------------------------------
# Code to test model setup and results for expected behavior
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(testthat)

# Functions:

check_transformations <- function(pomp_object) {
  # Function to check that parameters are correctly being transformed
  # params pomp_object: The pomp model object to be checked
  
  obj <- pomp_object
  x <- coef(obj, transform = TRUE)
  obj1 <- obj
  coef(obj1, transform = TRUE) <- x
  
  expect_true(all.equal(coef(obj), coef(obj1)),
              info = 'Parameters not correctly transformed')
}

check_params <- function(pomp_object) {
  # Function to check that setting parameter values works as expected
  # params pomp_object: The pomp model object to be checked
  
  select_parms <- sample(1:length(coef(pomp_object)), 3)
  coef(pomp_object, names(coef(pomp_object))[select_parms]) <- c(0.1, 0.2, 0.3)
  expect_true(all(coef(pomp_object)[select_parms] == c(0.1, 0.2, 0.3)))
}

check_correct_N_CONST <- function(pomp_object, true_n, n_sim = 10) {
  # Function to check that deterministic and stochastic simulations maintain correct N (when pop size constant)
  # params pomp_object: The pomp model object to be checked
  # params true_n: The expected population size
  # params n_sim: The number of stochastic simulations to run
  sim_determ <- trajectory(object = pomp_object, format = 'data.frame') %>%
    rowwise() %>%
    mutate(Ncheck = sum(c_across(X_SS:X_RR)),
           Ntrue = true_n)
  # sim_stoch <- simulation(object = pomp_object, nsim = n_sim, format = 'data.frame') %>%
  #   rowwise() %>%
  #   mutate(Ncheck = sum(c_across(X_SS:X_RR)),
  #          Ntrue = true_n)
  
  expect_true(all.equal(sim_determ$Ncheck, sim_determ$Ntrue))
  # expect_true(all.equal(sim_stoch$Ncheck, sim_stoch$Ntrue))
}

# check_obs_lessthan_samples <- function(pomp_object, n_sim = 10) {
#   # Function to check that simulated case numbers never exceed the total number of tests performed
#   # params pomp_object: The pomp model object to be checked
#   # params n_sim: The number of stochastic simulations to run
#   # returns: Plot of # of samples, observations, and simulated positive tests
#   
#   sim_stoch <- simulate(object = pomp_object, nsim = n_sim, format = 'data.frame')
#   sim_stoch <- sim_stoch %>%
#     as_tibble() %>%
#     left_join(dat_pomp, by = 'time') %>%
#     rename(sim_obs = obs.x,
#            obs = obs.y)
#   
#   expect_true(all(sim_stoch$sim_obs <= sim_stoch$n_samp))
#   
#   p_temp <- ggplot(data = sim_stoch) + geom_line(aes(x = time, y = n_samp), lwd = 1.5) +
#     geom_line(aes(x = time, y = sim_obs, group = .id, col = .id)) +
#     geom_point(aes(x = time, y = obs)) +
#     labs(x = 'Time (Weeks)', y = 'Simulated Observations') + theme_classic() +
#     scale_color_viridis(discrete = TRUE)
#   
#   return(p_temp)
# }

check_independent_dynamics <- function(pomp_object) {
  # Function to check that, when no interaction is specified, virus dynamics are independent
  # params pomp_object: The pomp model object to be checked
  
  p_mat <- parmat(params = coef(pomp_object), nrep = 3)
  p_mat['I10', ] <- c(1e-5, 1e-5, 0)
  p_mat['I20', ] <- c(1e-5, 0, 1e-5)
  
  sim_determ <- trajectory(object = pomp_object,
                           params = p_mat,
                           format = 'data.frame') %>%
    pivot_longer(cols = -c('time', '.id'), names_to = 'var_nm', values_to = 'val') %>%
    mutate(var_type = if_else(str_detect(var_nm, 'X_'), 'state', 'accum')) %>%
    filter(var_nm %in% c('H1', 'H2'))
  
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H1'],
                        sim_determ$val[sim_determ$.id == 2 & sim_determ$var_nm == 'H1']))
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H2'],
                        sim_determ$val[sim_determ$.id == 3 & sim_determ$var_nm == 'H2']))
  
  p_temp <- ggplot(data = sim_determ %>% filter(.id == 1), aes(x = time, y = 100 * val)) +
    geom_point(aes(color = var_nm)) +
    geom_line(data = sim_determ %>% filter(.id == 2, var_nm == 'H1'), color = 'pink') +
    geom_line(data = sim_determ %>% filter(.id == 3, var_nm == 'H2'), color = 'purple') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)') +
    theme_classic()
  
  return(p_temp)
}

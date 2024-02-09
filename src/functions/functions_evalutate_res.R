# ---------------------------------------------------------------------------------------------------------------------
# Functions to help evaluate model results
# ---------------------------------------------------------------------------------------------------------------------

run_sim <- function(pomp_object, seas, mle, shared_estpars, unit_estpars, model_type, obs_only, analysis) {
  # Function to generate deterministic or stochastic simulations from the model at the MLE
  # params pomp_object: The pomp model object used to run the model
  # params seas: The season to be simulated
  # params mle: The maximum likelihood estimates of all parameter values
  # params shared_estpars: A vector containing the names of all shared parameters
  # params unit_estpars: A vector containing the names of all season-specific parameters
  # params model_type: Should simulations be taken from the deterministic or the stochastic model?
  # params obs_only: Should average number of observed cases also be included? (model_type deterministic only)
  # params analysis: Should any special analyses be performed? (options: 'paf')
  # returns: Tibble of simulated values for flu and RSV in a given season
  
  # Get all estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get parameter values:
  pars <- mle[1, ] %>%
    select(all_of(shared_estpars),
           contains(seas))
  names(pars)[(length(names(pars)) - 6):length(names(pars))] <- unit_estpars
  
  coef(pomp_object, true_estpars) <- pars
  
  # Simple simulation or additional analysis?:
  if (analysis == 'paf') {
    
    if (model_type == 'deterministic') {
      
      if (!obs_only) {
        
        # Generate matrix of model parameters for fitting:
        param_mat_temp <- parmat(coef(pomp_object), nrep = 5)
        
        # Remove impact of flu from some parameter sets:
        param_mat_temp['theta_lambda1', 2] <- 1.0
        param_mat_temp['I10', 3] <- 0
        
        # Remove impact of RSV from some parameter sets:
        param_mat_temp['theta_lambda2', 4] <- 1.0
        param_mat_temp['I20', 5] <- 0
        
        # Simulate using deterministic model:
        traj_temp <- trajectory(pomp_object, params = param_mat_temp, format = 'data.frame') %>%
          mutate(season = yr) %>%
          select(time:season, H1:H2)
        
        # Check that removal of interaction works as expected:
        expect_true(all.equal(traj_temp %>% filter(.id == 2) %>% pull(H2), traj_temp %>% filter(.id == 3) %>% pull(H2)))
        expect_true(all.equal(traj_temp %>% filter(.id == 4) %>% pull(H1), traj_temp %>% filter(.id == 5) %>% pull(H1)))
        
        # Remove unneeded simulations:
        traj_temp <- traj_temp %>%
          filter(.id %in% c(1, 3, 5))
        
      }
    }
    
  } else {
    
    # Run simple simulation
    
    # Deterministic or stochastic?:
    if (model_type == 'deterministic') {
      
      # Get trajectory at MLE:
      traj_temp <- trajectory(pomp_object, format = 'data.frame') %>%
        mutate(season = seas) %>%
        select(time, season, H1:H2) %>%
        cbind(pomp_object@covar@table[1, ]) %>%
        cbind(pomp_object@covar@table[2, ]) %>%
        as_tibble()
      names(traj_temp)[5:6] <- c('i_ILI', 'n_T')
      
      if (obs_only) {
        
        # Calculate mean number of observed cases:
        rho1 <- as.numeric(pars['rho1'])
        rho2 <- as.numeric(pars['rho2'])
        alpha <- as.numeric(pars['alpha'])
        phi <- as.numeric(pars['phi'])
        
        rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H1 / traj_temp$i_ILI
        rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (traj_temp$time - phi))) * traj_temp$H2 / traj_temp$i_ILI
        
        rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
        rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
        
        expect_equal(nrow(traj_temp), length(rho1_w))
        expect_equal(nrow(traj_temp), length(rho2_w))
        
        traj_temp$rho1_w <- rho1_w
        traj_temp$rho2_w <- rho2_w
        
        traj_temp <- traj_temp %>%
          mutate(obs1 = rho1_w * n_T,
                 obs2 = rho2_w * n_T) %>%
          select(time:H2, obs1:obs2)
        
      }
      
    } else if (model_type == 'stochastic') {
      
      
      
    } else {
      
      stop('Model type neither deterministic nor stochastic!')
      
    }
    
  }
  
  # Return simulations:
  return(traj_temp)
  
}


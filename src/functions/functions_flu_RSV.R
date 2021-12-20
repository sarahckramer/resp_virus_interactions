# ---------------------------------------------------------------------------------------------------------------------
# Functions to assist with flu/RSV model
# ---------------------------------------------------------------------------------------------------------------------

prepare_data <- function(virus1_nm, virus2_nm, epiyear_val, dat, early_start = FALSE) {
  # Function to format French flu/RSV data for input into pomp object
  # param virus1_nm: Name of first virus of interest (string, 'flu_A', 'flu_B', or 'rsv')
  # param virus2_nm: Name of second virus of interest (string, 'flu_A', 'flu_B', or 'rsv')
  # param epiyear_val: season (numeric, for example 2010 for season 2009/2010)
  # param dat: Unformatted data (data frame or tibble)
  # param early_start: Start model run before data are available?
  # returns: List containing formatted data and data formatted specifically for pomp
  
  out <- vector(mode = 'list', length = 2)
  names(out) <- c('dat_full', 'dat_pomp')
  
  # Prepare data for all viruses during the given year:
  out[['dat_full']] <- dat %>%
    mutate(virus = fct_collapse(.f = type,
                                'rsv' = 'rsv',
                                'flu_B' = 'flu_b',
                                'flu_A' = c('flu_ah1', 'flu_ah3', 'flu_a_nst'))) %>%
    filter(epiyear == epiyear_val) %>%
    group_by(week_no, week_date, type, virus) %>%
    summarise(pop_tot = sum(pop),
              pop_eff = sum(pop[!is.na(ira_n)]),
              n_samp = sum(n_samp, na.rm = TRUE),
              n_pos = sum(n_pos, na.rm = TRUE),
              n_ARI = sum(ira_n, na.rm = TRUE)) %>%
    group_by(week_no, week_date, virus) %>%
    summarise(pop = unique(pop_tot),
              pop_eff = unique(pop_eff),
              n_ARI = unique(n_ARI),
              n_samp = unique(n_samp),
              n_pos = sum(n_pos)) %>%
    ungroup() %>%
    mutate(i_ARI = n_ARI / pop_eff,
           i_ARI_wrong = n_ARI / pop,
           time = as.numeric(week_date - min(week_date)) / 7 + 1) %>%
    arrange(week_date)
  
  # Format data to pass to pomp object:
  out[['dat_pomp']] <- out[['dat_full']] %>%
    filter(virus %in% c(virus1_nm, virus2_nm)) %>% 
    pivot_wider(names_from = 'virus', values_from = 'n_pos') %>%
    rename(c('n_P1' = all_of(virus1_nm), 'n_P2' = all_of(virus2_nm), 'n_T' = 'n_samp')) %>%
    arrange(time)
  
  stopifnot(all(out[['dat_pomp']]$n_P1 >= 0))
  stopifnot(all(out[['dat_pomp']]$n_P2 >= 0))
  stopifnot(all(out[['dat_pomp']]$n_T >= 0))
  stopifnot(all(out[['dat_pomp']]$i_ARI >= 0))
  
  # Add option to start run before week 40:
  if (early_start) {
    desired_start_wk <- 30
    n_rep <- length(desired_start_wk:39)
    
    na_df <- matrix(data = NA, nrow = n_rep, ncol = ncol(out[['dat_pomp']])) %>%
      as_tibble() %>%
      mutate(V1 = 39:desired_start_wk,
             V2 = min(out[['dat_pomp']]$week_date) - seq(7, by = 7, length.out = n_rep),
             V3 = unique(out[['dat_pomp']]$pop))
    names(na_df) <- names(out[['dat_pomp']])
    
    out[['dat_pomp']] <- out[['dat_pomp']] %>%
      bind_rows(na_df) %>%
      arrange(week_date) %>%
      mutate(time = as.numeric(week_date - min(week_date)) / 7 + 1)
  }
  
  return(out)
}

create_SITRxSITR_mod <- function(dat, Ri1_max = 3.0, Ri2_max = 3.0, delta_min = 7 / 60, debug_bool = F) {
  #Args: 
  # dat: ARI and virological data (data frame, first time point must be 1) 
  # Ri1_max: upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # Ri2_max: upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # delta_min: lower bound of 1 / refractory period duration (denominator is upper bound of duration)
  # (double, passed as global argument in the C script) 
  # debug_bool: should debugging info be printed? (boolean)
  
  # Read model C code:
  mod_code <- readLines('src/resp_interaction_model.c')
  components_nm <- c('globs', 'toest', 'fromest', 'dmeas', 'rmeas', 'rinit', 'skel', 'rsim')
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for (nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
      str_flatten(collapse = '\n')
    
    if(nm == 'globs') {
      components_l[[nm]] <- paste(components_l[[nm]], 
                                  sprintf('static int debug = %d; \nstatic double Ri1_max = %f; \nstatic double Ri2_max = %f; \nstatic double delta_min = %f;', 
                                          as.integer(debug_bool),
                                          as.numeric(Ri1_max),
                                          as.numeric(Ri2_max),
                                          as.numeric(delta_min)), 
                                  sep = '\n')
    }
    
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  # Create pomp object:
  po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
             times = 'time',
             t0 = 0,
             covar = covariate_table(dat[, c('time', 'i_ARI', 'n_T', 'temp', 'ah')], times = 'time'),
             accumvars = c('H1_tot', 'H2_tot', 'H1', 'H2'),
             obsnames = c('n_P1', 'n_P2'),
             statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                            'X_SI', 'X_II', 'X_TI', 'X_RI', 
                            'X_ST', 'X_IT', 'X_TT', 'X_RT',
                            'X_SR', 'X_IR', 'X_TR', 'X_RR', 
                            'H1_tot', 'H2_tot', 
                            'H1', 'H2'),
             paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                            'gamma1', 'gamma2', # 1 / average infectious periods
                            # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                            'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                            'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                            'rho1', 'rho2', # probs. infection leads to ARI consultation
                            'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                            'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                            'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                            'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                            'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                            'N', # population size
                            'I10', 'I20', # props. infectious at outbreak start
                            'R10', 'R20', 'R120'), # props. recovered at outbreak start
             params = c(Ri1 = 1.5, Ri2 = 2,
                        gamma1 = 7 / 5, gamma2 = 7 / 10, # or 4 for flu?
                        # delta = 7 / 5,
                        delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                        theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                        rho1 = 0.5, rho2 = 0.15,
                        alpha = 0, phi = 0,
                        theta_rho1 = 1.0, theta_rho2 = 1.0,
                        eta_temp1 = 0, eta_temp2 = 0,
                        eta_ah1 = 0, eta_ah2 = 0,
                        beta_sd1 = 0.1, beta_sd2 = 0.1,
                        N = unique(dat$pop),
                        I10 = 1e-5, I20 = 1e-5,
                        R10 = 0, R20 = 0, R120 = 0),
             globals = components_l[['globs']],
             dmeasure = components_l[['dmeas']],
             rmeasure = components_l[['rmeas']],
             skeleton = vectorfield(components_l[['skel']]),
             rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
             partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
             rinit = components_l[['rinit']]
  )
  
  return(po)
}

#######################################################################################################
# Create pomp object for the age-structured model of flu-RSV interaction
# Extension on Sarah's homogeneous model, with other parameters taken from Waterlow et al. Epidemics 2021
# All rates are PER WEEK, as in Sarah's homogeneous model
#######################################################################################################

CreateInteractionMod <- function(dat, nA = 5, debug_bool = F) {
  # Args: 
  # dat: data frame with temperature and AH data
  # nA: no of age classes
  # dt: time step (in days) for stochastic model (currently not used)
  # debug_bool: boolean, should messages be displayed to help debug? 
  
  # Extract C code from file
  mod_code <- readLines("src/age_structured_SA/c-interactionMod.c")
  components_nm <- c("globs", "skel", "rinit")
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for(nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0("start_", nm)):str_which(mod_code, paste0("end_", nm))] %>% 
      str_flatten(collapse = "\n")
    
    if(nm == "globs") {
      components_l[[nm]] <- paste(components_l[[nm]], 
                                  sprintf("static int nA = %d;\nstatic int debug = %d;", 
                                          nA, as.integer(debug_bool)), 
                                  sep = "\n")
    }
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  po <- pomp(data = data.frame(time = seq(from = 1, to = 52, by = 1), X = NA),
             times = "time",
             t0 = 0,
             covar = covariate_table(dat[, c('time', 'temp', 'ah')], times = 'time'),
             obsnames = "X",
             statenames = c(
               paste0("X_SS_", 1:nA),
               paste0("X_IS_", 1:nA),
               paste0("X_TS_", 1:nA),
               paste0("X_RS_", 1:nA),
               paste0("X_SI_", 1:nA),
               paste0("X_II_", 1:nA),
               paste0("X_TI_", 1:nA),
               paste0("X_RI_", 1:nA),
               paste0("X_ST_", 1:nA),
               paste0("X_IT_", 1:nA),
               paste0("X_TT_", 1:nA),
               paste0("X_RT_", 1:nA),
               paste0("X_SR_", 1:nA),
               paste0("X_IR_", 1:nA),
               paste0("X_TR_", 1:nA),
               paste0("X_RR_", 1:nA),
               paste0("H1tot_", 1:nA),
               paste0("H2tot_", 1:nA),
               paste0("H1_", 1:nA), 
               paste0("H2_", 1:nA)
             ),
             accumvars = c(
               paste0("H1tot_", 1:nA),
               paste0("H2tot_", 1:nA),
               paste0("H1_", 1:nA), 
               paste0("H2_", 1:nA)
             ),
             paramnames = c(paste0("CM_", 1:(nA * nA)), # Contact matrix
                            paste0("N_", 1:nA), # Age-specific populations
                            paste0("tau_", 1:nA), # Age-specific susceptibility to RSV infection
                            'b1', 'b2', # transmission coefficients
                            'gamma1', 'gamma2', # 1 / average infectious periods
                            'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                            'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                            'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                            'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                            'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                            'I10', 'I20', # props. infectious at outbreak start
                            'R10', 'R20', 'R120'), # props. recovered at outbreak start
             params = c(setNames(object = rep(0, nA * nA), nm = paste0("CM_", 1:(nA * nA))),
                        setNames(object = rep(0, nA), paste0("N_", 1:nA)),
                        setNames(object = c(1, 0.75, 0.65, 0.65, 0.65), paste0("tau_", 1:nA)),
                        b1 = 6.152, b2 = 3.823,
                        gamma1 = 7 / 5, gamma2 = 7 / 10, # or 4 for flu?
                        delta1 = 0.065, d2 = 2.1, 
                        theta_lambda1 = 0, theta_lambda2 = 0,
                        theta_rho1 = 1.0, theta_rho2 = 1.0,
                        eta_temp1 = -0.288, eta_temp2 = -0.287,
                        eta_ah1 = 0.456, eta_ah2 = 0.342,
                        I10 = 6.386e-05, 
                        I20 = 1.598e-05,
                        R10 = 0.0507, 
                        R20 = 0.108, 
                        R120 = 0.539),
             globals = components_l[["globs"]],
             skeleton = vectorfield(components_l[["skel"]]),
             rinit = components_l[["rinit"]]
             # cdir = "_help/", 
             # cfile = "codes"
  )
  return(po)
}



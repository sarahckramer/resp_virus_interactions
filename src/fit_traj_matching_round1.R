# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# Round 1: Fit each season individually (w/ or w/o interaction) to obtain good start values for Ri/initial conditions
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(nloptr)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
sens <- as.character(Sys.getenv("SENS")); print(sens)
fit_canada <- as.logical(Sys.getenv("FITCANADA")); print(fit_canada)

if (fit_canada) {
  yr <- c('s10-11', 's11-12', 's12-13', 's13-14')[(ceiling(jobid / no_jobs) - 1) %% 4 + 1]; print(yr)
  vir1 <- 'flu'
} else {
  yr <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')[(ceiling(jobid / no_jobs) - 1) %% 6 + 1]; print(yr)
  vir1 <- 'flu_h1_plus_b'
}
jobid <- (jobid - 1) %% no_jobs + 1; print(jobid)

# # Set parameters for local runs (temp):
# jobid <- 1
# no_jobs <- 10
# 
# vir1 <- 'flu_h1_plus_b'
# yr <- 's15-16'
# 
# sobol_size <- 500
# search_type <- 'broad'
# sens <- 'main'
# fit_canada <- FALSE

# Set parameters for run:
vir2 <- 'rsv'

time_max <- 1.75 # 5.0 # Maximal execution time (in hours)
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
d2_max <- 10.0

# # Fit for synthetic data from age-structured model?:
# age_structured <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Run trajectory matching
# Note: Separate runs for each season

# Load data and create pomp object:
source('src/resp_interaction_model.R')

# Check that sufficient epidemic activity:
if (exists('resp_mod')) {
  
  # Set start ranges for estimated parameters:
  if (sens == 'no_rsv_immune') {
    estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'rho1', 'rho2')
  } else {
    estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')
  }
  
  start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                            Ri2 = c(1.0, Ri_max2),
                            I10 = c(0, 1e-3),
                            I20 = c(0, 1e-3),
                            R10 = c(0, 0.3),
                            R20 = c(0, 0.3),
                            R120 = c(0, 0.3),
                            rho1 = c(0, 1.0),
                            rho2 = c(0, 1.0),
                            theta_lambda1 = c(0, 1.0),
                            theta_lambda2 = c(0, 1.0),
                            theta_rho1 = c(0, 1.0),
                            theta_rho2 = c(0, 1.0),
                            delta1 = c(7 / 60, 7),
                            d2 = c(7 / 60, 7),
                            alpha = c(0, 0.5),
                            phi = c(0, 52.25),
                            eta_temp1 = c(-0.5, 0.5),
                            eta_temp2 = c(-0.5, 0.5),
                            eta_ah1 = c(-0.5, 0.5),
                            eta_ah2 = c(-0.5, 0.5))
  start_range <- start_range[, estpars]
  
  if (search_type == 'broad') {
    start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                                 upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                                 nseq = sobol_size)
  }
  
  print(estpars)
  print(start_range)
  print(summary(start_values))
  
  # Create objective function for call to nloptr:
  obj_fun <- traj_objfun(data = resp_mod,
                         est = estpars,
                         partrans = resp_mod@partrans,
                         verbose = TRUE)
  
  # Set maximal execution time for each estimation:
  nmins_exec <- time_max * 60 / (sobol_size / no_jobs)
  print(sprintf("Max estimation time=%.1f min", nmins_exec))
  
  # Get unique identifiers:
  sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)
  
  # Loop through start values and perform trajectory matching:
  for (i in seq_along(sub_start)) {
    
    print(paste0('Estimation: ', sub_start[i]))
    
    # Get param start values:
    x0 <- as.numeric(start_values[sub_start[i], ])
    coef(resp_mod, estpars) <- x0
    x0_trans <- coef(resp_mod, estpars, transform = TRUE)
    print(-1 * obj_fun(x0_trans))
    
    # Run trajectory matching using subplex algorithm:
    # http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
    tic <- Sys.time()
    m <- try(
      nloptr(x0 = unname(x0_trans),
             eval_f = obj_fun,
             opts = list(algorithm = 'NLOPT_LN_SBPLX',
                         maxtime = 60.0 * nmins_exec,
                         maxeval = -1, # disabled
                         xtol_rel = -1, # disabled; default: 1e-4
                         print_level = 0))
    )
    toc <- Sys.time()
    etime <- toc - tic
    units(etime) <- 'mins'
    print(etime)
    
    # If estimation is successful, save results:
    if (!inherits(m, 'try-error')) {
      coef(resp_mod, estpars, transform = TRUE) <- m$solution
      
      # Collect all results:
      out <- list(allpars = coef(resp_mod),
                  estpars = coef(resp_mod, estpars),
                  ll = -m$objective,
                  conv = m$status,
                  message = m$message,
                  niter = m$iterations,
                  etime = as.numeric(etime))
      
      # Write to file:
      saveRDS(out,
              file = sprintf('results/res_%s_%s_%s_%d.rds',
                             vir1, vir2,
                             as.character(yr),
                             sub_start[i])
      )
      
      # Print results:
      print(out$ll)
      print(out$estpars, digits = 2)
      print(out$conv)
      print(out$message)
    }
    
  }
  
} else {
  print('No data for given season.')
}

print('Done.')

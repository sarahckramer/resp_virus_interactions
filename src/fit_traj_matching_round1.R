# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# Round 1: Fit each season without interaction to obtain good start values for Ri/initial conditions
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

# yr <- c(2006:2014)[(ceiling(jobid / no_jobs) - 1) %% 9 + 1]; print(yr)
# vir1 <- c('flu_A', 'flu_B')[ceiling(jobid / (no_jobs * 9))]; print(vir1)
yr <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')[(ceiling(jobid / no_jobs) - 1) %% 6 + 1]; print(yr)
vir1 <- c('flu_h1', 'flu_b')[ceiling(jobid / (no_jobs * 6))]; print(vir1)
jobid <- (jobid - 1) %% no_jobs + 1; print(jobid)

# # Set parameters for local runs (temp):
# jobid <- 1
# no_jobs <- 10
# 
# vir1 <- 'flu_b' # 'flu_A', 'flu_B'
# yr <- 's15-16'
# 
# sobol_size <- 500
# search_type <- 'broad'

# Set parameters for run:
vir2 <- 'rsv'

time_max <- 9.75 # Maximal execution time (in hours)
debug_bool <- FALSE

Ri_max1 <- 3.0
Ri_max2 <- 3.0
delta_min <- 7 / 60.0

# ---------------------------------------------------------------------------------------------------------------------

# Run trajectory matching
# Note: Separate runs for each season

# Load data and create pomp object:
source('src/resp_interaction_model.R')

# Check that sufficient epidemic activity:
if (exists('resp_mod')) {
  if (sum(resp_mod@data[1, ], na.rm = TRUE) > 100 & sum(resp_mod@data[2, ], na.rm = TRUE) > 100) {
    
    # Set start ranges for estimated parameters:
    # estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')
    estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2', 'theta_lambda1', 'delta')
    
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
                              delta = c(7 / 60, 7))
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
                           maxtime = 60.0,
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
    print('Insufficient epidemic activity')
  }
} else {
  print('No data for given season.')
}

print('Done.')

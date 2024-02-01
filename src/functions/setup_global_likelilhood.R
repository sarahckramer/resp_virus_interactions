# ---------------------------------------------------------------------------------------------------------------------
# Code to read in and prepare all pomp objects needed to quickly calculate the global log-likelihood
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(nloptr)

# Set parameters for run:
debug_bool <- FALSE
vir2 <- 'rsv'

seasons_hk <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')
seasons_can <- c('s10-11', 's11-12', 's12-13', 's13-14')

Ri_max1 <- 2.0
Ri_max2 <- 3.0
d2_max <- 10.0

if (!exists('sens')) {
  sens <- 'main'
}

if (sens == 'less_circ_h3') {
  seasons_hk <- c('s17-18', 's18-19')
}

# ---------------------------------------------------------------------------------------------------------------------

# Functions

create_obj_fxn <- function(po, estpars) {
  # Creates the objective function for a given season
  # param po: A pomp model object for a specific season
  # param estpars: A vector listing the parameters to be fit
  # returns: The negative log likelihood
  
  ofun <- traj_objfun(data = po, 
                      est = estpars, 
                      partrans = po@partrans,
                      verbose=TRUE)
  
  return(ofun)
}

calculate_global_loglik <- function(trans_vals) {
  # Calculates the log-likelihood for each season, and combines to yield a global log-likelihood
  # param trans_vals: Unnamed vector of transformed parameters; fxn only works if this is the only input?
  # returns: The global, negative log-likelihood
  
  # Add names to vector:
  if(is.null(names(trans_vals))) names(trans_vals) <- x0_trans_names
  
  # Split params into global, shared, and unit params:
  shared_hk_in <- as.data.frame(sapply(paste0('^', 'hk_'), grep, x = names(trans_vals)))
  unit_hk_in <- as.data.frame(sapply(paste0('^', seasons_hk, '_hk_'), grep, x = names(trans_vals)))
  names(unit_hk_in) <- seasons_hk
  
  shared_can_in <- as.data.frame(sapply(paste0('^', 'can_'), grep, x = names(trans_vals)))
  unit_can_in <- as.data.frame(sapply(paste0('^', seasons_can, '_can_'), grep, x = names(trans_vals)))
  names(unit_can_in) <- seasons_can
  
  global_params <- trans_vals[-c(unlist(shared_hk_in), unlist(shared_can_in), unique(unlist(unit_hk_in)), unique(unlist(unit_can_in)))]
  
  shared_hk_params <- trans_vals[unlist(shared_hk_in)]
  names(shared_hk_params) <- gsub('hk_', '', names(shared_hk_params))
  shared_hk_params <- c(global_params, shared_hk_params)
  
  shared_can_params <- trans_vals[unlist(shared_can_in)]
  names(shared_can_params) <- gsub('can_', '', names(shared_can_params))
  shared_can_params <- c(global_params, shared_can_params)
  
  unit_hk_params <- list()
  for (i in 1:length(seasons_hk)) {
    unit <- trans_vals[unit_hk_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_hk_params[[i]] <- unit
    }
  }
  
  unit_can_params <- list()
  for (i in 1:length(seasons_can)) {
    unit <- trans_vals[unit_can_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 3]
      unit_can_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season (HK):
  units_hk_ll <- rep(NA, length(seasons_hk))
  
  for (i in 1:length(seasons_hk)) {
    if(length(unit_hk_params) > 0) {
      params_temp <- c(shared_hk_params, unit_hk_params[[i]])
    } else {
      params_temp <- c(shared_hk_params)
    }
    
    units_hk_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Get -ll for each season (Canada):
  units_can_ll <- rep(NA, length(seasons_can))
  
  for (i in 1:length(seasons_can)) {
    if (length(unit_can_params) > 0) {
      params_temp <- c(shared_can_params, unit_can_params[[i]])
    } else {
      params_temp <- c(shared_can_params)
    }
    
    units_can_ll[i] <- obj_fun_list[[i + length(seasons_hk)]](params_temp)
  }
  
  # Calculate -ll for each location:
  hk_ll <- sum(units_hk_ll)
  can_ll <- sum(units_can_ll)
  
  # Calculate global -ll:
  glob_ll <- hk_ll + can_ll
  return(glob_ll)
}

transform_params <- function(orig_vals, po1, po2, seas1, seas2, params_all, params_global, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po1: pomp object with correct parameter transformations (Hong Kong)
  # param po2: pomp object with correct parameter transformations (Canada)
  # param seas1: Vector containing all seasons of interest (Hong Kong)
  # params seas2: Vector containing all seasons of interest (Canada)
  # param params_all: Names of all parameters to be estimated
  # param params_global: Names of all global parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_global <- orig_vals[which(params_global %in% params_all)]
  
  coef(po1, params_global) <- orig_vals_global
  trans_vals_global <- coef(po1, params_global, transform = TRUE)
  trans_vals <- trans_vals_global
  
  coef(po2, params_global) <- orig_vals_global
  
  for (i in c('hk', 'can')) {
    
    orig_vals_shared <- orig_vals[grep(paste0('^', i, '_'), params_all, value = TRUE)]
    names(orig_vals_shared) <- gsub(paste0(i, '_'), '', names(orig_vals_shared))
    
    if (i == 'hk') {
      
      coef(po1, names(orig_vals_shared)) <- orig_vals_shared
      trans_vals_shared <- coef(po1, names(orig_vals_shared), transform = TRUE)
      
    } else if (i == 'can') {
      
      coef(po2, names(orig_vals_shared)) <- orig_vals_shared
      trans_vals_shared <- coef(po2, names(orig_vals_shared), transform = TRUE)
      
    }
    
    names(trans_vals_shared) <- paste0(i, '_', names(orig_vals_shared))
    trans_vals <- c(trans_vals, trans_vals_shared)
    
  }
  
  po_save1 <- po1
  po_save2 <- po2
  
  for (i in 1:length(c(seas1, seas2))) {
    
    po1 <- po_save1
    po2 <- po_save2
    
    if (i <= length(seas1)) {
      
      orig_vals_unit <- orig_vals[grep(paste0('^', seas1[i], '_hk'), params_all, value = TRUE)]
      names(orig_vals_unit) <- gsub(paste0(seas1[i], '_hk_'), '', names(orig_vals_unit))
      
      coef(po1, names(orig_vals_unit)) <- orig_vals_unit
      trans_vals_unit <- coef(po1, names(orig_vals_unit), transform = TRUE)
      names(trans_vals_unit) <- paste0(seas1[i], '_hk_', names(orig_vals_unit))
      
    } else {
      
      orig_vals_unit <- orig_vals[grep(paste0('^', seas2[i - length(seas1)], '_can'), params_all, value = TRUE)]
      names(orig_vals_unit) <- gsub(paste0(seas2[i - length(seas1)], '_can_'), '', names(orig_vals_unit))
      
      coef(po2, names(orig_vals_unit)) <- orig_vals_unit
      trans_vals_unit <- coef(po2, names(orig_vals_unit), transform = TRUE)
      names(trans_vals_unit) <- paste0(seas2[i - length(seas1)], '_can_', names(orig_vals_unit))
      
    }
    
    trans_vals <- c(trans_vals, trans_vals_unit)
    
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po1, po2, seas1, seas2, params_all, params_global, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po1: pomp object with correct parameter transformations (Hong Kong)
  # param po2: pomp object with correct parameter transformations (Canada)
  # param seas1: Vector containing all seasons of interest (Hong Kong)
  # params seas2: Vector containing all seasons of interest (Canada)
  # param params_all: Names of all parameters to be estimated
  # param params_global: Names of all global parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_global <- trans_vals[which(params_global %in% params_all)]
  
  coef(po1, params_global, transform = TRUE) <- trans_vals_global
  orig_vals_global <- coef(po1, params_global)
  orig_vals <- orig_vals_global
  
  coef(po2, params_global, transform = TRUE) <- trans_vals_global
  
  for (i in c('hk', 'can')) {
    
    trans_vals_shared <- trans_vals[grep(paste0('^', i, '_'), params_all, value = TRUE)]
    names(trans_vals_shared) <- gsub(paste0(i, '_'), '', names(trans_vals_shared))
    
    if (i == 'hk') {
      
      coef(po1, names(trans_vals_shared), transform = TRUE) <- trans_vals_shared
      orig_vals_shared <- coef(po1, params_shared)
      
    } else if (i == 'can') {
      
      coef(po2, names(trans_vals_shared), transform = TRUE) <- trans_vals_shared
      orig_vals_shared <- coef(po2, params_shared)
      
    }
    
    names(orig_vals_shared) <- paste0(i, '_', names(trans_vals_shared))
    orig_vals <- c(orig_vals, orig_vals_shared)
    
  }
  
  po_save1 <- po1
  po_save2 <- po2
  
  for (i in 1:length(c(seas1, seas2))) {
    
    po1 <- po_save1
    po2 <- po_save2
    
    if (i <= length(seas1)) {
      
      trans_vals_unit <- trans_vals[grep(paste0('^', seas1[i], '_hk'), params_all, value = TRUE)]
      names(trans_vals_unit) <- gsub(paste0(seas1[i], '_hk_'), '', names(trans_vals_unit))
      
      coef(po1, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
      orig_vals_unit <- coef(po1, names(trans_vals_unit))
      
      names(orig_vals_unit) <- paste0(seas1[i], '_hk_', names(trans_vals_unit))
      
    } else {
      
      trans_vals_unit <- trans_vals[grep(paste0('^', seas2[i - length(seas1)], '_can'), params_all, value = TRUE)]
      names(trans_vals_unit) <- gsub(paste0(seas2[i - length(seas1)], '_can_'), '', names(trans_vals_unit))
      
      coef(po2, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
      orig_vals_unit <- coef(po2, names(trans_vals_unit))
      
      names(orig_vals_unit) <- paste0(seas2[i - length(seas1)], '_can_', names(trans_vals_unit))
      
    }
    
    orig_vals <- c(orig_vals, orig_vals_unit)
    
  }
  
  return(orig_vals)
}

# ---------------------------------------------------------------------------------------------------------------------

# Load pomp objects and prepare likelihood fxn

# Loop through years and construct pomp models:
po_list <- vector('list', length(c(seasons_hk, seasons_can)))
for (yr_index in 1:length(po_list)) {
  
  if (yr_index <= length(seasons_hk)) {
    
    yr <- seasons_hk[yr_index]
    fit_canada <- FALSE
    vir1 <- 'flu_h1_plus_b'
    
  } else {
    
    yr <- seasons_can[yr_index - length(seasons_hk)]
    fit_canada <- TRUE
    vir1 <- 'flu'
    
  }
  
  print(fit_canada)
  print(yr)
  
  # Load data and create pomp object:
  source('src/resp_interaction_model.R')
  
  # Check whether any data for given season:
  if (exists('resp_mod')) {
    
    # Add pomp object to list:
    po_list[[yr_index]] <- resp_mod
    
  }
  
  # Remove pomp object before repeating loop:
  rm(resp_mod)
  
}

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

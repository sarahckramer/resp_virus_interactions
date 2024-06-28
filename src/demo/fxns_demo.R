# ---------------------------------------------------------------------------------------------------------------------
# Functions necessary for fitting to several seasons simultaneously
# ---------------------------------------------------------------------------------------------------------------------

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
  
  # Split params into shared and unit params:
  unit_in <- as.data.frame(sapply(paste0('^', seasons, '_'), grep, x = names(trans_vals)))
  if (ncol(unit_in) == 1) {
    unit_in <- as.data.frame(matrix(data = c(unlist(unit_in)), nrow = 1))
  }
  names(unit_in) <- seasons
  shared_params <- trans_vals[-unique(unlist(unit_in))]
  
  unit_params <- list()
  for (i in 1:length(seasons)) {
    unit <- trans_vals[unit_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 2]
      unit_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season:
  units_ll <- rep(NA, length(seasons))
  
  for (i in 1:length(seasons)) {
    if (length(unit_params) > 0) {
      params_temp <- c(shared_params, unit_params[[i]])
    } else {
      params_temp <- c(shared_params)
    }
    
    units_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Calculate global -ll:
  glob_ll <- sum(units_ll)
  return(glob_ll)
}

transform_params <- function(orig_vals, po, seas, params_all, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_shared <- orig_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared) <- orig_vals_shared
  trans_vals_shared <- coef(po, params_shared, transform = TRUE)
  trans_vals <- trans_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    orig_vals_unit <- orig_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(orig_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(orig_vals_unit))
    
    coef(po, names(orig_vals_unit)) <- orig_vals_unit
    trans_vals_unit <- coef(po, names(orig_vals_unit), transform = TRUE)
    names(trans_vals_unit) <- paste0(seas[i], '_', names(orig_vals_unit))
    trans_vals <- c(trans_vals, trans_vals_unit)
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po, seas, params_all, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_shared <- trans_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared, transform = TRUE) <- trans_vals_shared
  orig_vals_shared <- coef(po, params_shared)
  orig_vals <- orig_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    trans_vals_unit <- trans_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(trans_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(trans_vals_unit))
    
    coef(po, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
    orig_vals_unit <- coef(po, names(trans_vals_unit))
    
    names(orig_vals_unit) <- paste0(seas[i], '_', names(trans_vals_unit))
    orig_vals <- c(orig_vals, orig_vals_unit)
  }
  
  return(orig_vals)
}

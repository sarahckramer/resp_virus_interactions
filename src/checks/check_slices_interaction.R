# ---------------------------------------------------------------------------------------------------------------------
# Slice over theta_lambda1/theta_lambda2/delta at MLEs
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# ---------------------------------------------------------------------------------------------------------------------

# Global estimates

# Set virus 1:
vir1 <- 'flu_B'

# Get list of results files:
res_files <- list.files(path = 'results/300721_CHECK_fluB_noGlobalParams/', full.names = TRUE)

# Read in results:
res_full = list()
for (i in seq_along(res_files)) {
  res_full[[i]] <- read_rds(res_files[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
pars_df <- lapply(res_full, getElement, 'estpars') %>%
  bind_rows() %>%
  bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
              unlist())
expect_true(nrow(pars_df) == length(res_files))
expect_true(all(is.finite(pars_df$loglik)))

# Keep only MLE:
pars_df <- pars_df %>%
  arrange(desc(loglik))

estpars <- names(pars_df)[1:35]
mle <- setNames(object = c(1.0, 1.0, 7/5, as.numeric(pars_df[1, 1:length(estpars)])),
                nm = c('theta_lambda1', 'theta_lambda2', 'delta', estpars))

# Clean up:
rm(res_files, res_full)

# Slice likelihood over theta_lambda1:
slices <- slice_design(center = mle,
                       theta_lambda1 = c(seq(from = 0, to = 1.0, by = 0.05),
                                         seq(from = 1.0, to = 2.0, by = 0.1)),
                       theta_lambda2 = c(seq(from = 0, to = 1.0, by = 0.05),
                                         seq(from = 1.0, to = 2.0, by = 0.1))) %>%#,
                       # delta = 7 / seq(from = 30, to = 1, by = -1)) %>% # delta when theta's=1 is not informative!
  mutate(ll = NA)

# Calculate likelihoods:
estpars <- names(mle)
shared_estpars <- c('theta_lambda1', 'theta_lambda2', 'delta')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
true_estpars <- c(shared_estpars, unit_estpars)

source('src/functions/setup_global_likelilhood.R')

for (i in 1:nrow(slices)) {
  x0 <- slices[i, 1:(length(estpars))]
  x0_trans <- transform_params(x0, resp_mod, seasons, estpars, shared_estpars)
  slices$ll[i] <- -1 * calculate_global_loglik(x0_trans)
}
rm(i, x0, x0_trans)

par(mfrow = c(1, 2), bty = 'l')
for (par in shared_estpars[1:2]) {
  slices_cur <- filter(slices, slice == par)
  plot(slices_cur[[par]], slices_cur$ll, type = 'l',
       xlab = par, ylab = 'Log-Likelihood',
       main = par)
}
rm(par, slices_cur)
# Slices suggest a positive impact of flu on RSV and about a null impact of RSV on flu, but note that this is at the
# MLE when no interaction is assumed.

# Either way, this seems to confirm (or at least support) the idea that theta_lambda2 will be much more identifiable,
# since it is stabilizing around 1.0 whereas theta_lambda1 is not, and in fact has similar ll values across a wide
# range of parameter values.

# ---------------------------------------------------------------------------------------------------------------------

# By individual season

# Specify parameters estimated:
estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

# Read in results:
mle_list <- read_rds('results/traj_match_round1_byvirseas_MLE.rds')
pars_top <- read_rds('results/traj_match_round1_byvirseas_TOP.rds')

# Limit to vir1 of interest:
mle_list <- mle_list[which(str_detect(names(mle_list), vir1))]

# Loop through seasons:
seasons <- 2006:2014
counter <- 1
for (yr in seasons) {
  
  vir_seas <- paste(vir1, yr, sep = '_')
  if (vir_seas %in% names(mle_list)) {
    
    # Get MLE:
    mle <- c(1.0, 1.0, 7 / 5, mle_list[[vir_seas]])
    names(mle)[1:3] <- c('theta_lambda1', 'theta_lambda2', 'delta')
    
    # Set up slices:
    slices <- slice_design(center = mle,
                           theta_lambda1 = c(seq(from = 0, to = 1.0, by = 0.05),
                                             seq(from = 1.0, to = 3.0, by = 0.1)),
                           theta_lambda2 = c(seq(from = 0, to = 1.0, by = 0.05),
                                             seq(from = 1.0, to = 3.0, by = 0.1))) %>%
      mutate(ll = NA)
    
    # Calculate slice likelihoods:
    obj_fun <- obj_fun_list[[counter]]
    
    for (i in 1:nrow(slices)) {
      x0 <- slices[i, 1:(length(estpars) + 3)]
      coef(resp_mod, names(mle)) <- slices[i, 1:length(names(mle))]
      x0_trans <- coef(resp_mod, names(mle), transform = TRUE)
      slices$ll[i] <- -1 * obj_fun(x0_trans)
    }
    
    # Plot slices:
    par(mfrow = c(1, 2), bty = 'l')
    for (par in shared_estpars[1:2]) {
      slices_cur <- filter(slices, slice == par)
      plot(slices_cur[[par]], slices_cur$ll, type = 'l',
           xlab = par, ylab = 'Log-Likelihood',
           main = yr)
    }
    rm(par, slices_cur)
    
    # Iterate counter:
    counter <- counter + 1
    
  }
}

# For theta_lambda2, values tend to fit around 1.0 or above; range of likelihoods very small for 2006/2009/2011, but
# wider for 2008/2013; these seem to be the seasons with a medium level of overlap between RSV and flu_B?

# For theta_lambda1, not well identified, and for most seasons higher yields higher likelihood

# Seems like some seasons are more informative for theta_lambda2 than others, and it definitely looks more identifiable
# when we consider all seasons at once

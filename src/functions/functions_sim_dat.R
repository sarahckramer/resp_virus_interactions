# ---------------------------------------------------------------------------------------------------------------------
# Functions to generate/assess synthetic data
# ---------------------------------------------------------------------------------------------------------------------

generate_param_lhs <- function(nm_pars, bound_pars, lhs_iter, po, adjust_init) {
  
  # # Set seed:
  # set.seed(seed)
  
  # Draw from ranges:
  parms <- lhs(lhs_iter, bound_pars) %>%
    as.data.frame() %>%
    as_tibble()
  names(parms) <- nm_pars
  
  # Correct where R10+R20+R120 > 0.99:
  # if (all(c('R10', 'R20', 'R120') %in% nm_pars)) {
  if (adjust_init) {
    
    parms <- parms %>%
      mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
      mutate(R10 = if_else(sum >= 1.0, R10 - ((sum - 0.9995) * (R10 / sum)), R10),
             R20 = if_else(sum >= 1.0, R20 - ((sum - 0.9995) * (R20 / sum)), R20),
             R120 = if_else(sum >= 1.0, R120 - ((sum - 0.9995) * (R120 / sum)), R120)) %>%
      # mutate(sum_new = I10 + I20 + R10 + R20 + R120) %>%
      select(-sum)
    
    check_sum <- parms %>%
      mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
      pull(sum) %>%
      max()
    
    while (check_sum > 1.0) {
      parms <- parms %>%
        mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
        mutate(R10 = if_else(sum >= 1.0, R10 - ((sum - 0.9995) * (R10 / sum)), R10),
               R20 = if_else(sum >= 1.0, R20 - ((sum - 0.9995) * (R20 / sum)), R20),
               R120 = if_else(sum >= 1.0, R120 - ((sum - 0.9995) * (R120 / sum)), R120)) %>%
        # mutate(sum_new = I10 + I20 + R10 + R20 + R120) %>%
        select(-sum)
      
      check_sum <- parms %>%
        mutate(sum = I10 + I20 + R10 + R20 + R120) %>%
        pull(sum) %>%
        max()
    }
    
  }
  
  # Set values in parameter input matrix:
  parms_mat <- parmat(params = coef(po), nrep = lhs_iter)
  for (ix in 1:length(nm_pars)) {
    parms_mat[nm_pars[ix], ] <- as.vector(pull(parms[, ix]))
  }
  
  # Set identifiers in parms:
  parms$id.parm <- 1:(dim(parms)[1])
  
  # Return parameter sets:
  return(list(parms, parms_mat))
  
}

run_and_format_simulations <- function(po, parms_mat, parms, sim_num) {
  
  # Run simulations:
  sim_dat <- simulate(object = po,
                      params = parms_mat,
                      nsim = sim_num,
                      format = 'data.frame')
  
  # Format results:
  sim_dat <- sim_dat %>%
    as_tibble() %>%
    select(time:.id, H1_tot:n_P2) %>%
    mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
           id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2)))
  
  # Calculate attack rates, peak timings, and differences:
  sim_metrics <- sim_dat %>%
    group_by(.id) %>%
    summarise(ar1 = sum(n_P1), ar2 = sum(n_P2),
              pt1 = which.max(n_P1) + 39, pt2 = which.max(n_P2) + 39,
              pt_diff = pt1 - pt2)
  # pt_diff describes how many weeks flu peak occurs AFTER the RSV peak (negative values mean flu came first)
  
  # Join with identifying info:
  sim_metrics <- sim_metrics %>%
    left_join(sim_dat %>%
                select(.id, id.parm, id.sim) %>%
                unique(),
              by = '.id')
  
  # And with parameter value info:
  sim_dat <- sim_dat %>% left_join(parms, by = 'id.parm')
  sim_metrics <- sim_metrics %>% left_join(parms, by = 'id.parm')
  
  # Return results:
  return(list(sim_dat, sim_metrics))
  
}

classify_realistic <- function(sim_dat, sim_metrics, lhs_iter) {
  
  # Create column to store info on whether param set produces realistic outbreaks or not:
  sim_dat$real <- NA
  sim_metrics$real <- NA
  
  # Identify if either or both viruses have no discernible outbreak in 50% or more of simulations:
  parm_sets_no_outbreak <- sim_metrics %>%
    mutate(no_outbreak1 = ar1 < 100,
           no_outbreak2 = ar2 < 100) %>%
    group_by(id.parm, no_outbreak1, no_outbreak2) %>%
    summarise(count_low_ar = n()) %>%
    filter(no_outbreak1 | no_outbreak2) %>%
    group_by(id.parm) %>%
    summarise(total_low = sum(count_low_ar)) %>%
    filter(total_low >= n_sim / 2) %>%
    pull(id.parm)
  
  # Identify parameter sets that can produce realistic outbreaks:
  print(lhs_iter - length(parm_sets_no_outbreak))
  
  ids_realistic <- sim_metrics %>%
    filter(pt_diff >= 0 & pt_diff < 14) %>%
    filter(ar2 < 500 & ar1 < 1500 & ar1 > 200) %>%
    filter(pt1 >= 53 & pt1 <= 63 &
             pt2 >= 46 & pt2 <= 55) %>%
    pull(id.parm) %>%
    unique()
  ids_realistic <- ids_realistic[!(ids_realistic %in% parm_sets_no_outbreak)]
  print(length(ids_realistic) / lhs_iter)
  
  sim_dat <- sim_dat %>% mutate(real = if_else(id.parm %in% ids_realistic, TRUE, FALSE))
  sim_metrics <- sim_metrics %>% mutate(real = if_else(id.parm %in% ids_realistic, TRUE, FALSE))
  
  # Return results:
  return(list(sim_dat, sim_metrics))
  
}

plot_res <- function(sim_met, plot_title = '', incl_R2 = TRUE) {
  
  # Plot parameter ranges/relationships by whether outbreaks are realistic vs. not:
  if (incl_R2) {
    sim_plot <- sim_met %>%
      filter(id.sim == 1) %>%
      mutate(R10.tot = R10 + R120,
             R20.tot = R20 + R120)
  } else {
    sim_plot <- sim_met %>%
      filter(id.sim == 1)
  }
  
  colors <- rep('gray80', n_lhs); sizes <- rep(1.0, n_lhs)
  colors[sim_plot$real] <- 'steelblue2'; sizes[sim_plot$real] <- 1.5
  
  cols_to_plot <- names(sim_plot %>%
                          select(-c(.id:id.sim, real)))
  pairs(sim_plot %>% select(all_of(cols_to_plot)), pch = 20, col = colors, cex = sizes, main = plot_title)
  
  # Limit to just realistic outbreaks and plot parameter ranges:
  sim_plot <- sim_plot %>%
    filter(real)
  # pairs(sim_plot %>% select(all_of(cols_to_plot)), pch = 20, col = 'steelblue2', cex = 1.25, main = plot_title)
  
  sim_plot <- sim_plot %>%
    select(-real) %>%
    pivot_longer(all_of(cols_to_plot), names_to = 'param', values_to = 'val')
  p1 <- ggplot(data = sim_plot) + geom_violin(aes(x = '', y = val), fill = 'steelblue2', alpha = 0.5) +
    facet_wrap(~ param, scales = 'free_y') + theme_classic() +
    labs(x = '', y = 'Param. Value', title = plot_title)
  print(p1)
  
}

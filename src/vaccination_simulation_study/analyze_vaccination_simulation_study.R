# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(patchwork)
library(viridis)

# List seasons:
seasons <- c('s13-14', 's15-16', 's16-17', 's17-18', 's18-19')

# Set vaccination coverage levels and time points:
vacc_cov_vec <- round(seq(0.05, 1.0, by = 0.05), digits = 2) # seq(0.1, 1.0, by = 0.1)
vacc_time_vec <- round(seq(0, 52, by = 1)) # seq(0, 52, by = 2)
 
# # Set vaccination efficacy against flu:
# vacc_eff <- 0.8
# 
# # Set parameters for run:
# vir1 <- 'flu_h1'
# vir2 <- 'rsv'
# 
# Ri_max1 <- 2.0
# Ri_max2 <- 3.0
# d2_max <- 10.0
# 
# debug_bool <- FALSE

# ---------------------------------------------------------------------------------------------------------------------

# Compile results and assess impact (Hong Kong / "subtropical" scenario)

# Read in all simulation results:
file_list <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'SUBTROPICAL', full.names = TRUE)

res_list <- vector('list', length(file_list))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_split(file_list[i], '_')[[1]][5])
  # res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_sub(file_list[i], 62, 67))
}
rm(i)

# Combine into single tibble:
res <- bind_rows(res_list) %>%
  as_tibble()

# Calculate measures of vaccine impact from simulation results:
res_metrics <- res %>%
  group_by(season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

# Want vacc_time relative to timing of flu/RSV peak?
# Get peak timing of flu and RSV without vaccination:
pt_novacc <- res %>%
  filter(.id == 1) %>%
  group_by(season, vacc_cov, vacc_time) %>%
  summarise(pt1 = which.max(H1), pt2 = which.max(H2)) %>%
  ungroup() %>%
  select(!vacc_cov) %>%
  unique()
expect_true(nrow(pt_novacc) == length(seasons) * length(vacc_time_vec))
expect_true(nrow(pt_novacc %>% select(!vacc_time) %>% unique()) == length(seasons))
rm(pt_novacc)

# Write results to file:
write_rds(res_metrics, 'results/vaccine_simulation_study/res_METRICS_simulation_study_SUBTROPICAL.rds')

# Limit results to coverage levels of 60% and lower:
res_metrics <- res_metrics %>%
  filter(vacc_cov <= 0.60)

# Plot results:
res_metrics <- res_metrics %>%
  pivot_longer(ar1_impact:ar2_impact, names_to = 'metric', values_to = 'val')

res_metrics_AVG <- res_metrics %>%
  group_by(vacc_cov, vacc_time, metric) %>%
  summarise(val = median(val))

upper_bound_ar <- ceiling(max(res_metrics_AVG$val[res_metrics_AVG$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10

facet_labs <- c('Influenza', 'RSV')
names(facet_labs) <- c('ar1_impact', 'ar2_impact')

p_lines_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
  geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', col = 'Vacc. Cov.')

p_heat_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = vacc_cov, fill = val)) +
  geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + #theme(axis.ticks = element_blank()) +
  scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
       fill = 'Impact')

print(p_lines_ar / p_heat_ar)
p_hk <- p_lines_ar / p_heat_ar
ggsave('results/vaccine_simulation_study/plots/ar_impact_AVG.svg', p_hk, width = 9.5, height = 6)

# Plot results by individual season:
pdf('results/vaccine_simulation_study/plots/vacc_impact_AR.pdf',
    width = 20, height = 10)
for (seas in seasons) {
  
  sim_temp <- res %>%
    filter(season == seas, .id == 1, vacc_cov == 1, vacc_time == 0) %>%
    select(time:H2) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val,
           Virus = if_else(Virus == 'H1', 'Flu', 'RSV'))
  res_temp <- res_metrics %>% filter(season == seas)
  
  p_outbreak <- ggplot(data = sim_temp, aes(x = time, y = val, color = Virus)) +
    geom_line() + geom_point() +
    theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  upper_bound_ar <- ceiling(max(res_temp$val[res_temp$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10
  
  p_lines_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
    geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', title = seas, col = 'Vacc. Cov.')
  
  p_heat_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = vacc_cov, fill = val)) +
    geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + #theme(axis.ticks = element_blank()) +
    scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                         breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
         title = seas, fill = 'Impact')
  
  if (upper_bound_ar == 1.0) {
    
    p_heat_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = vacc_cov, fill = val)) +
      geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
      theme_classic() +
      scale_fill_distiller(palette = 'RdBu', values = c(0, 1, 1.1), limits = c(0, upper_bound_ar),
                           breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
           title = seas, fill = 'Impact')
    
  }
  
  print(p_outbreak | (p_lines_ar / p_heat_ar))
  
}
dev.off()

# Plot all outbreaks by vaccination timing and coverage:
pdf('results/vaccine_simulation_study/plots/all_sims.pdf',
    width = 30, height = 11)
for (seas in seasons) {
  
  sim_temp_all <- res %>%
    filter(vacc_cov <= 0.6) %>%
    filter(season == seas & .id == 2) %>%
    select(-c(.id, season)) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val)
  
  p_outbreak_all <- ggplot(data = sim_temp_all, aes(x = time, y = val, col = Virus)) + geom_line() +
    facet_grid(vacc_cov ~ vacc_time) + theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  print(p_outbreak_all)
  
}
dev.off()

# For each scenario, calculate correlations to determine consistency between seasons:
vec_hk <- vector('list', length = length(unique(res$season)))
for (i in 1:length(unique(res_metrics$season))) {
  
  seas <- unique(res_metrics$season)[i]
  
  vec_hk[[i]] <- res_metrics %>%
    filter(season == seas,
           metric == 'ar2_impact') %>%
    arrange(vacc_cov, vacc_time) %>%
    pull(val)
  
}

corr_mat_hk <- matrix(NA, nrow = 5, ncol = 5)
for (i in 1:5) {
  for (j in 1:5) {
    corr_mat_hk[i, j] <- cor.test(vec_hk[[i]], vec_hk[[j]], method = 'kendall')$estimate
  }
}
rownames(corr_mat_hk) = colnames(corr_mat_hk) = unique(res_metrics$season)

# ---------------------------------------------------------------------------------------------------------------------

# Compile results and assess impact ("temperate" scenario)

# Read in all simulation results:
file_list <- list.files('results/vaccine_simulation_study/simulations/main/', pattern = 'TEMPERATE', full.names = TRUE)

res_list <- vector('list', length(file_list))
for (i in 1:length(res_list)) {
  res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_split(file_list[i], '_')[[1]][5])
  # res_list[[i]] <- read_rds(file_list[i]) %>% mutate(season = str_sub(file_list[i], 62, 67))
}
rm(i)

# Combine into single tibble:
res <- bind_rows(res_list) %>%
  as_tibble()

# Calculate measures of vaccine impact from simulation results:
res_metrics <- res %>%
  group_by(season, vacc_cov, vacc_time, .id) %>%
  summarise(ar1 = sum(H1), ar2 = sum(H2)) %>%
  ungroup() %>%
  group_by(season, vacc_cov, vacc_time) %>%
  summarise(ar1_impact = ar1[.id == 2] / ar1[.id == 1],
            ar2_impact = ar2[.id == 2] / ar2[.id == 1]) %>%
  ungroup()

# Write results to file:
write_rds(res_metrics, 'results/vaccine_simulation_study/res_METRICS_simulation_study_TEMPERATE.rds')

# Limit results to coverage levels of 60% and lower:
res_metrics <- res_metrics %>%
  filter(vacc_cov <= 0.60)

# Plot results:
res_metrics <- res_metrics %>%
  pivot_longer(ar1_impact:ar2_impact, names_to = 'metric', values_to = 'val')

res_metrics_AVG <- res_metrics %>%
  group_by(vacc_cov, vacc_time, metric) %>%
  summarise(val = median(val))

upper_bound_ar <- ceiling(max(res_metrics_AVG$val[res_metrics_AVG$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10

facet_labs <- c('Influenza', 'RSV')
names(facet_labs) <- c('ar1_impact', 'ar2_impact')

p_lines_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
  geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', col = 'Vacc. Cov.')

p_heat_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = vacc_cov, fill = val)) +
  geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
  theme_classic() + #theme(axis.ticks = element_blank()) +
  scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                       breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
       fill = 'Impact')

if (upper_bound_ar == 1.0) {
  
  p_heat_ar <- ggplot(data = res_metrics_AVG, aes(x = vacc_time, y = vacc_cov, fill = val)) +
    geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + #theme(axis.ticks = element_blank()) +
    scale_fill_distiller(palette = 'RdBu', values = c(0, 1, 1.1), limits = c(0, upper_bound_ar),
                         breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
         fill = 'Impact')
  
}

print(p_lines_ar / p_heat_ar)
p_temp <- p_lines_ar / p_heat_ar
ggsave('results/vaccine_simulation_study/plots/ar_impact_AVG_temperate.svg', p_temp, width = 9.5, height = 6)

# Plot results by individual season:
pdf('results/vaccine_simulation_study/plots/vacc_impact_AR_temperate.pdf',
    width = 20, height = 10)
for (seas in seasons) {
  
  sim_temp <- res %>%
    filter(season == seas, .id == 1, vacc_cov == 1, vacc_time == 0) %>%
    select(time:H2) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val,
           Virus = if_else(Virus == 'H1', 'Flu', 'RSV'))
  res_temp <- res_metrics %>% filter(season == seas)
  
  p_outbreak <- ggplot(data = sim_temp, aes(x = time, y = val, color = Virus)) +
    geom_line() + geom_point() +
    theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  upper_bound_ar <- ceiling(max(res_temp$val[res_temp$metric %in% c('ar1_impact', 'ar2_impact')]) * 10) / 10
  
  # facet_labs <- c('Influenza', 'RSV')
  # names(facet_labs) <- c('ar1_impact', 'ar2_impact')
  
  p_lines_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = val, col = vacc_cov, group = vacc_cov)) +
    geom_line() + facet_wrap(~ metric, scales = 'free_y', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + scale_color_viridis() + scale_y_continuous(limits = c(0, upper_bound_ar)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'AR_Vacc / AR_NoVacc', title = seas, col = 'Vacc. Cov.')
  
  p_heat_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = vacc_cov, fill = val)) +
    geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
    theme_classic() + #theme(axis.ticks = element_blank()) +
    scale_fill_distiller(palette = 'RdBu', values = c(0, 1 / upper_bound_ar, 1), limits = c(0, upper_bound_ar),
                         breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
         title = seas, fill = 'Impact')
  
  if (upper_bound_ar == 1.0) {
    
    p_heat_ar <- ggplot(data = res_temp, aes(x = vacc_time, y = vacc_cov, fill = val)) +
      geom_tile() + facet_wrap(~ metric, scales = 'free', ncol = 2, labeller = labeller(metric = facet_labs)) +
      theme_classic() + #theme(axis.ticks = element_blank()) +
      scale_fill_distiller(palette = 'RdBu', values = c(0, 1, 1.1), limits = c(0, upper_bound_ar),
                           breaks = c(0, 0.25, 0.5, 0.75, seq(1.0, upper_bound_ar, by = 0.25))) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      labs(x = 'Week of Instantaneous Vaccination', y = 'Vaccine Coverage (%)',
           title = seas, fill = 'Impact')
    
  }
  
  print(p_outbreak | (p_lines_ar / p_heat_ar))
  
}
dev.off()

# Plot all outbreaks by vaccination timing and coverage:
pdf('results/vaccine_simulation_study/plots/all_sims_temperate.pdf',
    width = 30, height = 11)
for (seas in seasons) {
  
  sim_temp_all <- res %>%
    filter(vacc_cov <= 0.6) %>%
    filter(season == seas & .id == 2) %>%
    select(-c(.id, season)) %>%
    pivot_longer(H1:H2, names_to = 'Virus', values_to = 'val') %>%
    mutate(val = 100 * val)
  
  p_outbreak_all <- ggplot(data = sim_temp_all, aes(x = time, y = val, col = Virus)) + geom_line() +
    facet_grid(vacc_cov ~ vacc_time) + theme_classic() + scale_color_brewer(palette = 'Set1') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)', title = seas)
  
  print(p_outbreak_all)
  
}
dev.off()

# For each scenario, calculate correlations to determine consistency between seasons:
vec_temp <- vector('list', length = length(unique(res$season)))
for (i in 1:length(unique(res_metrics$season))) {
  
  seas <- unique(res_metrics$season)[i]
  
  vec_temp[[i]] <- res_metrics %>%
    filter(season == seas,
           metric == 'ar2_impact') %>%
    arrange(vacc_cov, vacc_time) %>%
    pull(val)
  
}

corr_mat_temp <- matrix(NA, nrow = 5, ncol = 5)
for (i in 1:5) {
  for (j in 1:5) {
    corr_mat_temp[i, j] <- cor.test(vec_temp[[i]], vec_temp[[j]], method = 'kendall')$estimate
  }
}
rownames(corr_mat_temp) = colnames(corr_mat_temp) = unique(res_metrics$season)

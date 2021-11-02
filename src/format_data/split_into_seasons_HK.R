# ---------------------------------------------------------------------------------------------------------------------
# Divide Hong Kong data into individual "seasons"
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(testthat)

# Read in data:
dat_hk <- read_csv('data/formatted/dat_hk.csv')

# Limit (for now) to data of interest:
dat_hk_ORIG <- dat_hk
dat_hk <- dat_hk_ORIG %>%
  select(Time:n_rsv, n_rhino_est_rnd)

# Calculate proportion rather than number positive:
dat_hk <- dat_hk %>%
  mutate(prop_h1 = n_h1 / n_samp,
         prop_h3 = n_h3 / n_samp,
         prop_b = n_b / n_samp,
         prop_rsv = n_rsv / n_samp,
         prop_rhino = n_rhino_est_rnd / n_samp) %>%
  select(Time:Week, prop_h1:prop_rhino)

# Check for NAs:
summary(dat_hk)

# List viruses:
vir_list <- c('prop_h1', 'prop_h3', 'prop_b', 'prop_rsv', 'prop_rhino')

# Initiate lists to store results:
onset_list <- vector('list', length(vir_list))
end_list <- vector('list', length(vir_list))

# Explore onset thresholds for flu viruses:
baseline_perc <- 0.5

par(mfrow = c(3, 1))
for (vir_index in 1:3) {
  
  # Get virus:
  vir <- vir_list[vir_index]
  
  # Subset data:
  dat_temp <- dat_hk %>%
    select(Time:Week, all_of(vir)) %>%
    rename(prop_vir = vir)
  obs_i <- dat_temp$prop_vir
  
  # Get threshold value:
  baseline <- quantile(obs_i[obs_i != 0], probs = baseline_perc, na.rm = TRUE)
  
  # Find onsets and ends:
  onsets = ends = c()
  above.baseline.count <- 0
  below.baseline.count <- 0
  
  for (i in 1:length(obs_i)) {
    if (above.baseline.count < 3) { # start searching for onset
      if (obs_i[i] > baseline) {
        above.baseline.count <- above.baseline.count + 1
      } else {
        above.baseline.count <- 0
      }
    }
    
    if (above.baseline.count == 3) { # record onset
      onsets <- c(onsets, i - 2)
      above.baseline.count <- 10 # indicator that an onset has occurred
    }
    
    if (above.baseline.count >= 3) { # start searching for end
      if (obs_i[i] < baseline) {
        below.baseline.count <- below.baseline.count + 1
      } else {
        below.baseline.count <- 0
      }
      
      if (below.baseline.count == 2) { # record end
        ends <- c(ends, i - 1)
        above.baseline.count <- 0
        below.baseline.count <- 0
      }
    }
    
  }
  
  if (length(onsets) != length(ends)) {
    ends <- c(ends, length(obs_i))
  }
  
  # Remove where outbreaks are too small (<3x baseline?):
  if (length(onsets) > 0) {
    for (j in 1:length(onsets)) {
      obs_temp <- obs_i[onsets[j]:ends[j]]
      if (!any(obs_temp > baseline * 3)) {
        onsets[j] <- NA
        ends[j] <- NA
      }
    }
    
    onsets <- onsets[!is.na(onsets)]
    ends <- ends[!is.na(ends)]
  }
  
  # # Add buffer of 4 weeks before onset and after end:
  # onsets <- onsets - 4
  # onsets[onsets < 1] <- 1
  # 
  # ends <- ends + 4
  # ends[ends > length(obs_i)] <- length(obs_i)
  
  # Plot data and onsets/ends:
  plot(obs_i, pch = 20, type = 'b', xlab = '', ylab = 'Prop.')
  abline(h = baseline, lty = 2)
  abline(v = onsets, col = 'blue')
  abline(v = ends, col = 'green')
  
  # Store results:
  onset_list[[vir_index]] <- onsets
  end_list[[vir_index]] <- ends
  
}

# Explore onset thresholds for RSV/rhinovirus:
baseline_perc <- 0.3

par(mfrow = c(2, 1))
for (vir_index in 4:5) {
  
  # Get virus:
  vir <- vir_list[vir_index]
  
  # Subset data:
  dat_temp <- dat_hk %>%
    select(Time:Week, all_of(vir)) %>%
    rename(prop_vir = vir)
  obs_i <- dat_temp$prop_vir
  
  # Get threshold value:
  baseline <- quantile(obs_i[obs_i != 0], probs = baseline_perc, na.rm = TRUE)
  if (vir == 'prop_rhino') {
    baseline <- quantile(obs_i[7:length(obs_i)][obs_i[7:length(obs_i)] != 0], probs = baseline_perc, na.rm = TRUE)
  }
  
  # Find onsets and ends:
  onsets = ends = c()
  above.baseline.count <- 0
  below.baseline.count <- 0
  
  for (i in 1:length(obs_i)) {
    if (above.baseline.count < 3) { # start searching for onset
      if (obs_i[i] > baseline) {
        above.baseline.count <- above.baseline.count + 1
      } else {
        above.baseline.count <- 0
      }
    }
    
    if (above.baseline.count == 3) { # record onset
      onsets <- c(onsets, i - 2)
      above.baseline.count <- 10 # indicator that an onset has occurred
    }
    
    if (above.baseline.count >= 3) { # start searching for end
      if (obs_i[i] < baseline) {
        below.baseline.count <- below.baseline.count + 1
      } else {
        below.baseline.count <- 0
      }
      
      if (below.baseline.count == 2) { # record end
        ends <- c(ends, i - 1)
        above.baseline.count <- 0
        below.baseline.count <- 0
      }
    }
    
  }
  
  if (length(onsets) != length(ends)) {
    ends <- c(ends, length(obs_i))
  }
  
  # Remove where outbreaks are too small (<3x baseline?):
  if (length(onsets) > 0) {
    for (j in 1:length(onsets)) {
      obs_temp <- obs_i[onsets[j]:ends[j]]
      if (!any(obs_temp > baseline * 1.5)) {
        onsets[j] <- NA
        ends[j] <- NA
      }
    }
    
    onsets <- onsets[!is.na(onsets)]
    ends <- ends[!is.na(ends)]
  }
  
  # # Add buffer of 4 weeks before onset and after end:
  # onsets <- onsets - 4
  # onsets[onsets < 1] <- 1
  # 
  # ends <- ends + 4
  # ends[ends > length(obs_i)] <- length(obs_i)
  
  # Plot data and onsets/ends:
  plot(obs_i, pch = 20, type = 'b', xlab = '', ylab = 'Prop.')
  abline(h = baseline, lty = 2)
  abline(v = onsets, col = 'blue')
  abline(v = ends, col = 'green')
  
  # Store results:
  onset_list[[vir_index]] <- onsets
  end_list[[vir_index]] <- ends
  
}

# Clean up:
rm(dat_temp, above.baseline.count, below.baseline.count, baseline, baseline_perc, onsets, ends, i, j, obs_i, obs_temp, vir_index, vir)

# Now find overlap between flu "seasons" and other viruses:

# par(mfrow = c(2, 1))
# plot(dat_hk$prop_h1, type = 'b', pch = 20)
# abline(v = onset_list[[1]], col = 'blue')
# abline(v = end_list[[1]], col = 'green')
# plot(dat_hk$prop_rhino, type = 'b', pch = 20)
# abline(v = onset_list[[5]], col = 'blue')
# abline(v = end_list[[5]], col = 'green')

outbreak_df <- NULL
for (flu_index in 1:3) {
  
  onsets <- onset_list[[flu_index]]
  ends <- end_list[[flu_index]]
  
  for (resp_index in 4:5) {
    
    onsets_2 <- onset_list[[resp_index]]
    ends_2 <- end_list[[resp_index]]
    
    end_next <- 0
    
    for (i in length(onsets):1) {
      o <- onsets[i]; e <- ends[i]
      
      # want any ends that occur during or directly after the flu outbreak:
      if (any(ends_2 >= e)) {
        outbreak_start <- o
        outbreak_end <- max(e, ends_2[min(which(ends_2 >= e))])
        
        if (min(which(ends_2 >= e)) == end_next) {
          if(end_next > 1) {
            outbreak_end <- e
          } else {
            outbreak_end <- e
          }
        }
        
        end_next <- min(which(ends_2 >= e))
        
        outbreak_df <- rbind(outbreak_df,
                             c(vir_list[flu_index],
                               vir_list[resp_index],
                               outbreak_start,
                               outbreak_end))
      } #else {
      #   if (any(ends_2 > o & ends_2 <= e)) {
      #     outbreak_start <- o
      #     outbreak_end <- e
      #     
      #     outbreak_df <- rbind(outbreak_df,
      #                          c(vir_list[flu_index],
      #                            vir_list[resp_index],
      #                            outbreak_start,
      #                            outbreak_end))
      #   }
      # }
      
    }
  }
}

# Clean up:
rm(flu_index, resp_index, i, onsets, ends, onsets_2, ends_2, o, e, outbreak_start, outbreak_end, end_next)

# Format data frame:
outbreak_df <- outbreak_df %>%
  as_tibble() %>%
  rename('vir1' = 'V1',
         'vir2' = 'V2',
         'start' = 'V3',
         'end' = 'V4') %>%
  mutate(vir1 = str_sub(vir1, 6, str_length(vir1)),
         vir2 = str_sub(vir2, 6, str_length(vir2))) %>%
  mutate(across(.cols = !starts_with('vir'), as.numeric))

# Remove where outbreaks are shorter (or longer?) than given cutoff values:
outbreak_df <- outbreak_df %>%
  mutate(duration = end - start + 1) %>%
  filter(duration >= 25)

# Add buffer of 4 weeks before onsets/after ends:
outbreak_df <- outbreak_df %>%
  mutate(start = start - 4,
         end = end + 4) %>%
  mutate(start = if_else(start < 1, 1, start),
         end = ifelse(end > nrow(dat_hk_ORIG), nrow(dat_hk_ORIG), end))

# Format data for input to models:
dat_hk_pomp <- vector('list', length = 6)

counter <- 1
for (vir1_nm in unique(outbreak_df$vir1)) {
  for (vir2_nm in unique(outbreak_df$vir2)) {
    
    # Rename list element:
    names(dat_hk_pomp)[counter] <- paste(vir1_nm, vir2_nm, sep = '_')
    
    # Get start and end times for each outbreak:
    onsets <- outbreak_df %>% filter(vir1 == vir1_nm, vir2 == vir2_nm) %>% pull(start) %>% rev()
    ends <- outbreak_df %>% filter(vir1 == vir1_nm, vir2 == vir2_nm) %>% pull(end) %>% rev()
    
    # Get data of interest:
    dat_temp <- dat_hk_ORIG %>%
      select(Time:n_samp, contains(vir1_nm), contains(vir2_nm), GOPC:PMP.Clinics) %>%
      mutate(outbreak_number = NA)
    
    if (vir2_nm == 'rhino') {
      dat_temp <- dat_temp %>%
        select(!c(n_rhino_entero, n_rhino_est)) %>%
        rename('n_rhino' = 'n_rhino_est_rnd')
    }
    
    # Build data frame containing each outbreak:
    dat_out <- NULL
    for (i in 1:length(onsets)) {
      
      dat_i <- dat_temp %>%
        filter(Time >= onsets[i] & Time <= ends[i]) %>%
        mutate(outbreak_number = i,
               Time = Time - min(Time) + 1)
      dat_out <- rbind(dat_out, dat_i)
      
    }
    
    # Check that no NAs in outbreak number:
    dat_out %>% pull(outbreak_number) %>% is.na() %>% any() %>% expect_false()
    
    # Format column names:
    dat_out <- dat_out %>%
      rename('time' = 'Time',
             'n_T' = 'n_samp') %>%
      rename_with(~ 'n_P1', .cols = contains(paste0('_', vir1_nm))) %>%
      rename_with(~ 'n_P2', .cols = contains(vir2_nm))
    
    # Store data in list:
    dat_hk_pomp[[counter]] <- dat_out
    
    # Iterate counter:
    counter <- counter + 1
    
  }
}

# Clean up:
rm(counter, vir1_nm, vir2_nm, onsets, ends, dat_temp, dat_out, i, dat_i)

# Plot outbreaks for each combination of pathogens:
plot_list <- vector('list', length = length(dat_hk_pomp))
for (i in 1:length(dat_hk_pomp)) {
  dat_temp <- dat_hk_pomp[[i]]
  
  # p_temp <- ggplot(data = dat_temp %>% pivot_longer(n_P1:n_P2, names_to = 'vir', values_to = 'val')) +
  #   geom_line(aes(x = time, y = n_T), col = 'black', lwd = 1.25) +
  #   geom_line(aes(x = time, y = val, group = vir, col = vir)) +
  #   theme_classic() + scale_color_brewer(palette = 'Set1') +
  #   facet_wrap(~ outbreak_number)
  
  p_temp <- ggplot(data = dat_temp %>%
                     mutate(prop_P1 = n_P1 / n_T, prop_P2 = n_P2 / n_T) %>%
                     pivot_longer(prop_P1:prop_P2, names_to = 'vir', values_to = 'val')
  ) +
    geom_line(aes(x = time, y = val, group = vir, col = vir)) +
    theme_classic() +
    scale_color_brewer(palette = 'Set1') +
    facet_wrap(~ outbreak_number) +
    labs(x = 'Time (Weeks)', y = 'Proportion Positive', col = 'Virus', title = names(dat_hk_pomp)[i])
  plot_list[[i]] <- p_temp
  
}
print(plot_list)

# Clean up:
rm(i, dat_temp)

# Save plots to file:
pdf('results/plots/data_Hong_Kong_byOutbreak.pdf', width = 8, height = 5)
print(plot_list)
dev.off()

# Save formatted data:
write_rds(dat_hk_pomp, file = 'data/formatted/dat_hk_byOutbreak.rds')

# Clean up:
rm(list = ls())

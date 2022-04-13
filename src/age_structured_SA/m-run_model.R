#######################################################################################################
#  Run simulations of the age-structured model of flu-RSV interaction
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
source("f-CreateInteractionMod.R")
library(socialmixr)
debug_bool <- T
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))

# Get POLYMOD data and create contact matrix -------------------------------------------------------
# Take POLYMOD matrix from the UK, corrected for reciprocity 
age_limits <- c(0, 1, 5, 16, 65) # Age limits for age categories
data(polymod)
CM_all <- bake(file = "contact_matrix_UK_symm.rds", 
               expr = {
                 contact_matrix(polymod, 
                                countries = "United Kingdom",  
                                age.limits = age_limits,
                                symmetric = T)
               })

# Plot matrix
CM <- CM_all$matrix
rownames(CM) <- colnames(CM)
n_ages <- nrow(CM)
N <- CM_all$demography$population 

# age1 represents the reporter age, age2 the contact age
CM_long <- CM %>% 
  melt(varnames = c("age1", "age2"), value.name = "contacts")

# Plot matrix; x-axis: contact age, y-axis: reporter age
# CAUTION: rates are per day!
pl <- ggplot(data = CM_long, 
             mapping = aes(x = age2, y = age1, fill = contacts)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "red") + 
  theme_minimal() + 
  #theme(legend.position = "top") + 
  labs(x = "Contact age", y = "Reported age", fill = "Daily contacts")
print(pl)

# Create pomp object ------------------------------------------------------
interMod <- CreateInteractionMod(dat = data.frame(time = seq(from = -1, to = 53, by = 0.1), temp = 0, ah = 0), 
                                 nA = n_ages, 
                                 debug_bool = T)

# Set parameters
# CAUTION: rates of contact matrix are per day
coef(interMod, paste0("N_", 1:n_ages)) <- N
coef(interMod, paste0("CM_", 1:(n_ages * n_ages))) <- as.numeric(t(7 * CM))

base_pars <- coef(interMod)
init_pars <- base_pars[c("I10", "I20", "R10", "R20", "R120")]

# Check initial conditions
x0 <- rinit(interMod)
x0 <- data.frame(agecat = rownames(x0), 
                 N0 = as.numeric(x0))
print(sprintf("Model tot pop: %d, data: %d", sum(x0$N0), sum(N)))
print(sum(x0$N0)); print(sum(N))

init_vars <- c("X_IS", "X_SI", "X_RS", "X_SR", "X_RR")
for(i in seq_along(init_vars)) {
  print(sprintf("Target initial conditions: %.6f", init_pars[i]))
  print("Realized initial conditions: ")
  print(x0$N0[x0$agecat %in% paste0(init_vars[i], "_", 1:n_ages)] / N)
}

# Run simulation ----------------------------------------------------------
coef(interMod, names(base_pars)) <- unname(base_pars)
#coef(interMod, c("b1", "b2")) <- c(0.063, 0.043)
#coef(interMod, names(init_pars)) <- c(1e-4, 1e-4, 0, 0, 0)
tj <- trajectory(object = interMod, format = "data.frame")

# Add prevalence variables
for(i in 1:n_ages) {
  tj[[paste0("p1_", i)]] <- rowSums(tj[, paste0(c("X_IS_", "X_II_", "X_IT_", "X_IR_"), i)])
  tj[[paste0("p2_", i)]] <- rowSums(tj[, paste0(c("X_SI_", "X_II_", "X_TI_", "X_RI_"), i)])
}

# Pivot to long format 
tj_long <- tj %>% 
  pivot_longer(cols = -c("time", ".id"), names_to = "var_nm_full") %>% 
  mutate(var_type = ifelse(str_detect(var_nm_full, "X"), "state_var", ifelse(str_detect(var_nm_full, "H"), "obs_var", "other_var")),
         var_nm = str_remove(string = var_nm_full, pattern = "_[0-9]"),
         age_no = str_extract(string = var_nm_full, pattern = "_[0-9]")) %>%
  mutate(age_no = str_remove(age_no, "_") %>% as.factor()) %>% 
  group_by(time, .id, age_no) %>% 
  mutate(N = sum(value[var_type == "state_var"])) %>% 
  ungroup()

# Plot population size
pl <- ggplot(data = tj_long %>% select(time, age_no, N) %>% unique(), 
             mapping = aes(x = time, y = N, color = age_no)) + 
  geom_line()
print(pl)

# Plot prevalence of infection
pl <- ggplot(data = tj_long %>% filter(var_nm == "p1"), 
             mapping = aes(x = time, y = value / N, color = age_no)) + 
  geom_line()
print(pl)

#######################################################################################################
# End
#######################################################################################################
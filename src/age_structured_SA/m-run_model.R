#######################################################################################################
#  Run simulations of the age-structured model of flu-RSV interaction
#######################################################################################################

rm(list = ls())
source("s-base_packages.R")
source("f-CreateInteractionMod.R")
library(socialmixr)
library(rootSolve)
debug_bool <- T
theme_set(theme_bw())
par(bty = "l", las = 1, lwd = 2)
print(packageVersion("pomp"))

# Create contact matrix -------------------------------------------------------
# Contact data from Hong-Kong: https://zenodo.org/record/3874808
# 2021 mid-year population estimates from Hong-Kong: https://www.censtatd.gov.hk/en/web_table.html?id=1A
# Could not find pop <1 yr, assume it represents one fifth of the population aged 0-4 yr 
# Take POLYMOD matrix from the UK, corrected for reciprocity 
age_limits <- c(0, 1, 5, 16, 65) # Age limits for age categories
data(polymod)
CM_all <- contact_matrix(survey = get_survey(survey = "https://doi.org/10.5281/zenodo.3874808"), 
                         age.limits = age_limits,
                         survey.pop = data.frame(lower.age.limit = age_limits, population = c(45.8, 183.2, 578.8, 5153.8, 1451.5) * 1e3),
                         symmetric = T)
  
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

# Calculate reproduction number using next-generation matrix approach
Ri_fun <- function(beta, gamma, C_mat, tau_vec, r0_vec, N_vec) {
  NGM_mat <- matrix(data = 0, nrow = length(r0_vec), ncol = length(r0_vec))
  
  for(i in 1:nrow(NGM_mat)) {
    for(j in 1:ncol(NGM_mat)) {
      NGM_mat[i, j] <- beta * C_mat[i, j] / gamma * tau_vec[i] * (1 - r0_vec[i]) * N_vec[i] / N_vec[j]
    }
  }
  
  out <- NGM_mat %>% 
    eigen() %>% 
    getElement("values") %>% 
    max()
  return(out)
}

# Define function to find beta for a target Ri
find_beta <- function(beta, gamma, C_mat, tau_vec, r0_vec, N_vec, Ri) {
  Ri_fun(beta, gamma, C_mat, tau_vec, r0_vec, N_vec) - Ri
}

# Flu and RSV parameters
gamma_vec <- c(1 / 5, 1 / 10) # Recovery rates (per day)
tau_l <- list(rep(1, n_ages), c(1, 0.75, 0.65, 0.65, 0.65)) # Susceptibility
r0_l <- list(rep(0.5897, n_ages), rep(0.647, n_ages)) # Initial fraction recovered
R_vec <- c(1.805, 1.931) # Initial reproduction number
beta_vec <- c(0, 0)

for(i in 1:2) {
  
  # Find root of find_beta function
  tmp <- uniroot(f = find_beta, 
                 lower = 0, 
                 upper = 1, 
                 gamma = gamma_vec[i], 
                 C_mat = CM, 
                 tau_vec = tau_l[[i]],
                 r0_vec = r0_l[[i]], 
                 N_vec = N, 
                 Ri = R_vec[i])
  
  beta_vec[i] <- tmp$root
  print(sprintf("beta[%d] = %.2f", i, beta_vec[i]))
}

# Create pomp object ------------------------------------------------------
interMod <- CreateInteractionMod(dat = data.frame(time = seq(from = 0, to = 53, by = 0.1), temp = 0, ah = 0), 
                                 nA = n_ages, 
                                 debug_bool = F)

# Set parameters
# CAUTION: rates of contact matrix are per day
coef(interMod, paste0("N_", 1:n_ages)) <- N
coef(interMod, paste0("CM_", 1:(n_ages * n_ages))) <- as.numeric(t(7 * CM))
coef(interMod, c("b1", "b2")) <- beta_vec

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
#coef(interMod, c("theta_lambda1", "theta_lambda2")) <- 1
tj <- trajectory(object = interMod, 
                 format = "data.frame", 
                 ode_control = list(method = "lsoda"))

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
             mapping = aes(x = time, y = N / 1e6, color = age_no)) + 
  geom_line() + 
  labs(x = "Time (weeks)", y = "Population size (millions)", color = "Age")
print(pl)

# Plot prevalence of infection
var_to_plot <- c("H1", "H2")
pl <- ggplot(data = tj_long %>% filter(var_nm %in% var_to_plot), 
             mapping = aes(x = time, y = value / N, color = age_no)) + 
  geom_line() + 
  facet_wrap(~ var_nm) + 
  labs(x = "Time (weeks)", y = "Fraction", color = "Age", title = var_to_plot)
print(pl)

#######################################################################################################
# End
#######################################################################################################
# -----------------------------------------
# Run SEIR-SEI Vaccination Model for 1, 5, 10 years
# -----------------------------------------

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

# Define the model
diff_HV_det_logistic <- function(t, pop, para){
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Vh <- pop[5]
  Ddh <- pop[6]; Dmh <- pop[7]; Sv <- pop[8]; Ev <- pop[9]; Iv <- pop[10]; Dv <- pop[11]
  
  Nh <- Sh + Eh + Ih + Rh + Vh
  Nv <- Sv + Ev + Iv
  
  nu_t <- ifelse(t >= para$t_vax, para$nu_rate, 0)
  lambda_h <- para$f_h * para$alpha * para$p_h * Iv / Nh
  lambda_h_v <- (1 - para$epsilon_v) * lambda_h
  lambda_v <- para$f_h * para$alpha * para$p_v * Ih / Nh
  
  dSh  <- (1 - para$nu) * para$b_h * Nh - lambda_h * Sh - nu_t * Sh + para$omega * Rh + para$omega_v * Vh - para$mu_h * Sh
  dEh  <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh + lambda_h_v * Vh
  dIh  <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh  <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dVh  <- para$nu * para$b_h * Nh + nu_t * Sh - para$omega_v * Vh - lambda_h_v * Vh - para$mu_h * Vh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh + Vh)
  
  logistic_birth <- para$b_v * Nv * (1 - Nv / para$K_v)
  dSv <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv <- para$sigma_v * Ev - para$mu_v * Iv
  dDv <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dVh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# Wrapper function
ODE_SEIRSEI_logistic <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_det_logistic, parms = para, method = "ode45")
  return(as.data.frame(result))
}

# Parameters
base_para <- list(
  b_h = 1.1e-5, b_v = 0.25, f_h = 0.6, alpha = 0.8,
  p_h = 0.2, p_v = 0.5, sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  nu = 0.3, nu_rate = 0.04,
  epsilon_v = 0.7, omega_v = 1/730,
  K_v = 5000,
  t_vax = Inf # will be overwritten
)

ICs <- c(
  Sh = 500, Eh = 10, Ih = 30, Rh = 0, Vh = 0,
  Ddh = 0, Dmh = 0,
  Sv = 4000, Ev = 100, Iv = 50, Dv = 0
)

# Time horizons
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)

# Store all results
all_runs <- list()

# Loop through time durations
for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  para <- base_para
  
  # Step 1: Get baseline and first death day
  para$t_vax <- Inf
  baseline <- ODE_SEIRSEI_logistic(para, ICs, maxtime)
  first_death_day <- which(baseline$Ddh >= 1)[1] - 1
  if (is.na(first_death_day)) next
  cat(paste0("\n[", label, "] First death day: ", first_death_day, "\n"))
  
  # Step 2: Define timings
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  # Step 3: Run all scenarios
  results <- lapply(names(timings), function(scenario) {
    para$t_vax <- timings[[scenario]]
    sim <- ODE_SEIRSEI_logistic(para, ICs, maxtime)
    sim$scenario <- scenario
    sim$time_horizon <- label
    return(sim)
  })
  
  all_runs[[label]] <- do.call(rbind, results)
}

# Combine all results
all_data <- do.call(rbind, all_runs)

# ----------------------------------------
# Peak Infection Times
# ----------------------------------------
peak_infections <- all_data %>%
  group_by(scenario, time_horizon) %>%
  summarise(
    peak_Ih = max(Ih),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(peak_infections)


# Plot: Infections
ggplot(all_data, aes(x = time, y = Ih, color = scenario)) +
  geom_line() +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Infectious Humans Over Time", x = "Time (days)", y = "Ih") +
  theme_minimal()

# Plot: Cumulative Deaths
ggplot(all_data, aes(x = time, y = Ddh, color = scenario)) +
  geom_line() +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Cumulative Deaths Over Time", x = "Time (days)", y = "Deaths") +
  theme_minimal()

# SEIRV Plot
seir_long <- all_data %>%
  select(time, scenario, time_horizon, Sh, Eh, Ih, Rh, Vh) %>%
  pivot_longer(cols = c(Sh, Eh, Ih, Rh, Vh), names_to = "Compartment", values_to = "Count")

seir_long$Compartment <- factor(seir_long$Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh", "Vh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered", "Vaccinated"))

ggplot(seir_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line() +
  facet_grid(scenario ~ time_horizon) +
  labs(title = "SEIRV Dynamics by Scenario and Duration", x = "Time (days)", y = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Final death counts
final_deaths <- all_data %>%
  group_by(scenario, time_horizon) %>%
  filter(time == max(time)) %>%
  summarise(final_Ddh = round(max(Ddh)), .groups = "drop")

print(final_deaths)

# ---- SPLIT BY TIME HORIZON ----
data_1yr  <- all_data %>% filter(time_horizon == "1yr")
data_5yr  <- all_data %>% filter(time_horizon == "5yr")
data_10yr <- all_data %>% filter(time_horizon == "10yr")

# Unique time durations
horizons <- unique(all_data$time_horizon)

# Loop over each time horizon to generate plots
for (h in horizons) {
  data_subset <- all_data %>% filter(time_horizon == h)
  
  # Plot infections
  p1 <- ggplot(data_subset, aes(x = time, y = Ih, color = scenario)) +
    geom_line(size = 1.1) +
    labs(title = paste("Infectious Humans -", h),
         x = "Time (days)", y = "Ih", color = "Scenario") +
    theme_minimal()
  print(p1)
  
  # Plot deaths
  p2 <- ggplot(data_subset, aes(x = time, y = Ddh, color = scenario)) +
    geom_line(size = 1.1) +
    labs(title = paste("Cumulative Deaths -", h),
         x = "Time (days)", y = "Deaths", color = "Scenario") +
    theme_minimal()
  print(p2)
  
  # Plot SEIRV
  seir_subset <- data_subset %>%
    select(time, scenario, Sh, Eh, Ih, Rh, Vh) %>%
    pivot_longer(cols = c(Sh, Eh, Ih, Rh, Vh), names_to = "Compartment", values_to = "Count") %>%
    mutate(Compartment = factor(Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh", "Vh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered", "Vaccinated")))
  
  p3 <- ggplot(seir_subset, aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~ scenario) +
    labs(title = paste("SEIRV Dynamics -", h),
         x = "Time (days)", y = "Population", color = "Compartment") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p3)
}


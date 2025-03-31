library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Differential equations with insecticide effects ---
diff_HV_insecticide <- function(t, pop, para){
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]
  Ddh <- pop[5]; Dmh <- pop[6]; Sv <- pop[7]; Ev <- pop[8]; Iv <- pop[9]; Dv <- pop[10]
  
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  # Apply intervention only after t_intervene
  epsilon_t <- ifelse(t >= para$t_intervene, para$epsilon, 0)
  delta_v_t <- ifelse(t >= para$t_intervene, para$delta_v, 0)
  
  lambda_h <- para$f_h * (1 - epsilon_t) * para$alpha * para$p_h * Iv / Nh
  lambda_v <- para$f_h * (1 - epsilon_t) * para$alpha * para$p_v * Ih / Nh
  
  dSh <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Ih + Eh + Sh)
  
  logistic_birth <- para$b_v * Nv * (1 - Nv / para$K_v)
  dSv <- logistic_birth - lambda_v * Sv - (para$mu_v + delta_v_t) * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev - (para$mu_v + delta_v_t) * Ev
  dIv <- para$sigma_v * Ev - (para$mu_v + delta_v_t) * Iv
  dDv <- (para$mu_v + delta_v_t) * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# --- Solver wrapper ---
ODE_SEIRSEI_insecticide <- function(para, ICs, maxtime) {
  times <- seq(0, maxtime, 1)
  out <- ode(y = ICs, times = times, func = diff_HV_insecticide, parms = para, method = "ode45")
  df <- as.data.frame(out)
  colnames(df)[1] <- "time"
  return(df)
}

# --- Base parameters ---
base_para <- list(
  b_h = 1.1e-5, b_v = 0.25, f_h = 0.6, alpha = 0.8,
  p_h = 0.2, p_v = 0.5, sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  epsilon = 0.3, delta_v = 0.01,  # Intervention effect
  K_v = 5000,
  t_intervene = Inf              # Default: no intervention
)

ICs <- c(Sh = 500, Eh = 10, Ih = 30, Rh = 0,
         Ddh = 0, Dmh = 0, Sv = 4000, Ev = 100, Iv = 50, Dv = 0)

# --- Time Horizons ---
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_scenarios <- list()

# --- Simulate for each time horizon with scenario timings ---
for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  
  # Step 1: Run baseline (no intervention) to find first death
  para <- base_para
  para$t_intervene <- Inf
  baseline <- ODE_SEIRSEI_insecticide(para, ICs, maxtime)
  first_death_day <- which(baseline$Ddh >= 1)[1] - 1
  if (is.na(first_death_day)) next
  
  message(paste0("[", label, "] First disease-induced death at day: ", first_death_day))
  
  # Step 2: Define scenarios
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  # Step 3: Simulate each scenario
  for (scenario in names(timings)) {
    sim_para <- base_para
    sim_para$t_intervene <- timings[[scenario]]
    sim <- ODE_SEIRSEI_insecticide(sim_para, ICs, maxtime)
    sim$time_horizon <- label
    sim$scenario <- scenario
    all_scenarios[[paste(label, scenario, sep = "_")]] <- sim
  }
}

# Combine all runs
all_data <- bind_rows(all_scenarios)

# --- Plot 1: Infections Over Time ---
ggplot(all_data, aes(x = time, y = Ih, color = scenario)) +
  geom_line() +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Infectious Humans (Ih) Over Time with Insecticide Scenarios",
       x = "Time (days)", y = "Infectious Humans") +
  theme_minimal()

# --- Plot 2: Cumulative Deaths ---
ggplot(all_data, aes(x = time, y = Ddh, color = scenario)) +
  geom_line() +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Cumulative Deaths with Insecticide Scenarios",
       x = "Time (days)", y = "Deaths") +
  theme_minimal()

# --- Plot 3: SEIR Dynamics ---
seir_long <- all_data %>%
  select(time, scenario, time_horizon, Sh, Eh, Ih, Rh) %>%
  pivot_longer(cols = c(Sh, Eh, Ih, Rh), names_to = "Compartment", values_to = "Count")

seir_long$Compartment <- factor(seir_long$Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered"))

ggplot(seir_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line() +
  facet_grid(scenario ~ time_horizon, scales = "free_x") +
  labs(title = "SEIR Dynamics by Timing and Horizon",
       x = "Time (days)", y = "Human Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

# --- Final Summary Stats ---
final_summary <- all_data %>%
  group_by(time_horizon, scenario) %>%
  filter(time == max(time)) %>%
  summarise(
    Final_Ih = round(Ih),
    Final_Deaths = round(Ddh),
    Peak_Ih = round(max(Ih)),
    Peak_Day = time[which.max(Ih)],
    .groups = "drop"
  )

print(final_summary)

unique_horizons <- unique(all_data$time_horizon)

for (h in unique_horizons) {
  cat("Plotting for:", h, "\n")
  
  data_h <- all_data %>% filter(time_horizon == h)
  
  # 1. Infections over time
  p1 <- ggplot(data_h, aes(x = time, y = Ih, color = scenario)) +
    geom_line(size = 1) +
    labs(title = paste("Infectious Humans Over Time (", h, ")", sep = ""),
         x = "Time (days)", y = "Ih") +
    theme_minimal()
  print(p1)
  
  # 2. Cumulative deaths
  p2 <- ggplot(data_h, aes(x = time, y = Ddh, color = scenario)) +
    geom_line(size = 1) +
    labs(title = paste("Cumulative Deaths Over Time (", h, ")", sep = ""),
         x = "Time (days)", y = "Deaths") +
    theme_minimal()
  print(p2)
  
  # 3. SEIR Dynamics
  seir_h <- data_h %>%
    select(time, scenario, Sh, Eh, Ih, Rh) %>%
    pivot_longer(cols = c("Sh", "Eh", "Ih", "Rh"),
                 names_to = "Compartment", values_to = "Count")
  
  seir_h$Compartment <- factor(seir_h$Compartment,
                               levels = c("Sh", "Eh", "Ih", "Rh"),
                               labels = c("Susceptible", "Exposed", "Infectious", "Recovered"))
  
  p3 <- ggplot(seir_h, aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~ scenario) +
    labs(title = paste("SEIR Dynamics (", h, ")", sep = ""),
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p3)
}

# --- Final Summary Stats ---
summary_table <- all_data %>%
  group_by(time_horizon, scenario) %>%
  summarise(
    final_deaths = round(max(Ddh), 2),
    peak_Ih = round(max(Ih), 2),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(summary_table)

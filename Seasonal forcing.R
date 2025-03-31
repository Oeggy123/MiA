#########################################################
# MiA - Seasonal Forcing Malaria Model with Logistic Growth
#########################################################

library(deSolve)
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)

# Seasonal forcing parameters
T_season <- 365
phi <- 0
A <- 0.3     # Amplitude of biting rate variation
Ab <- 0.3    # Amplitude of vector birth rate variation

# Logistic + Seasonal ODE system
seasonal_logistic_model <- function(t, pop, para){
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Ddh <- pop[5]; Dmh <- pop[6]
  Sv <- pop[7]; Ev <- pop[8]; Iv <- pop[9]; Dv <- pop[10]
  
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  alpha_t <- para$alpha_bar * (1 + A * cos(2 * pi * t / T_season + phi))
  bv_t <- para$b_v_bar * (1 + Ab * cos(2 * pi * t / T_season + phi))
  
  lambda_h <- para$f_h * alpha_t * para$p_h * Iv / Nh
  lambda_v <- para$f_h * alpha_t * para$p_v * Ih / Nh
  
  dSh <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh)
  
  logistic_birth <- bv_t * Nv * (1 - Nv / para$K_v)
  dSv <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv <- para$sigma_v * Ev - para$mu_v * Iv
  dDv <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# ODE runner
run_seasonal_logistic_model <- function(para, ICs, maxtime){
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = seasonal_logistic_model, parms = para, method = "ode45")
  return(as.data.frame(result))
}

# Base parameters
base_para <- list(
  b_h = 1.1e-5, b_v_bar = 0.25, f_h = 0.6,
  alpha_bar = 0.8, p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  K_v = 5000
)

# Initial state
ICs <- c(
  Sh = 500, Eh = 10, Ih = 30, Rh = 0, Ddh = 0, Dmh = 0,
  Sv = 4000, Ev = 100, Iv = 50, Dv = 0
)

# Time horizons
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_results <- list()

# Run simulations
for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  para <- base_para
  sim <- run_seasonal_logistic_model(para, ICs, maxtime)
  sim$time_horizon <- label
  all_results[[label]] <- sim
}

# Combine into one dataframe
all_data <- bind_rows(all_results)

# ----------------------------
# Plot: Infections Over Time
# ----------------------------
ggplot(all_data, aes(x = time, y = Ih)) +
  geom_line(color = "red", linewidth = 1) +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Infectious Humans (Ih) with Seasonal Logistic Model",
       x = "Time (days)", y = "Infectious Humans") +
  theme_minimal()

# ----------------------------
# Plot: Cumulative Deaths
# ----------------------------
ggplot(all_data, aes(x = time, y = Ddh)) +
  geom_line(color = "darkred", linewidth = 1) +
  facet_wrap(~ time_horizon, scales = "free_x") +
  labs(title = "Cumulative Disease-Induced Deaths", x = "Time (days)", y = "Deaths") +
  theme_minimal()

# ----------------------------
# SEIR Dynamics
# ----------------------------
seir_long <- all_data %>%
  select(time, time_horizon, Sh, Eh, Ih, Rh) %>%
  pivot_longer(cols = c(Sh, Eh, Ih, Rh), names_to = "Compartment", values_to = "Count")

seir_long$Compartment <- factor(seir_long$Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered"))

ggplot(seir_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line() +
  facet_wrap(~ time_horizon) +
  labs(title = "SEIR Dynamics with Seasonal Logistic Growth", x = "Time (days)", y = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ----------------------------
# Summary Stats
# ----------------------------
final_stats <- all_data %>%
  group_by(time_horizon) %>%
  summarise(
    Final_Deaths = round(tail(Ddh, 1)),
    Peak_Ih = round(max(Ih)),
    Peak_Day = time[which.max(Ih)],
    .groups = "drop"
  )

print(final_stats)

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------
# Logistic SEIR-SEI Model with Bed Nets
# ------------------------------------------

diff_HV_bednets_logistic <- function(t, pop, para){
  # State variables
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]
  Ddh <- pop[5]; Dmh <- pop[6]; Sv <- pop[7]; Ev <- pop[8]; Iv <- pop[9]; Dv <- pop[10]
  
  # Totals
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  # Time-dependent bed net efficacy
  eta_t <- ifelse(t >= para$t_bednets, para$eta, 0)
  alpha_eff <- (1 - eta_t) * para$alpha
  
  # Force of infection
  lambda_h <- para$f_h * alpha_eff * para$p_h * Iv / Nh
  lambda_v <- para$f_h * alpha_eff * para$p_v * Ih / Nh
  
  # Human equations
  dSh  <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh  <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh  <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh  <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh)
  
  # Vector equations (logistic birth)
  logistic_birth <- para$b_v * Nv * (1 - Nv / para$K_v)
  dSv  <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv  <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv  <- para$sigma_v * Ev - para$mu_v * Iv
  dDv  <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# ODE solver
ODE_SEIRSEI_bednets_logistic <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_bednets_logistic, parms = para, method = "ode45")
  return(as.data.frame(result))
}

# ------------------------------------------
# Parameters and Initial Setup
# ------------------------------------------
para <- list(
  b_h = 1.1e-5, b_v = 0.25, f_h = 0.6, alpha = 0.8,
  p_h = 0.2, p_v = 0.5, sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  eta = 0.6,         # bed net efficacy
  t_bednets = Inf,   # intervention day (to be updated in scenarios)
  K_v = 5000         # carrying capacity
)

ICs <- c(
  Sh = 500, Eh = 10, Ih = 30, Rh = 0,
  Ddh = 0, Dmh = 0,
  Sv = 4000, Ev = 100, Iv = 50, Dv = 0
)

maxtime <- 2000

# ------------------------------------------
# Step 1: Run no-intervention to find first death
# ------------------------------------------
para$t_bednets <- Inf
baseline <- ODE_SEIRSEI_bednets_logistic(para, ICs, maxtime)
first_death_day <- which(baseline$Ddh >= 1)[1] - 1
cat("First disease-induced death occurs at day:", first_death_day, "\n")

# ------------------------------------------
# Step 2: Define timings
# ------------------------------------------
timings <- list(
  "Early (D - 30)" = max(0, first_death_day - 30),
  "On First Death (D)" = first_death_day,
  "Delayed (D + 30)" = first_death_day + 30,
  "Very Delayed (D + 90)" = first_death_day + 90
)

# ------------------------------------------
# Step 3: Run simulations
# ------------------------------------------
run_bednet_sim <- function(t_bednets_time, label) {
  para$t_bednets <- t_bednets_time
  sim <- ODE_SEIRSEI_bednets_logistic(para, ICs, maxtime)
  sim$scenario <- label
  return(sim)
}

scenarios <- lapply(names(timings), function(name) {
  run_bednet_sim(timings[[name]], name)
})

all_bednet_scenarios <- do.call(rbind, scenarios)

# ------------------------------------------
# Step 4: Plot infections over time
# ------------------------------------------
ggplot(all_bednet_scenarios, aes(x = time, y = Ih, color = scenario)) +
  geom_line(size = 1) +
  labs(title = "Impact of Bed Net Timing on Infections",
       x = "Time (days)", y = "Infectious Humans (Ih)", color = "Scenario") +
  theme_minimal()

# ------------------------------------------
# Step 5: Plot cumulative deaths
# ------------------------------------------
ggplot(all_bednet_scenarios, aes(x = time, y = Ddh, color = scenario)) +
  geom_line(size = 1) +
  labs(title = "Cumulative Disease-Induced Deaths",
       x = "Time (days)", y = "Deaths", color = "Scenario") +
  theme_minimal()

# ------------------------------------------
# Step 6: Plot full SEIR dynamics
# ------------------------------------------
seir_long <- all_bednet_scenarios %>%
  select(time, scenario, Sh, Eh, Ih, Rh) %>%
  pivot_longer(cols = c("Sh", "Eh", "Ih", "Rh"), names_to = "Compartment", values_to = "Count")

seir_long$Compartment <- factor(seir_long$Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh"),
                                labels = c("Susceptible (S)", "Exposed (E)", 
                                           "Infectious (I)", "Recovered (R)"))

ggplot(seir_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line(size = 1) +
  facet_wrap(~ scenario, ncol = 2) +
  labs(title = "SEIR Dynamics with Bed Net Interventions",
       x = "Time (days)", y = "Human Population",
       color = "Compartment") +
  theme_minimal() +
  theme(legend.position = "bottom")

# ------------------------------------------
# Step 7: Final cumulative deaths
# ------------------------------------------
end_deaths_bednets <- aggregate(Ddh ~ scenario, data = all_bednet_scenarios[all_bednet_scenarios$time == maxtime, ], max)
print(end_deaths_bednets)


# ------------------------------------------
# Time Horizons: 1 year, 5 years, 10 years
# ------------------------------------------
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_runs <- list()

for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  para$t_bednets <- Inf
  
  # Run baseline to find first disease-induced death
  baseline <- ODE_SEIRSEI_bednets_logistic(para, ICs, maxtime)
  first_death_day <- which(baseline$Ddh >= 1)[1] - 1
  if (is.na(first_death_day)) next
  cat(paste0("\n[", label, "] First death day: ", first_death_day, "\n"))
  
  # Define scenario timings
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  # Run scenarios
  results <- lapply(names(timings), function(scenario) {
    para$t_bednets <- timings[[scenario]]
    sim <- ODE_SEIRSEI_bednets_logistic(para, ICs, maxtime)
    sim$scenario <- scenario
    sim$time_horizon <- label
    return(sim)
  })
  
  all_runs[[label]] <- do.call(rbind, results)
}

# Combine all
all_bednet_data <- do.call(rbind, all_runs)


# Loop over each horizon to plot separately
for (h in unique(all_bednet_data$time_horizon)) {
  df <- filter(all_bednet_data, time_horizon == h)
  
  # Infections
  print(
    ggplot(df, aes(x = time, y = Ih, color = scenario)) +
      geom_line(size = 1.2) +
      labs(title = paste("Infectious Humans -", h),
           x = "Time (days)", y = "Infectious (Ih)") +
      theme_minimal()
  )
  
  # Deaths
  print(
    ggplot(df, aes(x = time, y = Ddh, color = scenario)) +
      geom_line(size = 1.2) +
      labs(title = paste("Cumulative Deaths -", h),
           x = "Time (days)", y = "Cumulative Deaths") +
      theme_minimal()
  )
  
  # SEIR dynamics
  seir_df <- df %>%
    select(time, scenario, Sh, Eh, Ih, Rh) %>%
    pivot_longer(cols = c("Sh", "Eh", "Ih", "Rh"), names_to = "Compartment", values_to = "Count") %>%
    mutate(Compartment = factor(Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered")))
  
  print(
    ggplot(seir_df, aes(x = time, y = Count, color = Compartment)) +
      geom_line() +
      facet_wrap(~ scenario) +
      labs(title = paste("SEIR Dynamics -", h),
           x = "Time (days)", y = "Population") +
      theme_minimal() +
      theme(legend.position = "bottom")
  )
}

# ------------------------------------------
# Final Summary: Deaths, Peak Infection, Day of Peak
# ------------------------------------------
summary_table <- all_bednet_data %>%
  group_by(scenario, time_horizon) %>%
  summarise(
    final_deaths = round(max(Ddh), 2),
    peak_Ih = round(max(Ih), 2),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(summary_table)

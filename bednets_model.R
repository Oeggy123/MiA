library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# Bed net seasonal logistic model
diff_HV_bednet_seasonal <- function(t, pop, para){
  alpha_t <- para$alpha_bar * (1 + para$A * cos(2 * pi * t / para$T + para$phi))
  b_v_t <- para$b_v_bar * (1 + para$A_b * cos(2 * pi * t / para$T + para$phi_b))
  
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Vh <- pop[5]
  Ddh <- pop[6]; Dmh <- pop[7]; Sv <- pop[8]; Ev <- pop[9]; Iv <- pop[10]; Dv <- pop[11]
  Nh <- Sh + Eh + Ih + Rh + Vh
  Nv <- Sv + Ev + Iv
  
  eta_t <- ifelse(t >= para$t_bednets, para$eta, 0)
  alpha_eff <- (1 - eta_t) * alpha_t
  
  lambda_h <- para$f_h * alpha_eff * para$p_h * Iv / Nh
  lambda_v <- para$f_h * alpha_eff * para$p_v * Ih / Nh
  
  logistic_birth <- b_v_t * Nv * (1 - Nv / para$K_v)
  
  dSh  <- (1 - para$nu) * para$b_h * Nh + para$omega * Rh + para$omega_v * Vh - lambda_h * Sh - para$mu_h * Sh
  dEh  <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh + (1 - para$epsilon_v) * lambda_h * Vh
  dIh  <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh  <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dVh  <- para$nu * para$b_h * Nh - para$omega_v * Vh - (1 - para$epsilon_v) * lambda_h * Vh - para$mu_h * Vh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Ih + Eh + Sh + Vh + Rh)
  dSv  <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv  <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv  <- para$sigma_v * Ev - para$mu_v * Iv
  dDv  <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dVh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

ODE_model <- function(para, ICs, maxtime){
  t_seq <- seq(0, maxtime, 1)
  out <- ode(y = ICs, times = t_seq, func = diff_HV_bednet_seasonal, parms = para, method = "ode45")
  as.data.frame(out)
}

# Parameters
base_para <- list(
  b_h = 1.1e-5, b_v_bar = 0.25, f_h = 0.6,
  alpha_bar = 0.8, p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  nu = 0, epsilon_v = 0, omega_v = 0,
  A = 0.3, T = 365, phi = 0,
  A_b = 0.2, phi_b = 0,
  eta = 0.5,
  t_bednets = Inf,
  K_v = 5000
)

ICs <- c(Sh = 500, Eh = 10, Ih = 30, Rh = 0, Vh = 0,
         Ddh = 0, Dmh = 0,
         Sv = 4000, Ev = 100, Iv = 50, Dv = 0)

time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_runs <- list()

for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  para <- base_para
  
  # Baseline run to find first death
  para$t_bednets <- Inf
  baseline <- ODE_model(para, ICs, maxtime)
  first_death_day <- which(baseline$Ddh >= 1)[1] - 1
  if (is.na(first_death_day)) next
  
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  scenarios <- lapply(names(timings), function(name){
    para$t_bednets <- timings[[name]]
    sim <- ODE_model(para, ICs, maxtime)
    sim$scenario <- name
    sim$time_horizon <- label
    sim
  })
  
  all_runs[[label]] <- do.call(rbind, scenarios)
}

all_bednet_data <- do.call(rbind, all_runs)

# Plot SEIRV by scenario and horizon
seir_long <- all_bednet_data %>%
  select(time, time_horizon, scenario, Sh, Eh, Ih, Rh) %>%
  pivot_longer(cols = c(Sh, Eh, Ih, Rh), names_to = "Compartment", values_to = "Count")

seir_long$Compartment <- factor(seir_long$Compartment,
                                levels = c("Sh", "Eh", "Ih", "Rh"),
                                labels = c("Susceptible", "Exposed", "Infectious", "Recovered"))

# Plot SEIRV dynamics separately for each time horizon
for (horizon in unique(seir_long$time_horizon)) {
  data_subset <- seir_long %>% filter(time_horizon == horizon)
  
  p <- ggplot(data_subset, aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~ scenario, ncol = 2, scales = "free_x") +
    labs(title = paste("Bed Net Scenarios with Seasonal Forcing -", horizon),
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)  # Show the plot 
}


summary_stats <- all_bednet_data %>%
  group_by(time_horizon, scenario) %>%
  summarise(
    final_Ddh = round(max(Ddh), 2),
    peak_Ih = round(max(Ih), 2),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(summary_stats)

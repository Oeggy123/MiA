library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------
# Seasonal Functions
# ----------------------------------------
alpha_t <- function(t, alpha_bar = 0.5, A = 0.3, T = 365, phi = 0) {
  alpha_bar * (1 + A * cos((2 * pi * t / T) + phi))
}

bv_t <- function(t, bv_bar = 0.13, A = 0.3, T = 365, phi = 0) {
  bv_bar * (1 + A * cos((2 * pi * t / T) + phi))
}

# ----------------------------------------
# Differential Equation System 
# ----------------------------------------
diff_chemo_migration <- function(t, pop, para) {
  with(as.list(c(pop, para)), {
    
    Nh <- Sh + Eh + Ih + Rh + Ph
    Nv <- Sv + Ev + Iv
    
    alpha <- alpha_t(t, alpha_bar = para$alpha_bar, A = para$A, T = para$T, phi = para$phi)
    bv <- bv_t(t, bv_bar = para$b_v_bar, A = para$A_b, T = para$T, phi = para$phi_b)
    
    lambda_h <- f_h * alpha * p_h * Iv / Nh
    lambda_v <- f_h * alpha * p_v * Ih / Nh
    
    psi_t <- ifelse(t >= para$t_chemo, para$psi_rate, 0)
    logistic_birth <- bv * Nv * (1 - Nv / para$K_v)
    
    dSh <- b_h * Nh - lambda_h * Sh + omega * Rh - mu_h * Sh -
      Mh_out * Sh + Mh_in * Sh_ext + omega_p * Ph
    
    dEh <- lambda_h * (1 - rho) * Sh - sigma_h * Eh - mu_h * Eh -
      Mh_out * Eh + Mh_in * Eh_ext
    
    dIh <- sigma_h * Eh - gamma * Ih - (mu_h + d + xi) * Ih -
      Mh_out * Ih + Mh_in * Ih_ext
    
    dRh <- gamma * Ih + xi * Ih - omega * Rh - mu_h * Rh -
      Mh_out * Rh + Mh_in * Rh_ext
    
    dPh <- psi_t * Mh_in * Sh_ext - lambda_h * (1 - rho) * Ph - mu_h * Ph - omega_p * Ph
    
    dDdh <- d * Ih
    
    dSv <- logistic_birth - lambda_v * Sv - mu_v * Sv +
      Mv_in * Sv_ext - Mv_out * Sv
    
    dEv <- lambda_v * Sv - sigma_v * Ev - mu_v * Ev +
      Mv_in * Ev_ext - Mv_out * Ev
    
    dIv <- sigma_v * Ev - mu_v * Iv +
      Mv_in * Iv_ext - Mv_out * Iv
    
    return(list(c(dSh, dEh, dIh, dRh, dPh, dDdh, dSv, dEv, dIv)))
  })
}

# ----------------------------------------
# Initial Conditions (add Ddh = 0)
# ----------------------------------------
ICs <- c(
  Sh = 500, Eh = 10, Ih = 20, Rh = 0, Ph = 0, Ddh = 0,
  Sv = 4000, Ev = 100, Iv = 50
)


# ODE solver
run_model <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  out <- ode(y = ICs, times = t_seq, func = diff_chemo_migration, parms = para, method = "ode45")
  as.data.frame(out)
}

# ----------------------------------------
# Parameters and Initial Conditions
# ----------------------------------------
base_para <- list(
  b_h = 1.1e-5, b_v_bar = 0.25,
  f_h = 0.6,                        
  alpha_bar = 0.8, p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  
  Mh_in = 0.2,    # ⬆️ heavy inflow
  Mh_out = 0.2,   # ⬆️ heavy outflow
  
  psi = 0.8, psi_rate = 0.8, t_chemo = Inf,
  rho = 0.7, omega_p = 0.001, xi = 0.005,
  
  Mv_in = 0.02, Mv_out = 0.02,
  Sh_ext = 500, Eh_ext = 50, Ih_ext = 20, Rh_ext = 100,
  Sv_ext = 100, Ev_ext = 10, Iv_ext = 5,
  
  A = 0.3, A_b = 0.2, T = 365, phi = 0, phi_b = 0,
  K_v = 5000
)


time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_runs <- list()

# ----------------------------------------
# Loop Through Time Horizons
# ----------------------------------------
for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  para <- base_para
  
  # Baseline run
  para$t_chemo <- Inf
  baseline <- run_model(para, ICs, maxtime)
  first_death_day <- which(baseline$Ih >= 1)[1] - 1
  if (is.na(first_death_day)) next
  
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  scenarios <- lapply(names(timings), function(name) {
    para$t_chemo <- timings[[name]]
    sim <- run_model(para, ICs, maxtime)
    sim$scenario <- name
    sim$time_horizon <- label
    sim
  })
  
  all_runs[[label]] <- do.call(rbind, scenarios)
}

all_chemo_data <- do.call(rbind, all_runs)

# ----------------------------------------
# Plotting
# ----------------------------------------
chemo_long <- all_chemo_data %>%
  select(time, time_horizon, scenario, Sh, Eh, Ih, Rh) %>%
  pivot_longer(cols = c(Sh, Eh, Ih, Rh), names_to = "Compartment", values_to = "Count") %>%
  mutate(Compartment = factor(Compartment, levels = c("Sh", "Eh", "Ih", "Rh")))

for (h in unique(chemo_long$time_horizon)) {
  df <- chemo_long %>% filter(time_horizon == h)
  
  p <- ggplot(df, aes(x = time, y = Count, color = Compartment)) +
    geom_line() +
    facet_wrap(~ scenario, scales = "free_x") +
    labs(title = paste("Chemo + Migration (", h, ")", sep = ""),
         x = "Time (days)", y = "Population") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(p)
  ggsave(paste0("RFigs/chemo_migration_seasonal_", h, ".pdf"), plot = p, width = 8, height = 6)
}

# ----------------------------------------
# Final Stats
# ----------------------------------------
summary_stats <- all_chemo_data %>%
  group_by(time_horizon, scenario) %>%
  summarise(
    final_Ddh = round(max(Ddh), 2),
    peak_Ih = round(max(Ih), 2),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(summary_stats)

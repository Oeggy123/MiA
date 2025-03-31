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

bv_t <- function(t, b_v_bar = 0.25, A = 0.2, T = 365, phi = 0) {
  b_v_bar * (1 + A * cos((2 * pi * t / T) + phi))
}

# ----------------------------------------
# Seasonal Logistic Metapopulation Model
# ----------------------------------------
metapop_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    np <- length(patches)
    
    Sh  <- state[patch_indices("Sh")]
    Eh  <- state[patch_indices("Eh")]
    Ih  <- state[patch_indices("Ih")]
    Rh  <- state[patch_indices("Rh")]
    Vh  <- state[patch_indices("Vh")]
    Ddh <- state[patch_indices("Ddh")]
    Dmh <- state[patch_indices("Dmh")]
    
    Sv  <- state[patch_indices("Sv")]
    Ev  <- state[patch_indices("Ev")]
    Iv  <- state[patch_indices("Iv")]
    Dv  <- state[patch_indices("Dv")]
    
    Nh <- Sh + Eh + Ih + Rh + Vh
    Nv <- Sv + Ev + Iv
    
    # Seasonal forcing
    alpha <- alpha_t(t, alpha_bar, A, T, phi)
    b_v   <- bv_t(t, b_v_bar, A_b, T, phi_b)
    
    # Vaccination trigger
    nu_t <- ifelse(t >= t_intv, nu_rate, 0)
    
    lambda_h  <- f_h * alpha * p_h * Iv / (Nh)
    lambda_hv <- (1 - epsilon_v) * lambda_h
    lambda_v  <- f_h * alpha * p_v * Ih / (Nh)
    
    logistic_birth <- b_v * Nv * (1 - Nv / K_v)
    
    # Humans
    dSh  <- (1 - nu) * b_h * Nh - lambda_h * Sh - nu_t * Sh + omega * Rh + omega_v * Vh - mu_h * Sh + inflow(Sh, M) - outflow(Sh, M)
    dEh  <- lambda_h * Sh + lambda_hv * Vh - sigma_h * Eh - mu_h * Eh + inflow(Eh, M) - outflow(Eh, M)
    dIh  <- sigma_h * Eh - gamma * Ih - (mu_h + d) * Ih + inflow(Ih, M) - outflow(Ih, M)
    dRh  <- gamma * Ih - omega * Rh - mu_h * Rh + inflow(Rh, M) - outflow(Rh, M)
    dVh  <- nu * b_h * Nh + nu_t * Sh - omega_v * Vh - lambda_hv * Vh - mu_h * Vh + inflow(Vh, M) - outflow(Vh, M)
    dDdh <- d * Ih
    dDmh <- mu_h * (Sh + Eh + Ih + Rh + Vh)
    
    # Vectors
    dSv <- logistic_birth - lambda_v * Sv - mu_v * Sv + inflow(Sv, Mv) - outflow(Sv, Mv)
    dEv <- lambda_v * Sv - sigma_v * Ev - mu_v * Ev + inflow(Ev, Mv) - outflow(Ev, Mv)
    dIv <- sigma_v * Ev - mu_v * Iv + inflow(Iv, Mv) - outflow(Iv, Mv)
    dDv <- mu_v * (Sv + Ev + Iv)
    
    return(list(c(dSh, dEh, dIh, dRh, dVh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
  })
}




patch_indices <- function(var) {
  start <- match(var, c("Sh", "Eh", "Ih", "Rh", "Vh", "Ddh", "Dmh", "Sv", "Ev", "Iv", "Dv"))
  return(((start - 1) * np + 1):(start * np))
}


inflow <- function(vec, M) as.vector(t(M) %*% vec)
outflow <- function(vec, M) as.vector(M %*% vec)

# ----------------------------------------
# Settings
# ----------------------------------------
np <- 3
patches <- 1:np
M <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.90
), nrow = np, byrow = TRUE)
Mv <- M

base_parameters <- list(
  b_h = 1.1e-5, b_v_bar = 0.25,
  f_h = 0.6, alpha_bar = 0.8,
  p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009,
  mu_h = 1.1e-5, mu_v = 0.1,
  d = 0.003, nu_rate = 0.04, epsilon_v = 0.8,
  nu = 0.3,
  omega_v = 1/730,
  K_v = 5000,
  A = 0.3, A_b = 0.2,
  T = 365, phi = 0, phi_b = 0,
  t_intv = Inf,
  np = np, patches = patches,
  M = M, Mv = Mv
)

init <- rep(0, 11 * np)
names(init) <- c(
  paste0("Sh", 1:np), paste0("Eh", 1:np), paste0("Ih", 1:np),
  paste0("Rh", 1:np), paste0("Vh", 1:np),
  paste0("Ddh", 1:np), paste0("Dmh", 1:np),
  paste0("Sv", 1:np), paste0("Ev", 1:np), paste0("Iv", 1:np),
  paste0("Dv", 1:np)
)


# Set initial values
init[c("Sh1", "Eh1", "Ih1", "Rh1", "Sv1", "Ev1", "Iv1")] <- c(1000, 20, 30, 0, 6000, 150, 100)
init[c("Sh2", "Eh2", "Ih2", "Rh2", "Sv2", "Ev2", "Iv2")] <- c(600, 0, 0, 0, 4000, 80, 40)
init[c("Sh3", "Eh3", "Ih3", "Rh3", "Sv3", "Ev3", "Iv3")] <- c(300, 0, 0, 0, 2000, 50, 25)

# Initialize new compartments to zero
init[paste0("Vh", 1:np)]   <- 0
init[paste0("Ddh", 1:np)]  <- 0
init[paste0("Dmh", 1:np)]  <- 0
init[paste0("Dv", 1:np)]   <- 0


# ----------------------------------------
# Simulation Loop
# ----------------------------------------
all_runs <- list()
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)

for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  times <- seq(0, maxtime, 1)
  para <- base_parameters
  para$t_intv <- Inf
  baseline <- ode(y = init, times = times, func = metapop_model, parms = para, method = "ode45")
  baseline_df <- as.data.frame(baseline)
  first_infection_day <- which(baseline_df$Ih1 >= 1)[1] - 1
  if (is.na(first_infection_day)) next
  
  timings <- list(
    "Early (D - 30)" = max(0, first_infection_day - 30),
    "On First Infection (D)" = first_infection_day,
    "Delayed (D + 30)" = first_infection_day + 30,
    "Very Delayed (D + 90)" = first_infection_day + 90
  )
  
  scenarios <- lapply(names(timings), function(name) {
    para$t_intv <- timings[[name]]
    sim <- ode(y = init, times = times, func = metapop_model, parms = para, method = "ode45")
    df <- as.data.frame(sim)
    df$scenario <- name
    df$time_horizon <- label
    df
  })
  
  all_runs[[label]] <- do.call(rbind, scenarios)
}

meta_df <- do.call(rbind, all_runs)

# ----------------------------------------
# Plotting
# ----------------------------------------
meta_long <- meta_df %>%
  select(time, scenario, time_horizon, Ih1, Ih2, Ih3) %>%
  pivot_longer(cols = starts_with("Ih"), names_to = "Patch", values_to = "Infected") %>%
  mutate(Patch = factor(Patch, levels = c("Ih1", "Ih2", "Ih3")))

for (h in unique(meta_long$time_horizon)) {
  df <- meta_long %>% filter(time_horizon == h)
  p <- ggplot(df, aes(x = time, y = Infected, color = Patch)) +
    geom_line() +
    facet_wrap(~ scenario, scales = "free_x") +
    labs(title = paste("Metapopulation Infection -", h),
         x = "Time (days)", y = "Infected Humans") +
    theme_minimal()
  print(p)
  ggsave(paste0("RFigs/meta_infection_", h, ".pdf"), p, width = 8, height = 6)
}

# ----------------------------------------
# Summary Stats
# ----------------------------------------
summary_stats <- meta_df %>%
  group_by(scenario, time_horizon) %>%
  summarise(
    total_deaths_patch1 = round(max(Ddh1)),
    total_deaths_patch2 = round(max(Ddh2)),
    total_deaths_patch3 = round(max(Ddh3)),
    
    peak_day_patch1 = time[which.max(Ih1)],
    peak_day_patch2 = time[which.max(Ih2)],
    peak_day_patch3 = time[which.max(Ih3)],
    
    peak_Ih_patch1 = round(max(Ih1)),
    peak_Ih_patch2 = round(max(Ih2)),
    peak_Ih_patch3 = round(max(Ih3)),
    
    first_death_patch1 = min(time[which(Ddh1 >= 1)]),
    first_death_patch2 = min(time[which(Ddh2 >= 1)]),
    first_death_patch3 = min(time[which(Ddh3 >= 1)]),
    
    .groups = "drop"
  )


print(summary_stats)

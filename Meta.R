library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------
# Seasonal Logistic Metapopulation Model
# ----------------------------------------
metapop_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    np <- length(patches)
    
    Sh <- state[patch_indices("Sh")]
    Eh <- state[patch_indices("Eh")]
    Ih <- state[patch_indices("Ih")]
    Rh <- state[patch_indices("Rh")]
    Sv <- state[patch_indices("Sv")]
    Ev <- state[patch_indices("Ev")]
    Iv <- state[patch_indices("Iv")]
    
    Nh <- Sh + Eh + Ih + Rh
    Nv <- Sv + Ev + Iv
    
    nu_t <- ifelse(t >= t_intv, nu_rate, 0)
    
    lambda_h <- f_h * alpha * p_h * Iv / Nh
    lambda_v <- f_h * alpha * p_v * Ih / Nh
    
    logistic_birth <- b_v * Nv * (1 - Nv / K_v)
    
    dSh <- b_h * Nh - lambda_h * Sh - nu_t * Sh + omega * Rh - mu_h * Sh + inflow(Sh, M) - outflow(Sh, M)
    dEh <- lambda_h * (Sh + (1 - epsilon_v) * nu_t * Sh) - sigma_h * Eh - mu_h * Eh + inflow(Eh, M) - outflow(Eh, M)
    dIh <- sigma_h * Eh - gamma * Ih - (mu_h + d) * Ih + inflow(Ih, M) - outflow(Ih, M)
    dRh <- gamma * Ih - omega * Rh - mu_h * Rh + inflow(Rh, M) - outflow(Rh, M)
    
    dSv <- logistic_birth - lambda_v * Sv - mu_v * Sv + inflow(Sv, Mv) - outflow(Sv, Mv)
    dEv <- lambda_v * Sv - sigma_v * Ev - mu_v * Ev + inflow(Ev, Mv) - outflow(Ev, Mv)
    dIv <- sigma_v * Ev - mu_v * Iv + inflow(Iv, Mv) - outflow(Iv, Mv)
    
    return(list(c(dSh, dEh, dIh, dRh, dSv, dEv, dIv)))
  })
}

# Indexing helpers
patch_indices <- function(var) {
  start <- match(var, c("Sh", "Eh", "Ih", "Rh", "Sv", "Ev", "Iv"))
  return(((start - 1) * np + 1):(start * np))
}

inflow <- function(vec, M) as.vector(t(M) %*% vec)
outflow <- function(vec, M) as.vector(M %*% vec)

# Settings
np <- 3
patches <- 1:np
M <- matrix(c(0.9, 0.05, 0.05, 0.05, 0.9, 0.05, 0.05, 0.05, 0.9), nrow = np, byrow = TRUE)
Mv <- M

# Parameters (consistent with previous models)
base_parameters <- list(
  b_h = 1.1e-5, b_v = 0.25,
  f_h = 0.3, alpha = 0.5,
  p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 0.01, omega = 0.0009,
  mu_h = 1.1e-5, mu_v = 0.1,
  d = 0.001, nu_rate = 0.04, epsilon_v = 0.8,
  K_v = 5000,
  t_intv = Inf,
  np = np, patches = patches,
  M = M, Mv = Mv
)

# Initial conditions (matching previous models)
init <- rep(0, 7 * np)
names(init) <- c(paste0("Sh", 1:np), paste0("Eh", 1:np), paste0("Ih", 1:np), paste0("Rh", 1:np),
                 paste0("Sv", 1:np), paste0("Ev", 1:np), paste0("Iv", 1:np))

# Set default values across all patches
init[paste0("Sh", 1:np)] <- 500
init[paste0("Eh", 1:np)] <- 0
init[paste0("Ih", 1:np)] <- 0
init[paste0("Rh", 1:np)] <- 0
init[paste0("Sv", 1:np)] <- 4000
init[paste0("Ev", 1:np)] <- 100
init[paste0("Iv", 1:np)] <- 50

# Seed Patch 1 with infections (10 exposed, 20 infected)
init["Eh1"] <- 10
init["Ih1"] <- 20
init["Sh1"] <- init["Sh1"] - 30  # Adjust to conserve population

# Time horizons
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_runs <- list()

# Simulation loop
for (label in names(time_horizons)) {
  maxtime <- time_horizons[[label]]
  times <- seq(0, maxtime, 1)
  para <- base_parameters
  
  para$t_intv <- Inf
  baseline <- ode(y = init, times = times, func = metapop_model, parms = para)
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
    sim <- ode(y = init, times = times, func = metapop_model, parms = para)
    df <- as.data.frame(sim)
    df$scenario <- name
    df$time_horizon <- label
    df
  })
  
  all_runs[[label]] <- do.call(rbind, scenarios)
}

meta_df <- do.call(rbind, all_runs)

# Plotting (Patch as factor ensures all appear in legend/plot)
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

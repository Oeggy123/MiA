library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# ----------------------------------------
# Seasonal Forcing Functions
# ----------------------------------------
alpha_t <- function(t, alpha_bar = 0.5, A = 0.3, T = 365, phi = 0) {
  alpha_bar * (1 + A * cos((2 * pi * t / T) + phi))
}

bv_t <- function(t, bv_bar = 0.25, A = 0.3, T = 365, phi = 0) {
  bv_bar * (1 + A * cos((2 * pi * t / T) + phi))
}

# ----------------------------------------
# Seasonal Logistic Metapopulation Model with Death Tracking
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
    Dh <- state[patch_indices("Dh")]
    
    Nh <- Sh + Eh + Ih + Rh
    Nv <- Sv + Ev + Iv
    
    alpha <- alpha_t(t, alpha_bar, A, T, phi)
    b_v <- bv_t(t, b_v_bar, A_b, T, phi_b)
    
    lambda_h <- f_h * alpha * p_h * Iv / (Nh + 1e-6)
    lambda_v <- f_h * alpha * p_v * Ih / (Nh + 1e-6)
    
    logistic_birth <- b_v * Nv * (1 - Nv / K_v)
    
    dSh <- b_h * Nh - lambda_h * Sh + omega * Rh - mu_h * Sh + inflow(Sh, M) - outflow(Sh, M)
    dEh <- lambda_h * Sh - sigma_h * Eh - mu_h * Eh + inflow(Eh, M) - outflow(Eh, M)
    dIh <- sigma_h * Eh - gamma * Ih - (mu_h + d) * Ih + inflow(Ih, M) - outflow(Ih, M)
    dRh <- gamma * Ih - omega * Rh - mu_h * Rh + inflow(Rh, M) - outflow(Rh, M)
    dDh <- d * Ih
    
    dSv <- logistic_birth - lambda_v * Sv - mu_v * Sv + inflow(Sv, Mv) - outflow(Sv, Mv)
    dEv <- lambda_v * Sv - sigma_v * Ev - mu_v * Ev + inflow(Ev, Mv) - outflow(Ev, Mv)
    dIv <- sigma_v * Ev - mu_v * Iv + inflow(Iv, Mv) - outflow(Iv, Mv)
    
    return(list(c(dSh, dEh, dIh, dRh, dSv, dEv, dIv, dDh)))
  })
}

# Indexing helpers
patch_indices <- function(var) {
  start <- match(var, c("Sh", "Eh", "Ih", "Rh", "Sv", "Ev", "Iv", "Dh"))
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

# Parameters with seasonal forcing
base_parameters <- list(
  b_h = 1.1e-5, b_v_bar = 0.25,
  f_h = 0.6, alpha_bar = 0.8,
  p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009,
  mu_h = 1.1e-5, mu_v = 0.1,
  d = 0.003,
  K_v = 5000,
  A = 0.3, A_b = 0.2,
  T = 365, phi = 0, phi_b = 0,
  np = np, patches = patches,
  M = M, Mv = Mv
)

# Initial conditions
init <- rep(0, 8 * np)
names(init) <- c(paste0("Sh", 1:np), paste0("Eh", 1:np), paste0("Ih", 1:np), paste0("Rh", 1:np),
                 paste0("Sv", 1:np), paste0("Ev", 1:np), paste0("Iv", 1:np), paste0("Dh", 1:np))

init["Sh1"] <- 1000; init["Eh1"] <- 20; init["Ih1"] <- 30; init["Sv1"] <- 6000; init["Ev1"] <- 150; init["Iv1"] <- 100
init["Sh2"] <- 600;  init["Sv2"] <- 4000; init["Ev2"] <- 80;  init["Iv2"] <- 40
init["Sh3"] <- 300;  init["Sv3"] <- 2000; init["Ev3"] <- 50;  init["Iv3"] <- 25

# ----------------------------------------
# Run Simulations for Multiple Horizons
# ----------------------------------------
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_outputs <- list()
all_stats <- list()

for (label in names(time_horizons)) {
  times <- seq(0, time_horizons[[label]], 1)
  out <- ode(y = init, times = times, func = metapop_model, parms = base_parameters, method = "ode45")
  out_df <- as.data.frame(out)
  out_df$time_horizon <- label
  all_outputs[[label]] <- out_df
  
  # Stats per patch
  stats <- data.frame(Patch = character(), Deaths = numeric(), Peak_Infection = numeric(), Peak_Day = numeric())
  for (i in 1:np) {
    deaths <- tail(out_df[[paste0("Dh", i)]], 1)
    peak_val <- max(out_df[[paste0("Ih", i)]])
    peak_day <- out_df$time[which.max(out_df[[paste0("Ih", i)]])]
    stats <- rbind(stats, data.frame(
      Patch = paste0("Patch ", i),
      Deaths = round(deaths),
      Peak_Infection = round(peak_val),
      Peak_Day = peak_day
    ))
  }
  stats$Horizon <- label
  all_stats[[label]] <- stats
}

# ----------------------------------------
# Combine for Plotting
# ----------------------------------------
combined_df <- bind_rows(all_outputs)
long_df <- combined_df %>%
  select(time, time_horizon, starts_with("Ih")) %>%
  pivot_longer(cols = starts_with("Ih"), names_to = "Patch", values_to = "Infected") %>%
  mutate(Patch = factor(Patch, levels = paste0("Ih", 1:np)))

# ----------------------------------------
# Plot Infected Over Time
# ----------------------------------------
for (h in names(time_horizons)) {
  df <- filter(long_df, time_horizon == h)
  p <- ggplot(df, aes(x = time, y = Infected, color = Patch)) +
    geom_line(size = 1) +
    labs(title = paste("Infected Over Time -", h),
         x = "Days", y = "Infected") +
    theme_minimal()
  print(p)
}

# ----------------------------------------
# Final Summary Table
# ----------------------------------------
final_stats <- bind_rows(all_stats)
print(final_stats)

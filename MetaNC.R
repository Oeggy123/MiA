library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Metapopulation SEIR-SEI Model with Death Tracking ---
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

patch_indices <- function(var) {
  start <- match(var, c("Sh", "Eh", "Ih", "Rh", "Sv", "Ev", "Iv", "Dh"))
  return(((start - 1) * np + 1):(start * np))
}
inflow <- function(vec, M) as.vector(t(M) %*% vec)
outflow <- function(vec, M) as.vector(M %*% vec)

# --- Shared Settings ---
np <- 3
patches <- 1:np
M <- matrix(c(
  0.9, 0.05, 0.05,
  0.05, 0.9, 0.05,
  0.05, 0.05, 0.90
), nrow = np, byrow = TRUE)
Mv <- M

base_parameters <- list(
  b_h = 1.1e-5, b_v = 0.25,
  f_h = 0.6, alpha = 0.8,
  p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009,
  mu_h = 1.1e-5, mu_v = 0.1,
  d = 0.003,
  K_v = 5000,
  np = np, patches = patches,
  M = M, Mv = Mv
)

init <- rep(0, 8 * np)
names(init) <- c(paste0("Sh", 1:np), paste0("Eh", 1:np), paste0("Ih", 1:np),
                 paste0("Rh", 1:np), paste0("Sv", 1:np), paste0("Ev", 1:np),
                 paste0("Iv", 1:np), paste0("Dh", 1:np))

init["Sh1"] <- 1000; init["Eh1"] <- 20; init["Ih1"] <- 30; init["Sv1"] <- 6000; init["Ev1"] <- 150; init["Iv1"] <- 100
init["Sh2"] <- 600;  init["Sv2"] <- 4000; init["Ev2"] <- 80;  init["Iv2"] <- 40
init["Sh3"] <- 300;  init["Sv3"] <- 2000; init["Ev3"] <- 50;  init["Iv3"] <- 25

# --- Loop through durations ---
horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
results_list <- list()
stats_list <- list()

for (label in names(horizons)) {
  times <- seq(0, horizons[[label]], 1)
  out <- ode(y = init, times = times, func = metapop_model, parms = base_parameters, method = "ode45")
  df <- as.data.frame(out)
  df$time_horizon <- label
  results_list[[label]] <- df
  
  # Summary stats for each patch
  stats <- data.frame(Patch = character(), Deaths = numeric(), Peak_Infection_Day = numeric(), First_Death_Day = numeric())
  for (i in 1:np) {
    deaths <- tail(df[[paste0("Dh", i)]], 1)
    peak_day <- df$time[which.max(df[[paste0("Ih", i)]])]
    first_death_day <- df$time[which(df[[paste0("Dh", i)]] >= 1)[1]]
    stats <- rbind(stats, data.frame(
      Patch = paste0("Patch ", i),
      Deaths = round(deaths),
      Peak_Infection_Day = peak_day,
      First_Death_Day = first_death_day
    ))
  }
  stats$Time_Horizon <- label
  stats_list[[label]] <- stats
}

# --- Combine for plotting ---
all_results <- bind_rows(results_list)
long_df <- all_results %>%
  select(time, time_horizon, starts_with("Ih")) %>%
  pivot_longer(cols = starts_with("Ih"), names_to = "Patch", values_to = "Infected") %>%
  mutate(Patch = factor(Patch, levels = paste0("Ih", 1:np)))

# --- Plot ---
for (h in names(horizons)) {
  df <- filter(long_df, time_horizon == h)
  p <- ggplot(df, aes(x = time, y = Infected, color = Patch)) +
    geom_line(size = 1.2) +
    labs(title = paste("Infectious Humans Over Time -", h),
         x = "Days", y = "Infected") +
    theme_minimal()
  print(p)
}

# --- Summary Table ---
final_summary <- bind_rows(stats_list)
print(final_summary)

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# Modified model with seasonal forcing, insecticide, logistic growth, and intervention scenarios
diff_HV_insecticide_seasonal <- function(t, pop, para){
  # Seasonal effects
  alpha_t <- para$alpha_bar * (1 + para$A * cos(2 * pi * t / para$T + para$phi))
  b_v_t <- para$b_v_bar * (1 + para$A_b * cos(2 * pi * t / para$T + para$phi_b))
  
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Vh <- pop[5]
  Ddh <- pop[6]; Dmh <- pop[7]; Sv <- pop[8]; Ev <- pop[9]; Iv <- pop[10]; Dv <- pop[11]
  
  Nh <- Sh + Eh + Ih + Rh + Vh
  Nv <- Sv + Ev + Iv
  
  # Intervention start
  delta_t <- ifelse(t >= para$t_insecticide, para$delta_v, 0)
  
  lambda_h <- para$f_h * alpha_t * para$p_h * Iv / Nh
  lambda_v <- para$f_h * alpha_t * para$p_v * Ih / Nh
  
  dSh <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dVh <- 0  # Placeholder
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Ih + Eh + Sh + Rh)
  
  logistic_birth <- b_v_t * Nv * (1 - Nv / para$K_v)
  dSv <- logistic_birth - lambda_v * Sv - (para$mu_v + delta_t) * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev - (para$mu_v + delta_t) * Ev
  dIv <- para$sigma_v * Ev - (para$mu_v + delta_t) * Iv
  dDv <- (para$mu_v + delta_t) * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dVh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

ODE_model_seasonal_insecticide <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_insecticide_seasonal, parms = para, method = "ode45")
  df <- as.data.frame(result)
  colnames(df)[1] <- "time"
  return(df)
}

# Parameters
base_para <- list(
  b_h = 1.1e-5, b_v_bar = 0.25, f_h = 0.6,
  alpha_bar = 0.8, p_h = 0.2, p_v = 0.5,
  sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  nu = 0, epsilon_v = 0, omega_v = 0,
  A = 0.3, T = 365, phi = 0, A_b = 0.2, phi_b = 0,
  delta_v = 0.05,
  t_insecticide = Inf,
  K_v = 5000
)

ICs <- c(Sh = 500, Eh = 10, Ih = 30, Rh = 0, Vh = 0,
         Ddh = 0, Dmh = 0, Sv = 4000, Ev = 100, Iv = 50, Dv = 0)

# Time horizons
time_horizons <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)
all_data <- list()

for (h in names(time_horizons)) {
  maxtime <- time_horizons[[h]]
  
  # Baseline to find first death
  base_para$t_insecticide <- Inf
  baseline <- ODE_model_seasonal_insecticide(base_para, ICs, maxtime)
  first_death_day <- which(baseline$Ddh >= 1)[1] - 1
  if (is.na(first_death_day)) next
  
  timings <- list(
    "Early (D - 30)" = max(0, first_death_day - 30),
    "On First Death (D)" = first_death_day,
    "Delayed (D + 30)" = first_death_day + 30,
    "Very Delayed (D + 90)" = first_death_day + 90
  )
  
  runs <- lapply(names(timings), function(scen) {
    p <- base_para
    p$t_insecticide <- timings[[scen]]
    out <- ODE_model_seasonal_insecticide(p, ICs, maxtime)
    out$scenario <- scen
    out$time_horizon <- h
    return(out)
  })
  
  all_data[[h]] <- bind_rows(runs)
}

final_data <- bind_rows(all_data)

# Split plots
for (h in unique(final_data$time_horizon)) {
  cat("Plotting for horizon:", h, "\n")
  df <- final_data %>% filter(time_horizon == h)
  
  # Infections
  print(ggplot(df, aes(x = time, y = Ih, color = scenario)) +
          geom_line() + labs(title = paste("Infections (", h, ")"), x = "Day", y = "Ih") + theme_minimal())
  
  # Deaths
  print(ggplot(df, aes(x = time, y = Ddh, color = scenario)) +
          geom_line() + labs(title = paste("Cumulative Deaths (", h, ")"), x = "Day", y = "Deaths") + theme_minimal())
  
  # SEIR
  seir_long <- df %>%
    select(time, scenario, Sh, Eh, Ih, Rh) %>%
    pivot_longer(cols = c("Sh", "Eh", "Ih", "Rh"), names_to = "Compartment", values_to = "Count")
  
  print(ggplot(seir_long, aes(x = time, y = Count, color = Compartment)) +
          geom_line() + facet_wrap(~ scenario) +
          labs(title = paste("SEIR Dynamics (", h, ")"), x = "Day", y = "Count") +
          theme_minimal() + theme(legend.position = "bottom"))
}

summary_stats <- final_data %>%
  group_by(time_horizon, scenario) %>%
  summarise(
    final_Ddh = round(max(Ddh), 2),
    peak_Ih = round(max(Ih), 2),
    peak_day = time[which.max(Ih)],
    .groups = "drop"
  )

print(summary_stats)

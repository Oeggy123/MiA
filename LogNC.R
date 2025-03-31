##############################################################
# Logistic Malaria SEIR-SEI Model: 1yr, 5yr, 10yr runs
##############################################################

library(deSolve)
library(pracma)

# -----------------------------
# Logistic SEIR-SEI Model
# -----------------------------
diff_HV_logistic <- function(t, pop, para){
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Ddh <- pop[5]; Dmh <- pop[6]
  Sv <- pop[7]; Ev <- pop[8]; Iv <- pop[9]; Dv <- pop[10]
  
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  lambda_h <- para$f_h * para$alpha * para$p_h * Iv / Nh
  lambda_v <- para$f_h * para$alpha * para$p_v * Ih / Nh
  
  # Human dynamics
  dSh <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh)
  
  # Logistic mosquito birth
  logistic_birth <- para$b_v * Nv * (1 - Nv / para$K_v)
  
  dSv <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv <- para$sigma_v * Ev - para$mu_v * Iv
  dDv <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# ODE Solver
ODE_model_logistic <- function(para, ICs, maxtime){
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_logistic, parms = para, method = "ode45")
  return(as.data.frame(result))
}

# -----------------------------
# Parameter Setup
# -----------------------------
para <- list(
  b_h = 1.1e-5, b_v = 0.25, f_h = 0.6, alpha = 0.8,
  p_h = 0.2, p_v = 0.5, sigma_h = 1/10, sigma_v = 1/15,
  gamma = 1/20, omega = 0.0009, d = 0.003,
  mu_h = 1.1e-5, mu_v = 0.1,
  K_v = 5000
)

# Initial Conditions
ICs <- c(
  Sh = 500, Eh = 10, Ih = 30, Rh = 0,
  Ddh = 0, Dmh = 0,
  Sv = 4000, Ev = 100, Iv = 50, Dv = 0
)

# Timeframes to run
time_spans <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)

# -----------------------------
# Loop through durations
# -----------------------------
for (label in names(time_spans)) {
  
  maxtime <- time_spans[[label]]
  cat("\n------ Running Logistic Model:", label, "------\n")
  
  # Run model
  Classes <- ODE_model_logistic(para, ICs, maxtime)
  
  # Plot
  pdf(paste0("RFigs/logistic_nocontrol_", label, ".pdf"), width = 7, height = 6)
  plot(Classes$time, Classes$Ih, type = "l", col = "red", xlab = "Time (days)", ylab = "Population", 
       ylim = c(0, 500), xlim = c(0, maxtime), lwd = 2)
  lines(Classes$time, Classes$Sh, col = "black", lwd = 2)
  lines(Classes$time, Classes$Eh, col = "darkorange", lwd = 2)
  lines(Classes$time, Classes$Rh, col = "darkblue", lwd = 2)
  
  legend("topright", legend = c("S", "E", "I", "R"), col = c("black", "darkorange", "red", "darkblue"), lty = 1, lwd = 2)
  dev.off()
  
  # Output stats
  total_deaths <- round(tail(Classes$Ddh, 1))
  peak_index <- which.max(Classes$Ih)
  peak_time <- Classes$time[peak_index]
  peak_value <- round(Classes$Ih[peak_index])
  
  cat(paste("Total deaths after", label, ":", total_deaths, "\n"))
  cat(paste("Peak infection on day:", peak_time, "\n"))
  cat(paste("Peak number of infected individuals:", peak_value, "\n"))
}

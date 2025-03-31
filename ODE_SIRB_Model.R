#########################################################
# Malaria SEIR-SEI Model with Bed Nets Only (No Vaccination)
#########################################################

library(deSolve)
library(pracma)

# Differential equations function
diff_HV_bednets <- function(t, pop, para){
  
  # Extract population compartments
  Sh <- pop[1]
  Eh <- pop[2]
  Ih <- pop[3]
  Rh <- pop[4]
  Ddh <- pop[5]
  Dmh <- pop[6]
  Sv <- pop[7]
  Ev <- pop[8]
  Iv <- pop[9]
  Dv <- pop[10]
  
  # Total populations
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  # Time-dependent bed net efficacy (optional)
  eta_t <- ifelse(t >= para$t_bednets, para$eta, 0)
  
  # Adjusted biting rate
  alpha_eff <- (1 - eta_t) * para$alpha
  
  # Forces of infection
  lambda_h <- para$f_h * alpha_eff * para$p_h * Iv / Nh
  lambda_v <- para$f_h * alpha_eff * para$p_v * Ih / Nh
  
  # Human equations
  dSh  <- para$b_h * Nh - lambda_h * Sh + para$omega * Rh - para$mu_h * Sh
  dEh  <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh
  dIh  <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh  <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh)
  
  # Vector equations
  dSv  <- para$b_v * Nv - lambda_v * Sv - para$mu_v * Sv
  dEv  <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv  <- para$sigma_v * Ev - para$mu_v * Iv
  dDv  <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# ODE wrapper function
ODE_SEIRSEI_bednets <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_bednets, parms = para, method = "ode45")
  return(as.data.frame(result))
}

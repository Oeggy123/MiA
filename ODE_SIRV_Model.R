diff_HV_det_logistic <- function(t, pop, para){
  
  # Extract population compartments
  Sh <- pop[1]
  Eh <- pop[2]
  Ih <- pop[3]
  Rh <- pop[4]
  Vh <- pop[5]
  Ddh <- pop[6]
  Dmh <- pop[7]
  Sv <- pop[8]
  Ev <- pop[9]
  Iv <- pop[10]
  Dv <- pop[11]
  
  # Total populations
  Nh <- Sh + Eh + Ih + Rh + Vh
  Nv <- Sv + Ev + Iv
  
  # Time-dependent vaccination of susceptibles
  nu_t <- ifelse(t >= para$t_vax, para$nu_rate, 0)
  
  # Force of infection
  lambda_h <- para$f_h * para$alpha * para$p_h * Iv / Nh
  lambda_h_v <- (1 - para$epsilon_v) * lambda_h
  lambda_v <- para$f_h * para$alpha * para$p_v * Ih / Nh
  
  # Human dynamics
  dSh  <- (1 - para$nu) * para$b_h * Nh - lambda_h * Sh - nu_t * Sh + para$omega * Rh + para$omega_v * Vh - para$mu_h * Sh
  dEh  <- lambda_h * Sh - para$sigma_h * Eh - para$mu_h * Eh + lambda_h_v * Vh
  dIh  <- para$sigma_h * Eh - para$gamma * Ih - (para$mu_h + para$d) * Ih
  dRh  <- para$gamma * Ih - para$omega * Rh - para$mu_h * Rh
  dVh  <- para$nu * para$b_h * Nh + nu_t * Sh - para$omega_v * Vh - lambda_h_v * Vh - para$mu_h * Vh
  dDdh <- para$d * Ih
  dDmh <- para$mu_h * (Sh + Eh + Ih + Rh + Vh)
  
  # Logistic mosquito birth
  logistic_birth <- para$b_v * Nv * (1 - Nv / para$K_v)
  
  # Vector dynamics with logistic growth
  dSv  <- logistic_birth - lambda_v * Sv - para$mu_v * Sv
  dEv  <- lambda_v * Sv - para$sigma_v * Ev - para$mu_v * Ev
  dIv  <- para$sigma_v * Ev - para$mu_v * Iv
  dDv  <- para$mu_v * (Sv + Ev + Iv)
  
  return(list(c(dSh, dEh, dIh, dRh, dVh, dDdh, dDmh, dSv, dEv, dIv, dDv)))
}

# Wrapper to solve the ODE system with logistic vectors
ODE_SEIRSEI_logistic <- function(para, ICs, maxtime) {
  t_seq <- seq(0, maxtime, by = 1)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_det_logistic, parms = para, method = "ode45")
  return(as.data.frame(result))
}
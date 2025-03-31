#########################################################

#MiA

#########################################################

############################################################################
#Define a function that sets out the differential equations to be solved
############################################################################

diff_HV_det <- function(t, pop, para){
  
  Sh <- pop[1]; Eh <- pop[2]; Ih <- pop[3]; Rh <- pop[4]; Ddh <- pop[5]
  Sv <- pop[6]; Ev <- pop[7]; Iv <- pop[8]
  
  Nh <- Sh + Eh + Ih + Rh
  Nv <- Sv + Ev + Iv
  
  lambda_h <- para$f_h * para$alpha * para$p_h * Iv / Nh
  lambda_v <- para$f_h * para$alpha * para$p_v * Ih / Nh
  
  dSh <- -lambda_h * Sh + para$omega * Rh
  dEh <- lambda_h * Sh - para$sigma_h * Eh
  dIh <- para$sigma_h * Eh - para$gamma * Ih - para$d * Ih
  dRh <- para$gamma * Ih - para$omega * Rh
  dDdh <- para$d * Ih
  
  dSv <- -lambda_v * Sv
  dEv <- lambda_v * Sv - para$sigma_v * Ev
  dIv <- para$sigma_v * Ev
  
  return(list(c(dSh, dEh, dIh, dRh, dDdh, dSv, dEv, dIv)))
}

##############################################################################
#Define the ODE_SIR_inf function 
##############################################################################

ODE_SEIRSEI_model <- function(para, ICs, maxtime) {
  
  #Time points for simulation
  t_seq <- seq(0, maxtime, by = 1)
  
  #Solve the ODEs using the ode function from deSolve package (ode 45 suitable for our well-behaved system)
  result <- ode(y = ICs, times = t_seq, func = diff_HV_det, parms = para, method = "ode45")
  
  #Convert the result to a data frame
  Classes <- data.frame(result)
  
  return(Classes)
}

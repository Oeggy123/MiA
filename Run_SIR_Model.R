##############################################################

#MiA

##############################################################

library(deSolve) #Access deSolve (for ode solving) 
library(pracma) #Access pracma (for trapezoidal integration)
source("ODE_SIR_Model.R")

##############################################################
#Model setup and run
##############################################################

#Define model parameters
para <- list(b_h = 1.1e-5, b_v = 0.25, f_h = 0.6, alpha = 0.8,
             p_h = 0.2, p_v = 0.5, sigma_h = 1/10, sigma_v = 1/15,
             gamma = 1/20, omega = 0.0009, d = 0.003,
             mu_h = 1.1e-5, mu_v = 0.1)

#Define initial conditions
ICs <- c("Sh" = 500 , "Eh" = 10, "Ih" = 30, "Rh" = 0, "Ddh" = 0, "Sv" = 4000, "Ev" = 100, "Iv" = 50)

#Define time to run model for
maxtime <- 2000

#Run model 
Classes <- ODE_SEIRSEI_model(para, ICs, maxtime)

##############################################################
#Plotting I and total deaths
##############################################################

#produce an SIR plot with different colours
pdf("RFigs/nocontrol.pdf", width = 7, height = 6)
plot(Classes$t, Classes$Ih, type = "l", col = "red", xlab = "Time (days)", ylab = "Population", 
     ylim = c(0, 500), xlim = c(0, maxtime), lwd = 2)
lines(Classes$t, Classes$Sh, type = "l", col = "black", lty = 1, lwd = 2)
lines(Classes$t, Classes$Eh, type = "l", col = "darkorange", lty = 1, lwd = 2)
lines(Classes$t, Classes$Rh, type = "l", col = "darkblue", lty = 1, lwd = 2)

#legend for plot
legend("topright", legend = c("S", "E", "I", "R"), col = c("black", "darkorange", "red","darkblue"), lty = 1, lwd = 2)

dev.off()

# Calculate total deaths using trapezoidal integration
total_deaths <- round(tail(Classes$Ddh, 1))


# Print total deaths
print(paste("Total deaths:", total_deaths))

# Find peak infection time and value
peak_index <- which.max(Classes$Ih)
peak_time <- Classes$t[peak_index]
peak_value <- Classes$Ih[peak_index]

# Print results
print(paste("Peak infection occurs at day:", peak_time))
print(paste("Peak number of infected individuals:", round(peak_value)))






# Simulation durations (in days)
time_spans <- list("1yr" = 365, "5yr" = 1825, "10yr" = 3650)

# Loop over durations
for (label in names(time_spans)) {
  
  maxtime <- time_spans[[label]]
  cat("\n----- Running", label, "simulation -----\n")
  
  # Run model
  Classes <- ODE_SEIRSEI_model(para, ICs, maxtime)
  
  # Plot
  pdf(paste0("RFigs/nocontrol_", label, ".pdf"), width = 7, height = 6)
  plot(Classes$t, Classes$Ih, type = "l", col = "red", xlab = "Time (days)", ylab = "Population",
       ylim = c(0, 500), xlim = c(0, maxtime), lwd = 2)
  lines(Classes$t, Classes$Sh, col = "black", lwd = 2)
  lines(Classes$t, Classes$Eh, col = "darkorange", lwd = 2)
  lines(Classes$t, Classes$Rh, col = "darkblue", lwd = 2)
  legend("topright", legend = c("S", "E", "I", "R"),
         col = c("black", "darkorange", "red", "darkblue"), lty = 1, lwd = 2)
  dev.off()
  
  # Output stats
  total_deaths <- round(tail(Classes$Ddh, 1))
  peak_index <- which.max(Classes$Ih)
  peak_time <- Classes$t[peak_index]
  peak_value <- Classes$Ih[peak_index]
  
  cat(paste("Total deaths after", label, ":", total_deaths, "\n"))
  cat(paste("Peak infection occurs at day:", peak_time, "\n"))
  cat(paste("Peak number of infected individuals:", round(peak_value), "\n"))
}


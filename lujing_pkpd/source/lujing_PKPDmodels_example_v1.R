# Header----
require(data.table)
require(PKPDmodels)
require(ggplot2)

# Look-up a PKPD model
# Source: 
PKPDmodels::PKexpr(admin = "bolus",
                   dosage = "sd",
                   cpt = 2)

# Initial conditions/true parameter values----
# Assume initial dose was 
init <- list(A = 10,
             ka = 0.5,
             B = 5,
             kb = 1)

# Simulate data----
dt1 <- data.table(Time = rep(1:10, 
                             3))

dt1$Conc <- jitter(init$A*exp(-init$ka*dt1$Time) + init$B*exp(-init$kb*dt1$Time),
                   amount = 0.5)

# Estimate the parameters with Nonlinear Least Squares (nls)----
m1 <- nls(Conc ~ A*exp(-ka*Time) + B*exp(-kb*Time),
          data = dt1,
          start = init,
          weights = rep(1,
                        nrow(dt1)),
          control = nls.control(maxiter = 1000, 
                                minFactor = 1/4096),
          na.action = na.omit)
s1 <- summary(m1)
s1
s1$parameters

# Predicted data----
tmp <- data.table(Predicted = predict(m1,
                                      newdata = list(Time = unique(dt1$Time))),
                  Time = unique(dt1$Time))
tmp

# Plot 
plot(dt1$Conc ~ dt1$Time)
points(tmp$Predicted  ~ tmp$Time,
       pch = 21,
       bg = "red")
lines(tmp$Predicted  ~ tmp$Time)

# Area under predicted concentraion curve
# i.e. total amount of drug passed through the system
tmp

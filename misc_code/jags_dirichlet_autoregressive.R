# Header ------------------------------------------------------------------

# Dirichlet Autoregressive model of order 1
# Andrew Parnell

# Some JAGS code to fit a Dirichlet AR(1) model

# Some boiler plate code to clear the workspace, and load in required packages
rm(list=ls()) # Clear the workspace
library(R2jags)
library(bayesm)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y(t) = compositional vector of repsonses at time t, t = 1,...,T
# y_r(t) = compositional value for species r, r = 1, ..., R
# alpha_r = overall mean parameter for species r
# beta_r = autocorrelation/autoregressive (AR) parameter for species r
# a_r(t) = dirichlet parameter for species r at time t

# Likelihood
# y[t, 1:R] ~ ddirch(a[t, 1:R])

# Second level log(a[t, r]) = alpha_r + beta_r * log(a[t-1, r])

# Priors
# alpha_r ~ dnorm(0,100)
# beta_r ~ dunif(-1,1) # If you want the process to be stable/stationary

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
# First AR1
set.seed(123)
T = 100
t_seq = 1:T
R = 4
alpha = rep(1, R)
beta = runif(R, 0.2, 0.8)
y = a =  matrix(NA, nrow = T, ncol = 4)
a[1,] = exp(alpha)
y[1,] = rdirichlet(a[1,])
for(t in 2:T) {
  for(r in 1:R) a[t,r] = exp( alpha[r] + beta[r] * log(a[t-1, r]))
  y[t,] = rdirichlet(a[t,])
}
# plot
plot(t_seq,y[,1],type='l', ylim = c(0, 1))
lines(t_seq,y[,2])
lines(t_seq,y[,3])
lines(t_seq,y[,4])

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
# This code is for a general AR(p) model

model_code = '
model
{
  for (r in 1:R) {
    log(a[1, r]) <- alpha[r]
  }

  # Likelihood
  for (t in 2:T) {
    y[t, 1:R] ~ ddirch(a[t, 1:R])
    for (r in 1:R) {
      a[t, r] <- exp(log_a[t, r])
      log_a[t, r] <- alpha[r] + beta[r] * y[t-1, r]#log_a[t-1, r]
    }
  }

  # Priors
  for (r in 1:R) {
    alpha[r] ~ dnorm(0, 10^-2)
    beta[r] ~ dnorm(0, 10^-2)
  }
}
'

# Set up the data
model_data = list(T = T, R = R, y = y)

# Choose the parameters to watch
model_parameters =  c("alpha","beta","a")

init_fun = function() {
  list('alpha' = rep(1, R), 'beta' = rep(0, R))
}

# Run the model
model_run = jags(data = model_data,
                 #inits = init_fun,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code))

# Simulated results -------------------------------------------------------

# Check the output - are the true values inside the 95% CI?
# Also look at the R-hat values - they need to be close to 1 if convergence has been achieved
print(model_run)
print(model_run_2) # Note: phi is correct but in the wrong order

# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory
hadcrut = read.csv('https://raw.githubusercontent.com/andrewcparnell/tsme_course/master/data/hadcrut.csv')
head(hadcrut)
with(hadcrut,plot(Year,Anomaly,type='l'))

# Look at the ACF/PACF
acf(hadcrut$Anomaly)
pacf(hadcrut$Anomaly)

# Set up the data
real_data = with(hadcrut,
                 list(T = nrow(hadcrut),
                      y = hadcrut$Anomaly,
                      p = 1))

# Run the model
real_data_run = jags(data = real_data,
                     parameters.to.save = model_parameters,
                     model.file=textConnection(model_code),
                     n.chains=4,
                     n.iter=1000,
                     n.burnin=200,
                     n.thin=2)

# Plot output
print(real_data_run) # Very high degree of autocorrelation

# Plot some of the fitted values (also known as one-step-ahead predictions)
post = print(real_data_run)
alpha_mean = post$mean$alpha
phi_mean = post$mean$phi

# Create fitted values
fitted_values = alpha_mean + phi_mean * real_data$y[1:(nrow(hadcrut)-1)]

# Create fitted line
with(hadcrut, plot(Year, Anomaly, type='l'))
with(hadcrut, lines(Year[2:nrow(hadcrut)], fitted_values, col='red'))
# Why does this look strange?

# Create some predictions off into the future
T_future = 2050
future_values = rep(NA, T_future-max(hadcrut$Year))
future_values[1] = alpha_mean + phi_mean * real_data$y[nrow(hadcrut)]
for (i in 2:length(future_values)) {
  future_values[i] = alpha_mean + phi_mean * future_values[i-1]
}

# Plot these all together
with(hadcrut,
     plot(Year,
          Anomaly,
          type='l',
          xlim=c(min(hadcrut$Year),T_future),
          ylim=range(c(hadcrut$Anomaly,future_values))))
lines(((max(hadcrut$Year)+1):T_future),future_values,col='red')
# See - no global warming!

# Other tasks -------------------------------------------------------------

# 1) Try changing the values of phi2 for the simulated AR(p) model. What happens to the time series when some of these values get bigger?
# 2) Above we have only fitted the HadCrut data with an AR(1) model. You might like to try and fit it with AR(2), AR(3), etc models and see what happens to the fits
# 3) (Harder) See if you can create the fitted values by sampling from the posterior distribution of alpha and phi, and plotting an envelope/ensemble of lines, just like in the linear regression example

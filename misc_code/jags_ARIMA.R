# Header ------------------------------------------------------------------

# ARIMA (AutoRegressive Integrated Moving Average) fit in JAGS
# Andrew Parnell

# This model is a combination of the two previous files called jags_autoregressive and jags_moving_average. The differencing part (the I in integrated) takes place outside of JAGS, so we're really just fitting an ARMA model

# Some boiler plate code to clear the workspace, set the working directory, and load in required packages
rm(list=ls()) # Clear the workspace
setwd("~/GitHub/tsme_course/")
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# This is for a general ARIMA(p,d,q) model
# Notation:
# y(t) = response variable at time t, t=1,...,T
# alpha = mean parameter
# eps_t = residual at time t
# theta = MA parameters
# phi = AR parameters
# sigma = residual standard deviation
# d = number of first differences
# p and q = number of autoregressive and moving average components respecrively
# Likelihood:
# y ~ N(alpha + phi[1] * y[t-1] + ... + phi[p] * y[y-p] + theta_1 ept_{t-1} + ... + theta_q eps_{t-q}, sigma^2)
# Priors
# alpha ~ N(0,100)
# phi ~ N(0,100) - need to be a bit careful with these if you want the process to remain stable
# theta ~ N(0,100)
# sigma ~ unif(0,10)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
p = 1 # Number of autoregressive terms
d = 0 # Number of differences
q = 1 # Numner of MA terms
T = 100
sigma = 1
alpha = 0
set.seed(123)
theta = runif(q)
phi = sort(runif(p),decreasing=TRUE)
y = rep(NA,T)
y[1:q] = rnorm(q,0,sigma)
eps = rep(NA,T)
eps[1:q] = y[1:q] - alpha
for(t in (q+1):T) {
  ar_mean = sum( phi * y[(t-1):(t-p)] )
  ma_mean = sum( theta * eps[(t-q):(t-1)] )
  y[t] = rnorm(1, mean = alpha + ar_mean + ma_mean, sd = sigma)
  eps[t] = y[t] - alpha - ma_mean - ar_mean
}
plot(1:T,y,type='l')

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Set up residuals
  for(t in 1:q) {
    eps[t] <- y[t] - alpha
  }
  # Likelihood
  for (t in (q+1):T) {
    y[t] ~ dnorm(alpha + ar_mean[t] + ma_mean[t], tau)
    ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
    ar_mean[t] <- inprod(phi, y[(t-p):(t-1)])
    eps[t] <- y[t] - alpha - ar_mean[t] - ma_mean[t]
  }
  
  # Priors
  alpha ~ dnorm(0.0,0.01)
  for (i in 1:q) {
    theta[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:p) {
    phi[i] ~ dnorm(0.0,0.01)
  }
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
}
'

# Set up the data
model_data = list(T = T, y = y, q = 1, p = 1)

# Choose the parameters to watch
model_parameters =  c("alpha","theta","phi","sigma")

# Run the model
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run) # Parameters theta/phi/sigma should match the true value

# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory
hadcrut = read.csv('data/hadcrut.csv')
head(hadcrut)
with(hadcrut,plot(Year,Anomaly,type='l'))

# Look at the ACF/PACF
acf(hadcrut$Anomaly)
pacf(hadcrut$Anomaly)

# Set up the data
real_data = with(hadcrut,
                 list(T = nrow(hadcrut),
                      y = hadcrut$Anomaly,
                      q = 1,
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
theta_mean = post$mean$theta
phi_mean = post$mean$phi

# Create fitted values
eps_fit = y_fit = rep(NA,real_data$T)
eps_fit[1:real_data$q] = y[1:real_data$q] - alpha_mean
y_fit[1:real_data$q] = alpha_mean
for (t in (real_data$q+1):real_data$T) {
  ar_mean = sum( phi_mean * real_data$y[(t-1):(t-real_data$p)] )
  ma_mean = sum( theta_mean * eps_fit[(t-real_data$q):(t-1)] )
  eps_fit[t] = real_data$y[t] - alpha_mean - ar_mean - ma_mean
  y_fit[t] = alpha_mean + ar_mean + ma_mean
}

# Create fitted line
with(hadcrut, plot(Year, Anomaly, type='l'))
with(hadcrut, lines(Year, y_fit, col='red'))
# Not a bad fit!

# Create some predictions off into the future - this time do it within jags
# A neat trick - just increase T and add on NAs into y!
T_future = 20 # Number of future data points
real_data_future = with(hadcrut,
                        list(T = nrow(hadcrut) + T_future,
                             y = c(hadcrut$Anomaly, rep(NA,T_future)),
                             q = 1,
                             p = 1))

# Just watch y now
model_parameters =  c("y")

# Run the model
real_data_run_future = jags(data = real_data_future,
                     parameters.to.save = model_parameters,
                     model.file=textConnection(model_code),
                     n.chains=4,
                     n.iter=1000,
                     n.burnin=200,
                     n.thin=2)

# Print out the above
print(real_data_run_future)

# Get the future values
y_all = real_data_run_future$BUGSoutput$sims.list$y
# If you look at the above object you'll see that the first columns are all identical because they're the data
y_all_mean = apply(y_all,2,'mean')
year_all = c(hadcrut$Year,(max(hadcrut$Year)+1):(max(hadcrut$Year)+T_future))

# Plot these all together
plot(year_all,
     y_all_mean,
     type='n')
lines(year_all,y_all_mean,col='red')
with(hadcrut,lines(Year,Anomaly))

# Other tasks -------------------------------------------------------------

# 1) Calculate the estimated residuals from the AIMRA(1,0,1) fit and plot them. Just like a linear regression they should form a random pattern and be approx normally distributed. Check whether this is so using hist, qqplot, and acf
# 2) If you look at the fits for the ARIMA(1,0,1) model you'll see that phi is approximately 1. Try creating a model instead which is ARIMA(0,1,1) - this is a moving average model on the differences (go back to the jags_moving_average file to fit this mdoel). Do the fits look any differences
# 3) (harder) Try experimenting with the order of the random walk and look at the effect on the fits. Try changing the length of future predictions or adding in confidence intervals (Hint: see the code in the jags_linear_regression file for how to produce confidence intervals from the output; it involves changing the all to apply in the lines above)



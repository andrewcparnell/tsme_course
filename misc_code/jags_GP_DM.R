# Header ------------------------------------------------------------------

# Gaussian process models in JAGS
# Andrew Parnell

# This file fits a Gaussian Process (GP) regression model to data in JAGS, and produces predictions/forecasts

# Some boiler plate code to clear the workspace and load in required packages
rm(list=ls())
library(R2jags)
library(MASS) # Useful for mvrnorm function

# Maths -------------------------------------------------------------------

# Notation:
# y(t): Response variable at time t, defined on continuous time
# y: vector of all observations
# alpha: Overall mean parameter
# sigma: residual standard deviation parameter (sometimes known in the GP world as the nugget)
# rho: decay parameter for the GP autocorrelation function
# tau: GP standard deviation parameter

# Likelihood:
# y ~ MVN(Mu, Sigma)
# where MVN is the multivariate normal distribution and
# Mu[t] = alpha
# Sigma is a covariance matrix with:
# Sigma_{ij} = tau^2 * exp( -rho * (t_i - t_j)^2 ) if i != j
# Sigma_{ij} = tau^2 + sigma^2 if i=j
# The part exp( -rho * (t_i - t_j)^2 ) is known as the autocorrelation function

# Prior
# alpha ~ N(0,100)
# sigma ~ U(0,10)
# tau ~ U(0,10)
# rho ~ U(0.1, 5) # Need something reasonably informative here

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
T = 10 # can take to T = 100 (default) but fitting gets really slow ...
alpha = 0 # default is 0
sigma = 0.01 # default is 0.01
tau = 1
rho = 1
set.seed(123)
t = sort(runif(T))
Sigma = sigma^2 * diag(T) + tau^2 * exp( - rho * outer(t,t,'-')^2 )
y = mvrnorm(1,rep(alpha,T), Sigma)
plot(t,y)

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  y ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:T) {
    Mu[i] <- alpha
    Sigma[i,i] <- pow(sigma, 2) + pow(tau, 2)
  
    for(j in (i+1):T) {
      Sigma[i,j] <- pow(tau, 2) * exp( - rho * pow(t[i] - t[j], 2) )
      Sigma[j,i] <- Sigma[i,j]
    }
  }
  
  alpha ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10) # default dunif(0,10)
  tau ~ dunif(0, 10)
  rho ~ dunif(0.1, 5)
  
} 
'
  
# Set up the data
model_data = list(T = T, y = y, t = t)
  
# Choose the parameters to watch
model_parameters =  c("alpha", "sigma", "tau", "rho")
  
# Run the model - can be slow
model_run = jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file=textConnection(model_code),
                   n.chains=4, # Number of different starting positions
                   n.iter=1000, # Number of iterations
                   n.burnin=200, # Number of iterations to remove at start
                   n.thin=2) # Amount of thinning
  

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run)

# Now create some predictions of new values at new times t^new
# These are bsed on the formula:
# y^new | y ~ N( Mu^new + Sigma_new^T solve(Sigma, y - Mu), Sigma_* - Sigma_new^T solve(Sigma, Sigma_new)
# where 
# Mu^new[t] = alpha  
# Sigma_new[i,j] = tau^2 * exp( -rho * (t^new_i - t_j)^2 )
# Sigma_*[i,j] = tau^2 * exp( -rho * (t^new_i - t^new_j)^2 ) if i != j

# First look at parameters
alpha = model_run$BUGSoutput$sims.list$alpha
tau = model_run$BUGSoutput$sims.list$tau
sigma = model_run$BUGSoutput$sims.list$sigma
rho = model_run$BUGSoutput$sims.list$rho
par(mfrow = c(2,2))
hist(alpha, breaks=30)
hist(tau, breaks=30)
hist(sigma, breaks=30)
hist(rho, breaks=30)

par(mfrow=c(1,1))

# Now create predictions
T_new = 20
t_new = seq(0,1,length=T_new)
Mu = rep(mean(alpha), T)
Mu_new = rep(mean(alpha), T_new)
Sigma_new = mean(tau)^2 * exp( -mean(rho) * outer(t, t_new, '-')^2 )
Sigma_star = mean(sigma)^2 * diag(T_new) + mean(tau)^2 * exp( - mean(rho) * outer(t_new,t_new,'-')^2 )
Sigma = mean(sigma)^2 * diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
pred_mean = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)

# Plot output
plot(t,y)
points(t_new, pred_mean, col='red')
lines(t_new, pred_mean, col='red')

pred_low = pred_mean - 1.95 * sqrt(diag(pred_var))
pred_high = pred_mean + 1.95 * sqrt(diag(pred_var))
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)

# Take samples from the posterior distribution to show the possible solutions.
for(i in 1500:1600){
  Mu = rep(alpha[i], T)
  Mu_new = rep(alpha[i], T_new)
  Sigma_new = tau[i]^2 * exp( -rho[i] * outer(t, t_new, '-')^2 )
  Sigma_star = sigma[i]^2 * diag(T_new) + tau[i]^2 * exp( - rho[i] * outer(t_new,t_new,'-')^2 )
  Sigma = sigma[i]^2 * diag(T) + tau[i]^2 * exp( - rho[i] * outer(t,t,'-')^2 )
  
  # Use fancy equation to get predictions
  pred_samp = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
  #pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)
  lines(t_new, pred_samp, col='grey')
}



# Fit a a sine wave with a trend and noise

T = 20
t = 1:T
y <- rnorm(T,t/5 + sin(seq(from = 0, to = 10, length.out=T)), 0.3)
plot(t,y)

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  y ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:T) {
  Mu[i] <- alpha
  Sigma[i,i] <- pow(sigma, 2) + pow(tau, 2)
  
  for(j in (i+1):T) {
  Sigma[i,j] <- pow(tau, 2) * exp( - rho * pow(t[i] - t[j], 2) )
  Sigma[j,i] <- Sigma[i,j]
  }
  }
  
  alpha ~ dnorm(0, 0.01)
  sigma ~ dunif(0, 10) # default dunif(0,10)
  tau ~ dunif(0, 1)
  rho ~ dunif(0.1, 5)
  
} 
  '
  
# Set up the data
model_data = list(T = T, y = y, t = t)
  
# Choose the parameters to watch
model_parameters =  c("alpha", "sigma", "tau", "rho")
  
# Run the model - can be slow
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run)

alpha = model_run$BUGSoutput$sims.list$alpha
tau = model_run$BUGSoutput$sims.list$tau
sigma = model_run$BUGSoutput$sims.list$sigma
rho = model_run$BUGSoutput$sims.list$rho
par(mfrow = c(2,2))
hist(alpha, breaks=30)
hist(tau, breaks=30)
hist(sigma, breaks=30)
hist(rho, breaks=30)
par(mfrow=c(1,1))

# Now create predictions
T_new = 30
t_new = seq(0,20,length=T_new)
Mu = rep(mean(alpha), T)
Mu_new = rep(mean(alpha), T_new)
Sigma_new = mean(tau)^2 * exp( -mean(rho) * outer(t, t_new, '-')^2 )
Sigma_star = mean(sigma)^2 * diag(T_new) + mean(tau)^2 * exp( - mean(rho) * outer(t_new,t_new,'-')^2 )
Sigma = mean(sigma)^2 * diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
pred_mean = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)

# Plot output
plot(t,y)
lines(t_new, pred_mean, col='red')

pred_low = pred_mean - 1.95 * sqrt(diag(pred_var))
pred_high = pred_mean + 1.95 * sqrt(diag(pred_var))
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)


# doesn't work yet ...

T = 20
t = 1:T
y <- rnorm(T,sin(seq(from = 0, to = 10, length.out=T)), 0.1)
plot(t,y)

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  y ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:T) {
  Mu[i] <- alpha
  Sigma[i,i] <- pow(tau, 2)
  
  for(j in (i+1):T) {
  Sigma[i,j] <- pow(tau, 2) * exp( - rho * pow(t[i] - t[j], 2) )
  Sigma[j,i] <- Sigma[i,j]
  }
  }
  
  alpha ~ dnorm(0, 0.01)
  tau ~ dunif(0, 10)
  rho ~ dunif(0.1, 5)
  
} 
  '
  
# Set up the data
model_data = list(T = T, y = y, t = t)

# Choose the parameters to watch
model_parameters =  c("alpha", "tau", "rho")

# Run the model - can be slow
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

alpha = model_run$BUGSoutput$sims.list$alpha
tau = model_run$BUGSoutput$sims.list$tau
rho = model_run$BUGSoutput$sims.list$rho
par(mfrow = c(2,2))
hist(alpha, breaks=30)
hist(tau, breaks=30)
hist(sigma, breaks=30)
hist(rho, breaks=30)

par(mfrow=c(1,1))

# Now create predictions
T_new = 30
t_new = seq(0,20,length=T_new)
Mu = rep(mean(alpha), T)
Mu_new = rep(mean(alpha), T_new)
Sigma_new = mean(tau)^2 * exp( -mean(rho) * outer(t, t_new, '-')^2 )
Sigma_star = diag(T_new) + mean(tau)^2 * exp( - mean(rho) * outer(t_new,t_new,'-')^2 )
Sigma = diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
pred_mean = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)

# Plot output
plot(t,y)
points(t_new, pred_mean, col='red')
lines(t_new, pred_mean, col='red')

pred_low = pred_mean - 1.95 * sqrt(diag(pred_var))
pred_high = pred_mean + 1.95 * sqrt(diag(pred_var))
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)


# Real example ------------------------------------------------------------


model_code = '
model
{
  # Likelihood
  y ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)
  
  # Set up mean and covariance matrix
  for(i in 1:T) {
  Mu[i] <- alpha
  Sigma[i,i] <- pow(sigma, 2) + pow(tau, 2)
  
  for(j in (i+1):T) {
  Sigma[i,j] <- pow(tau, 2) * exp( - rho * pow(t[i] - t[j], 2) )
  Sigma[j,i] <- Sigma[i,j]
  }
  }
  
  alpha ~ dnorm(919, 100) # default dnorm(0,0.01)
  sigma ~ dunif(0, 300) # default dunif(0,10)
  tau ~ dunif(0, 10)
  rho ~ dunif(0.1, 5)
  
} 
  '


# Data wrangling and jags code to run the model on a real data set in the data directory
t = seq(from = start(Nile)[1], to = end(Nile)[1])
y = as.vector(Nile)
T = length(Nile)

alpha = 919 # default is 0
sigma = 0.01 # default is 0.01
tau = 1
rho = 1

real_data = list(T = T,
                 y = y,
                 t = t
                 )

# Choose the parameters to watch
real_parameters =  c("alpha", "sigma", "tau", "rho")

# Run the model - can be slow
real_run = jags(data = real_data,
                 parameters.to.save = real_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

alpha = real_run$BUGSoutput$sims.list$alpha
tau = real_run$BUGSoutput$sims.list$tau
sigma = real_run$BUGSoutput$sims.list$sigma
rho = real_run$BUGSoutput$sims.list$rho

# Now create predictions
T_new = 50
t_new = seq(1871,1970,length=T_new)
Mu = rep(mean(alpha), T)
Mu_new = rep(mean(alpha), T_new)
Sigma_new = mean(tau)^2 * exp( -mean(rho) * outer(t, t_new, '-')^2 )
Sigma_star = mean(sigma)^2 * diag(T_new) + mean(tau)^2 * exp( - mean(rho) * outer(t_new,t_new,'-')^2 )
Sigma = mean(sigma)^2 * diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
pred_mean = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)

# Plot output
plot(t,y)
points(t_new, pred_mean,type = 'l', col='red')
lines(t_new, pred_mean, col='red')

pred_low = pred_mean - 1.95 * sqrt(diag(pred_var))
pred_high = pred_mean + 1.95 * sqrt(diag(pred_var))
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)

# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks



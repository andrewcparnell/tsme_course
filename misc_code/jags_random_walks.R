# Header ------------------------------------------------------------------

# Random walk models (first and second difference)
# Andrew Parnell

# In this code we fit some random walk type models to data

# Some boiler plate code to clear the workspace, set the working directory, and load in required packages
rm(list=ls()) # Clear the workspace
setwd("~/GitHub/tsme_course/")
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y(t) = response variable at time t, t = 1,...,T
# mu: Optional drift parameter
# sigma = residual standard deviation

# Likelihood:
# Order 1: y(t) - y(t-1) ~ N(mu,sigma^2)
# Order 2: y(t) - 2y(t-1) + y(t-2) ~ N(mu,sigma^2)
# Prior:
# sigma ~ unif(0,100) - vague
# sigma ~ dgamma(a,b) ~ informative with good values for a and b
# mu ~ dnorm(0,100) - vague

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
set.seed(123)
T = 100
sigma = 1
mu = 0
t = 1:T
y = cumsum(rnorm(T,mu,sigma))
y2 = cumsum(cumsum(rnorm(T,mu,sigma)))
plot(t,y)
plot(t,y2)

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
# Note: running hte differencing offline here as part of the data step
model_code = '
model
{
  # Likelihood
  for (t in 2:N_T) {
    z[t] ~ dnorm(z[t-1] + mu, tau)
  }
  
  # Priors
  mu ~ dnorm(0.0,0.01) 
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,100.0) 
}
'

# Set up the data
order = 1
model_data = list(z = diff(y, differences=order), N_T = T-order)

# Choose the parameters to watch
model_parameters =  c("mu","sigma")  

# Run the model
model_run = jags(data = model_data, 
                 parameters.to.save = model_parameters, 
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning
stop()

# Try the order 2 version
order = 2
model_data_2 = list(z = diff(y2, differences=order), N_T = T-order)

model_run_2 = jags(data = model_data_2, 
                   parameters.to.save = model_parameters, 
                   model.file=textConnection(model_code),
                   n.chains=4, # Number of different starting positions
                   n.iter=1000, # Number of iterations
                   n.burnin=200, # Number of iterations to remove at start
                   n.thin=2) # Amount of thinning

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run)
print(model_run_2) # Slight over-estimation of sigma here

# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory
hadcrut = read.csv('data/hadcrut.csv')
head(hadcrut)
with(hadcrut,plot(Year,Anomaly,type='l'))

# Set up the data
order = 1
real_data = with(hadcrut,
                 list(T = nrow(hadcrut), 
                      z = hadcrut$Anomaly, 
                      p = 1))

# Run the model
real_data_run = jags(data = real_data, 
                     parameters.to.save = model_parameters, 
                     model.file=textConnection(model_code),
                     n.chains=4, 
                     n.iter=1000,
                     n.burnin=200,
                     n.thin=2)



# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks



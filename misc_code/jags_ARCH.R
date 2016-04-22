# Header ------------------------------------------------------------------

# Autoregressive conditional heteroskesticity (ARCH) models
# Andrew Parnell

# An ARCH model is just like an AR model but with the AR component applied to the variance instead. This script just contains an ARCH(1) model

# Some boiler plate code to clear the workspace, set the working directory, and load in required packages
rm(list=ls()) # Clear the workspace
setwd("~/GitHub/tsme_course/")
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_t = response variable at time t=1,...,T
# alpha = overall mean
# sigma_t = residual standard deviation at time t
# gamma_1 = mean of variance term
# gamma_2 = AR component of variance
# Likelihood:
# y_t ~ N( alpha, sigma_t^2)
# sigma_t^2 = gamma_1 + gamma_2 * sigma_{t-1}^2
# Priors
# gamma_1 ~ unif(0,10) - needs to be positive
# gamma_2 ~ unif(0,1) - ditto, and usually <1 too
# alpha ~ N(0,100) - vague

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
T = 100
alpha = 1
gamma_1 = 1
gamma_2 = 0.4
set.seed(123)
sigma = rep(0,length=T)
for(t in 2:T) sigma[t] = sqrt(gamma_1 + gamma_2 * sigma[t-1]^2)
y = rnorm(T,mean=alpha,sd=sigma)
plot(1:T,y,type='l')

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  for (t in 1:T) {
    y[t] ~ dnorm(alpha, tau[t])
    tau[t] <- 1/pow(sigma[t], 2)
  }
  sigma[1] <- 1
  for(t in 2:T) {  
    sigma[t] <- sqrt( gamma_1 + gamma_2 * pow(sigma[t-1], 2) )
  }
  
  # Priors
  alpha ~ dnorm(0.0, 0.01)
  gamma_1 ~ dunif(0, 10)
  gamma_2 ~ dunif(0, 1)
}
'

# Set up the data
model_data = list(T = T, y = y)

# Choose the parameters to watch
model_parameters =  c("gamma_1","gamma_2","alpha")

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
plot(model_run)
print(model_run)

# Real example ------------------------------------------------------------

# Run the ARCH(1) model on the ice core data set
ice = read.csv('data/GISP2_20yr.csv')
head(ice)
with(ice,plot(Age,Del.18O,type='l'))
# Try plots of differences
with(ice,plot(Age[-1],diff(Del.18O,differences=1),type='l'))
with(ice,plot(Age[-(1:2)],diff(Del.18O,differences=2),type='l'))

# Try this on the last 30k years
ice2 = subset(ice,Age<=30000)
table(diff(ice2$Age))
with(ice2,plot(Age[-1],diff(Del.18O),type='l'))

# Set up the data
real_data = with(ice2,
                 list(T = nrow(ice2) - 1, y = diff(Del.18O)))

# Save the sigma's the most interesting part!
model_parameters = c('sigma')

# Run the model - requires longer to converge
real_data_run = jags(data = real_data,
                     parameters.to.save = model_parameters,
                     model.file=textConnection(model_code),
                     n.chains=4,
                     n.iter=10000,
                     n.burnin=2000,
                     n.thin=8)

print(real_data_run)

# Plot the sigma outputs
sigma_med = apply(real_data_run$BUGSoutput$sims.list$sigma,2,'quantile',0.5)
sigma_low = apply(real_data_run$BUGSoutput$sims.list$sigma,2,'quantile',0.025)
sigma_high = apply(real_data_run$BUGSoutput$sims.list$sigma,2,'quantile',0.975)

plot(ice2$Age[-(1:10)],sigma_med[-c(1:9)],type='l')


# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks
# 1) Try playing with the values of gamma_1 and gamma_2 in the simulated data above. See if you can create some really crazy patterns (e.g. try gamma_2>1)
# 2) 
# 3) (harder) The above model is only an ARCH(1) model. See if you can simulate from and then fit an ARCH(2) version.

# Header ------------------------------------------------------------------

# Change point modelling in JAGS for Anna
# Andrew Parnell

# This file fits a random effects/longitudinal change point model based on SNP
# frequencies 

# Some boiler plate code to clear the workspace, and load in required packages
rm(list=ls())
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y(j,t) = frequency observed at generation t (t=1, .., 7) for replicate j (j=1, .., 10). This is going to be a percentage between 0 and 100, sometimes hitting 100
# alpha_j = mean frequency value at the change point
# beta_j = slope of the frequency before the change point
# t_j = change point generation value 

# Likelihood:
# Top level likelihood is always:
# y(j,t) ~ normal(mu[j,t], sigma^2)

# Change point pattern
# mu[j,t] = alpha_j + beta_j * (t-t_j) if t<t_j, or alpha_j if t>= t_j

# To achieve this kind of model in jags we use the step function which works via:
# step(x) = 1 if x>0 or 0 otherwise.
# We can use it to pick out which side of the change point(s) we're on

# Priors
# alpha_j ~ normal(0, 100)
# beta_j ~ normal(0, 100)
# sigma ~ uniform(0, 100)
# t_j ~ uniform(0, 7) # Restrict the change point to the range of the data

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model

# DCPR-1 model
T = 7 # Number of generations
J = 10 # Number of replicates
sigma = 1
set.seed(123)
alpha = runif(J, 40, 97)
beta = runif(J, 5, 10)
t_j = runif(J, 0, 6) # Time of change points for each replicate
t = 1:T
mu = y = matrix(NA, ncol = J, nrow = T) 
for(j in 1:J) {
  mu[t<t_j[j],j] = alpha[j] + beta[j] * (t[t<t_j[j]] - t_j[j])
  mu[t>=t_j[j],j] = alpha[j]
  y[,j] = rnorm(T, mu[,j], sigma)
}
plot(t,y[,1], type='n', ylim = c(0, 100))
for(j in 1:J) {
  lines(t, y[,j], col=j)  
}

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data

# Code
model_code_anna = "
model
{
  # Likelihood
  for(i in 1:T) {
    for(j in 1:J) {
      y[i,j] ~ dnorm(mu[i,j], tau)
      mu[i,j] <- alpha[j] + beta[K[i,j]]*(t[i]-cp[j])
      K[i,j] <- 1 + step(t[i] - cp[j])
    }
  }
  
  # Priors
  for(j in 1:J) {
    alpha[j] ~ dnorm(alpha_0, tau_alpha)
    beta[j] ~ dnorm(beta_0, tau_beta)
    cp[j] ~ dunif(t_min, t_max)
  }
  alpha_0 ~ dnorm(0, 0.001)
  beta_0 ~ dnorm(0, 0.001)
  
  tau <- 1/pow(sigma, 2)
  tau_alpha <- 1/pow(sigma_alpha, 2)
  tau_beta <- 1/pow(sigma_beta, 2)
  sigma ~ dunif(0, 100)
  sigma_alpha ~ dunif(0, 100)
  sigma_beta ~ dunif(0, 100)
  }
"

# Choose the parameters to watch
model_parameters_anna =  c("cp", "alpha_0", "alpha")

# Create data
model_data_anna = list(T = T,
                       J = J,
                       y = y,
                       t = 1:T,
                       t_min = 0.9,
                       t_max = 6.9)

# Run the model
model_run_anna = jags(data = model_data_anna,
                        parameters.to.save = model_parameters_anna,
                        model.file=textConnection(model_code_anna),
                        n.chains=4,
                        n.iter=1000,
                        n.burnin=200,
                        n.thin=2)

print(model_run_anna)

cp_mean = model_run_anna$BUGSoutput$mean$cp
abline(v=cp_mean, col=1:10, lty='dotted')
hist(model_run_anna$BUGSoutput$sims.list$alpha_0, breaks=30)


# Header ------------------------------------------------------------------

# Multivariate Autoregressive models, commonly known as Vector AR (VAR) models
# Andrew Parnell

# The VAR model is a multivariate extension of the standard AR(p) model. In this code I just fit the VAR(1) model but it is easily extended to VAR(p)

# Some boiler plate code to clear the workspace and load in required packages
rm(list=ls())
library(R2jags)
library(MASS) # Used to generate MVN samples
library(MCMCpack)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation
# y_t = multivariate response variable at time t, t=1,...,T. Each y_t is a vetor of dimension k
# A = k-vector of constants
# Phi = k by k matrix of AR coefficients - these are our most important parameters
# e_t = k-vector of residuals
# Sigma = k by k matrix of residual variance and co-variances

# Likelihood
# y_t = A + Phi * y_{t-1} + e_t with e_t ~ MVN(0, Sigma)
# or
# y_t ~ MVN(A + Phi * y_{t-1}, Sigma)

# Prior
# A[k] ~ normal(0, 100)
# Phi[j,k] ~ normal(0, 100)
# Sigma ~ Inverse Wishart(I, k+1)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
T = 100
k = 2
Sigma = matrix(c(1, 0.2, 0.2, 1), 2, 2)
Phi = matrix(c(0.6, 0.2, 0.2, 0.8), 2, 2)
A = matrix(c(0, 2), 2, 1)
y = matrix(NA, T, k)
y[1,] = A
set.seed(123)
for(t in 2:T) {
  y[t,] = mvrnorm(1, A + Phi %*% y[t-1,], Sigma)
}

# Plot the output
par(mfrow = c(2, 1))
plot(1:T, y[,1], type = 'l')
plot(1:T, y[,2], type = 'l')
par(mfrow = c(1, 1))

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Likelihood
  for (t in 2:T) {
    y[t, ] ~ dmnorm(mu[t, ], Sigma.Inv)
    mu[t, 1:k] <- A + Phi %*% y[t-1,]
  }

  # Priors
  for(i in 1:k) {
    A[i] ~ dnorm(0, 0.01)
    for(j in 1:k) {
      Phi[i,j] ~ dnorm(0, 0.01)
    }
  }
  Sigma.Inv ~ dwish(I, k+1)
}
'

# Set up the data
model_data = list(T = T, k = k, y = y, I = diag(k))

# Choose the parameters to watch
model_parameters =  c("A", "Phi", "Sigma.Inv")

# Run the model
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=10000, # Number of iterations
                 n.burnin=2000, # Number of iterations to remove at start
                 n.thin=8) # Amount of thinning

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run) # Results look pretty good

# Real example ------------------------------------------------------------

# Can we fit a vector AR model to both the sea level and global temperature
# series?




# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks



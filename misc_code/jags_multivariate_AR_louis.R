# Header ------------------------------------------------------------------

# Multivariate Autoregressive model with explanatory variables, now known as Vector AR with eXplanatory variables (VARX) model
# Andrew Parnell

# Some boiler plate code to clear the workspace and load in required packages
rm(list=ls())
library(R2jags)
library(MASS) # Used to generate MVN samples

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation
# y_t = multivariate response variable at time t, t=1,...,T. Each y_t is a vector of dimension k
# A = k-vector of constants
# Phi = k by k matrix of AR coefficients
# Beta = k-vector of regression coefficients
# e_t = k-vector of residuals
# Sigma = k by k matrix of residual variance and co-variances

# Likelihood
# y_t = A + Phi * y_{t-1} + Beta * x_t + e_t with e_t ~ MVN(0, Sigma)
# or
# y_t ~ MVN(A + Phi * y_{t-1} + Beta * x_t , Sigma)

# Prior
# A[k] ~ normal(0, 100)
# Beta[k] ~ normal(0, 100)
# Phi[j,k] ~ normal(0, 100)
# Sigma ~ Inverse Wishart(I, k+1)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
T = 40
k = 2
Sigma = matrix(c(1, 0.2, 0.2, 1), 2, 2)
Phi = matrix(c(0.6, 0.2, 0.2, 0.8), 2, 2)
A = matrix(c(0, 2), 2, 1)
Beta = matrix(c(2, 0, 0, 3), 2, 2)
y = matrix(NA, T, k)
y[1,] = A
set.seed(123)
x = matrix(runif(T*k), ncol=k, nrow=T)
for(t in 2:T) {
  y[t,] = mvrnorm(1, A + Phi %*% y[t-1,] + Beta %*% x[t,], Sigma)
}

# Plot the output
par(mfrow = c(2, 1))
plot(1:T, y[,1], type = 'l')
plot(1:T, y[,2], type = 'l')
par(mfrow = c(1, 1))

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code_louis = '
model
{
  # Likelihood
  for (t in 2:T) {
    y[t, ] ~ dmnorm(mu[t, ], Sigma.Inv)
    mu[t, 1:k] <- A + Phi %*% y[t-1,] + Beta %*% x[t,]
  }
  Sigma.Inv ~ dwish(I, k+1)
  Sigma <- inverse(Sigma.Inv)  

  # Priors
  for(i in 1:k) {
    A[i] ~ dnorm(0, 0.01)
    Phi[i,i] ~ dunif(-1, 1)
    Beta[i,i] ~ dnorm(0, 0.01)
    for(j in (i+1):k) {
      Phi[i,j] ~ dunif(-1,1)
      Phi[j,i] ~ dunif(-1,1)
      Beta[i,j] <- 0
      Beta[j,i] <- 0
    }
  }
}
'

# Set up the data
model_data_louis = list(T = T, k = k, y = y, x = x, I = diag(k))

# Choose the parameters to watch
model_parameters_louis =  c("A", "Phi", "Sigma", "Beta")

# Run the model
model_run_louis = jags(data = model_data_louis,
                 parameters.to.save = model_parameters_louis,
                 model.file=textConnection(model_code_louis),
                 n.chains=4, # Number of different starting positions
                 n.iter=10000, # Number of iterations
                 n.burnin=2000, # Number of iterations to remove at start
                 n.thin=8) # Amount of thinning

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
print(model_run_louis) # Results look pretty good


# Future predictions ------------------------------------------------------

T_future = 10
x_future = matrix(runif((T_future)*k), ncol=k, nrow=T_future)
x_all = rbind(x, x_future)

# Set up the data
model_data_louis_f = list(T = T + T_future, k = k, y = rbind(y, matrix(NA, ncol=k, nrow=T_future)), x = x_all, I = diag(k))

# Choose the parameters to watch
model_parameters_louis_f =  c("y")

# Run the model
model_run_louis_f = jags(data = model_data_louis_f,
                       parameters.to.save = model_parameters_louis_f,
                       model.file=textConnection(model_code_louis),
                       n.chains=4, # Number of different starting positions
                       n.iter=10000, # Number of iterations
                       n.burnin=2000, # Number of iterations to remove at start
                       n.thin=8) # Amount of thinning

# Get future predictions and plot
y_future = model_run_louis_f$BUGSoutput$mean$y

T_all = 1:(T+T_future)

# Plot the output
par(mfrow = c(2, 1))
plot(T_all, y_future[,1], type = 'l')
plot(T_all, y_future[,2], type = 'l')
par(mfrow = c(1, 1))

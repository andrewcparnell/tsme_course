# Header ------------------------------------------------------------------

# Change point modelling for Igor in JAGS
# Andrew Parnell

# This file implements a single continuous change point model in JAGS with covariates. The slope before the change point should be positive, and the slope after should be negative

# Some boiler plate code to clear the workspace, and load in required packages
rm(list=ls())
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y(i,t) = response variable for observation i at times t. there may be multiple observations at each time
# mu[t] = true mean value of y taking account of multiple observations
# cp = time of change point 
# gamma[1,k] is the effects ofexplanatory variable k x_k(t) before the change point
# gamma[2,k] is the effects of explanatory variable k x_k(t) after the change point
# alpha = intercept term
# beta_k = slope value for period k
# sigma_y = standard deviation between observations at the same time point
# sigma_mu = overall residual standard deviation

# Likelihood:
# Top level likelihood is always:
# y(i,t) ~ normal(mu[t], sigma_y^2)
# mu[t] ~ normal(b[t], sigma_mu^2) 
# with b[t] = alpha + beta[1] * (t - cp) + gamma[1,] %*% x[t,] if t < cp, or mu[t] = alpha + beta[2] * (t - cp)  + gamma[2,] * x[t,] if t>=cp
# Note that this is a clever way of expresssing the model as it means that alpha is the mean of y at the change point

# Priors
# alpha ~ normal(0, 100)
# gamma[k] ~ normal(0, 100)
# beta[1] ~ uniform(0, 100)
# beta[2] ~ uniform(-100, 0)
# sigma_y ~ uniform(0, 100)
# sigma_mu ~ uniform(0, 100)
# cp ~ uniform(t_min, t_max) # Restrict the change point to the range of the data

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model

# CCPR-1 model
T = 100
set.seed(123)
num_obs_per_time = sample(1:5, T, replace=TRUE) # Num obs per time point
t_seq = 1:T
t_seq_2 = rep(t_seq, times = num_obs_per_time)
sigma_y = 5
sigma_b = 2
alpha = 50
beta = c(0.5,-0.8) # Slopes before and after change point
# Explanatory variables
K = 3
x = matrix(runif(T*K), ncol=K, nrow=T)
gamma = matrix(rnorm(K*2, 0, 10), ncol=2, nrow=K)
cp = runif(1, 0, T) # Time of change point
b = rep(NA,T)
b[t_seq<cp] = alpha + beta[1] * (t_seq[t_seq<cp] - cp) + x[t_seq<cp,]%*%gamma[,1,drop=FALSE]
b[t_seq>=cp] = alpha + beta[2] * (t_seq[t_seq>=cp] - cp) + x[t_seq>=cp,]%*%gamma[,2,drop=FALSE]
b = b + rnorm(T, 0, sigma_b)
y = rnorm(length(t_seq_2), b[t_seq_2], sigma_y)
plot(t_seq_2,y)
lines(t_seq[t_seq<cp], b[t_seq<cp] -  x[t_seq<cp,]%*%gamma[,1,drop=FALSE], col='red')
lines(t_seq[t_seq>=cp], b[t_seq>=cp]  - x[t_seq>=cp,]%*%gamma[,2,drop=FALSE], col='red')
lines(t_seq[t_seq<cp], b[t_seq<cp], col='blue')
lines(t_seq[t_seq>=cp], b[t_seq>=cp], col='blue')
abline(v=cp,col='green')

# Store this all together for later use
model_data_igor = list(t=t_seq, 
                       y=y, 
                       x=x, 
                       T=T, 
                       K=K,
                       T_big = length(t_seq_2), 
                       b_select = t_seq_2, 
                       t_min = min(t_seq), 
                       t_max = max(t_seq))

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data

# Code for CCPR-1
model_code_igor="
model
{
  # Likelihood
  for(i in 1:T_big) {
    y[i] ~ dnorm(mu[i], tau_y)
    mu[i] <- b[b_select[i]]
  }
  for(j in 1:T) {
    b[j] ~ dnorm(alpha + beta[J[j]]*(t[j]-cp) + x[j,]%*%gamma[,J[j]], tau_b)
    # This is the clever bit - only pick out the right change point when above cp
    J[j] <- 1 + step(t[j] - cp)
  }

  # Priors
  alpha ~ dnorm(0.0, 0.01)
  for(k in 1:K) {
    for(i in 1:2) {
      gamma[k,i] ~ dnorm(0.0, 0.01)
    }
  } 
  beta[1] ~ dunif(0, 100)
  beta[2] ~ dunif(-100, 0)
  cp ~ dunif(t_min, t_max)

  tau_y <- 1/pow(sigma_y, 2)
  sigma_y ~ dunif(0, 100)
  tau_b <- 1/pow(sigma_b, 2)
  sigma_b ~ dunif(0, 100)
}
"

# Choose the parameters to watch
model_parameters_igor =  c("cp", "alpha", "beta", "sigma_y", "sigma_b", "gamma")

# Run the model
model_run_igor = jags(data = model_data_igor,
                        parameters.to.save = model_parameters_igor,
                        model.file=textConnection(model_code_igor),
                        n.chains=4,
                        n.iter=1000,
                        n.burnin=200,
                        n.thin=2)

print(model_run_igor)

# Key things are: slope before and after the change point
slopes = model_run_igor$BUGSoutput$sims.list$beta
par(mfrow=c(1,2))
hist(slopes[,1], main='Before cp', breaks=30)
hist(slopes[,2], main='After cp', breaks=30)
par(mfrow=c(1,1))
# The time of the change point
cp = model_run_igor$BUGSoutput$sims.list$cp
hist(cp, breaks=30, main='Time of change point')
# And the effect of the explanatory variables
explan = model_run_igor$BUGSoutput$sims.list$gamma
par(mfrow=c(K,2))
for(j in 1:K) {
    hist(explan[,j,1], main='Before cp', breaks=30)
    hist(explan[,j,2], main='After cp', breaks=30)
}
par(mfrow=c(1,1))

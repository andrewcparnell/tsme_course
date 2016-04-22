# Header ------------------------------------------------------------------

# A simple Bayesian Fourier model to produce a periodogram
# Andrew Parnell

# This model creates a periodogram of the data and applies to the Lynx data set example

# Some boiler plate code to clear the workspace, set the working directory, and load in required packages
rm(list=ls()) # Clear the workspace
#setwd("~/GitHub/tsme_course/")
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_t = Response variable at time t, t=1,...,T
# alpha = Overall mean parameters
# beta_k = cosine associated frequency coefficient k = 1,..., K frequencies in total
# gamma_k = sine associated frequency coefficient
# f_k = frequency value k
# sigma = residual standard deviation

# Likelihood:
# y_t ~ N( mu_t, sigma^2)
# mu_t = sum_k beta_k * cos ( 2 * pi * t * f_k) + gamma_k * sin ( 2 * pi * t * f_k )
# K and f_k are data and are set in advance

# Priors
# alpha ~ normal(0, 100)
# beta_k ~ normal(0, 100)
# gamma_k ~ normal(0, 100)
# sigma ~ uniform(0, 100)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
T = 100
K = 20
sigma = 1
alpha = 0
set.seed(123)
f = seq(0,2*pi,length=K)
beta = rnorm(K, 0, 1)
gamma = rnorm(K, 0, 1)
X = outer(2 * pi * 1:T, f, '*') # This creates a clever matrix of 2 * pi * t * f_k for every combination of t and f_k
mu = alpha + cos(X) %*% beta + sin(X) %*% gamma
y = rnorm(T, mu, sigma)
plot(1:T, y, type='l')

# Look at the acf/pacf
acf(y)
pacf(y)


# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data


# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc

# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory


# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks



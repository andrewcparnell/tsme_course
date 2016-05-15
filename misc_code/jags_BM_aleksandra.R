# Header ------------------------------------------------------------------

# Brownian motion in jags - version for Alexandra's data with different resolutions
# Andrew Parnell

# The general task here is to compare two time series, one with high resolution
# observations (the 'high-res' series) over a short time period, and the other
# with low resolution observation ('low-res' series) over a much longer time period

# This code proceeds by:
# 1. Assuming that both data sets share the same underlying parameters and fitting a model to both series. Here this is a Brownian Motion with drift
# 2. Using the NA trick to predict both series on a regular grid (here called t_grid)
# 3. Comparing the high-res series with the low-res series for the time grid
# The standard maths used can be found in jags_BM.R

# Some boiler plate code to clear the workspace, and load in required packages
rm(list=ls()) # Clear the workspace
library(R2jags)

# Simulate data -----------------------------------------------------------

# Some R code to simulate two data sets

# Set the seed to get reproducible data
set.seed(123)

# First simulate the high res data
# Assume 1000 years of data and 100 data points
T_hr = 100 # hr for high res
alpha = 0 # Drift parameter
sigma_bm = 1 # Variabiliy due to the random walk/brownian motion
sigma_hr = 0.3 # Extra variability at each observation
t_hr_low = 0 # Lower limit for high res series
t_hr_high = 1000 # Upper limit for high res series
t_hr = sort(runif(T_hr, t_hr_low, t_hr_high)) # Times for high res series
y_hr = rep(NA, T_hr)
y_hr[1] = rnorm(1, mean = alpha, sd = sigma_hr) # Simulate the first obs
for(i in 2:T_hr) y_hr[i] = y_hr[i-1] + rnorm(1, alpha * (t_hr[i] - t_hr[i-1]), sigma_bm * sqrt(t_hr[i] - t_hr[i-1])) + rnorm(1, 0, sigma_hr) # Simulate the later observations based on the random walk and the extra variability

# Now simulate the low-res data
# Now 10k years of data and 100 data points
T_lr = 100 # lr for low res
t_lr_low = 0 # Lower limit of low-res data
t_lr_high = 10000 # Upper time limit
sigma_lr = 2.1 # Extra variability at each observation
t_lr = sort(runif(T_lr, t_lr_low, t_lr_high)) # Simulate times
y_lr = rep(NA, T_lr) # Set up y obs
y_lr[1] = rnorm(1, mean = alpha, sd = sigma_lr) # Simulate first one
for(i in 2:T_lr) y_lr[i] = y_lr[i-1] + rnorm(1, alpha * (t_lr[i] - t_lr[i-1]), sigma_lr * sqrt(t_lr[i] - t_lr[i-1])) + rnorm(1, 0, sigma_lr) # Simulate remaining obs based on BM plus extra variation

# Plot them together
par(mfrow=c(2,1))
plot(t_hr, y_hr, type = 'l', xlim=range(c(t_lr, t_hr)))
plot(t_lr, y_lr, type = 'l', xlim=range(c(t_lr, t_hr)))
par(mfrow=c(1,1))

# Run on high-res data ----------------------------------------------------

# Jags code to fit the model to both data sets 
model_code = '
model
{
  # Likelihood - for high res data
  for(i in 1:T_hr_all) {
    y_hr_all[i] ~ dnorm(mean_hr[i], sigma_hr^-2)
  }
  # This is the random walk/brownian motion bit
  mean_hr[1] ~ dnorm(alpha, 0.01)
  for (i in 2:T_hr_all) {
    mean_hr[i] ~ dnorm( alpha * (t_hr_all[i] - t_hr_all[i-1]) + mean_hr[i-1], tau_hr_all[i] )
    tau_hr_all[i] <- 1/( pow(sigma_bm,2) * (t_hr_all[i] - t_hr_all[i-1]) )
  }

  # Likelihood for low res data
  for(i in 1:T_lr_all) {
    y_lr_all[i] ~ dnorm(mean_lr[i], sigma_lr^-2)
  }
  # This is the random walk/brownian motion bit
  mean_lr[1] ~ dnorm(alpha, 0.01)
  for (i in 2:T_lr_all) {
    mean_lr[i] ~ dnorm( alpha * (t_lr_all[i] - t_lr_all[i-1]) + mean_lr[i-1], tau_lr_all[i] )
    tau_lr_all[i] <- 1/( pow(sigma_bm,2) * (t_lr_all[i] - t_lr_all[i-1]) )
  }

  # Priors for parameters
  alpha ~ dnorm(0.0,0.01)
  sigma_bm ~ dunif(0.0,100.0)
  sigma_lr ~ dunif(0.0,100.0)
  sigma_hr ~ dunif(0.0,100.0)
}
'

# Create the overall grid
t_grid = seq(0, t_hr_high, by = 10)
y_hr_grid = y_lr_grid = rep(NA, length(t_grid)) # Use the NA trick

# Now put everything together 
t_hr_all_raw = c(t_hr, t_grid)
t_lr_all_raw = c(t_lr, t_grid)
y_hr_all_raw = c(y_hr, y_hr_grid)
y_lr_all_raw = c(y_lr, y_lr_grid)

# Use the order command to put the NAs inside the vectors in the correct time order
order_hr = order(t_hr_all_raw)
order_lr = order(t_lr_all_raw)

# Set up the data - use the order created previously 
model_data = list(T_hr_all = length(t_hr_all_raw),
                  T_lr_all = length(t_lr_all_raw),
                  t_hr_all = t_hr_all_raw[order_hr],
                  t_lr_all = t_lr_all_raw[order_lr],
                  y_hr_all = y_hr_all_raw[order_hr],
                  y_lr_all = y_lr_all_raw[order_lr])

# Choose the parameters to watch
model_parameters =  c("alpha",
                      "sigma_bm",  
                      "sigma_lr",
                      "sigma_hr",
                      "y_hr_all", 
                      "y_lr_all")

# Run the model
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

# Check output for convergence
print(model_run)
plot(model_run)

# Compare plots -----------------------------------------------------------

# Compare the predictions from the two versions on top of each other

# First get the output of the high res data
y_pred_hr_all = model_run$BUGSoutput$sims.list$y_hr_all
pick_hr_NA = which(is.na(y_hr_all_raw[order_hr]))

# Create 95% CIs and median
y_pred_hr_low = apply(y_pred_hr_all[,pick_hr_NA],2,'quantile',0.025)
y_pred_hr_med = apply(y_pred_hr_all[,pick_hr_NA],2,'quantile',0.5)
y_pred_hr_high = apply(y_pred_hr_all[,pick_hr_NA],2,'quantile',0.975)

# Do the same for low res
y_pred_lr_all = model_run$BUGSoutput$sims.list$y_lr_all
pick_lr_NA = which(is.na(y_lr_all_raw[order_lr]))

y_pred_lr_low = apply(y_pred_lr_all[,pick_lr_NA],2,'quantile',0.025)
y_pred_lr_med = apply(y_pred_lr_all[,pick_lr_NA],2,'quantile',0.5)
y_pred_lr_high = apply(y_pred_lr_all[,pick_lr_NA],2,'quantile',0.975)

# Now plot together
y_range = range(c(y_pred_hr_low, y_pred_hr_high, y_pred_lr_low, y_pred_lr_high))
plot(t_grid, y_pred_lr_med, type = 'n', ylim = y_range, ylab = 'y')
lines(t_grid, y_pred_lr_low, lty = 2)
lines(t_grid, y_pred_lr_med)
lines(t_grid, y_pred_lr_high, lty = 2)
lines(t_grid, y_pred_hr_low, lty = 2, col='red')
lines(t_grid, y_pred_hr_med, col='red')
lines(t_grid, y_pred_hr_high, lty = 2, col='red')

# Add a legend
legend('topleft', 
       legend = c('High-res series (with 95% CI)', 
                  "Low-res series (with 95% CI)"),
       lty=1, 
       bty = 'n',
       col=c('black', 'red'))



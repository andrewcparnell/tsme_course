# Header ------------------------------------------------------------------


# AR-RE model - AutoRegression with eXplanatory variables
# and Nested Random Effects
# Andrew Parnell

# Some jags code for fitting an ARIMAX-RE model.
# For the simpler ARIMAX model, see jags_ARIMAX.R
# For the simpler ARIMA model, see jags_ARIMA.R
# For the even simpler MA model, see jags_moving_average.R
# and for the similarly simple AR model, see jags_autoregressive.R
# Throughout this code I assume no differencing, so it is really an ARMAX model

# Some boiler plate code to clear the workspace, and load in required packages
rm(list=ls()) # Clear the workspace
library(R2jags)

# Maths -------------------------------------------------------------------

# Description of the Bayesian model fitted in this file
# Notation:
# y_{ij}(t) = response variable for region j and site i at time t
# i=1,...,N_region, j=1,...,N_site, t=1,...,T
# mu = mean parameter
# sigma_y = residual standard deviation
# sigma_t = time series standard deviation
# alpha = mean parameter
# phi = AR parameters
# b1 = random effect for region, with associated sd sigma_region
# b2 = random effect for site within region, with associated sd sigma_site_reg
# d = number of first differences
# p = number of autoregressive components
# We do the differencing outside the model so let z[t] = diff(y, differnces = d)
# k = number of explanatory variables
# beta = regression parameters
# x_j = explanatory variables, a T by k by N_region array

# Likelihood:
# y[site[t],region[t],t] = mu[t] + b1[region[t]] + b2[site[t],region[t]] + gamma[site[t],region[t],1:k] %*% x[region[t],t,1:k] + e[i,j,t]
# e[i,j,t] ~ N(0, sigma_y^2)
# mu[t] = alpha + phi*mu[t-1] + eps[t], eps[t] ~ N(0, sigma_mu^2)
# b1 ~ N(0, sigma_region)
# b2 ~ N(0, sigma_region_site)
# Priors - all vague here
# alpha ~ N(0,100)
# phi ~ N(0,100)
# gamma ~ N(0, 100)
# sigma_y ~ unif(0,10)
# sigma_mu ~ unif(0,10)
# beta ~ N(0,100)
# sigma_region ~ unif(0,10)
# sigma_site_region ~ unif(0,10)

# Simulate data -----------------------------------------------------------

# Some R code to simulate data from the above model
p = 1 # Number of autoregressive terms
k = 2 # Number of explanatory variables
T = 200
N_site = 10
N_region = 3
set.seed(123)
region = sort(sample(1:N_region, size=T, replace=TRUE))
site = sort(sample(1:N_site, size=T, replace=TRUE))
sigma_y = 1
sigma_t = 1
sigma_region = 0.5
sigma_site_region = 0.8
alpha = 0
gamma = matrix(rnorm(N_region*N_site, 0, 5), ncol=N_region,nrow=N_site)
phi = sort(runif(p),decreasing=TRUE)
y = rep(NA,T)
x = array(rnorm(T*k*N_region),dim = c(N_region,T,k))
b1 = rnorm(N_region, 0, sigma_region)
b2 = matrix(rnorm(N_region*N_site, 0, sigma_region),ncol=N_region,N_site)
mu = rep(NA, T)
mu[1] = 0
stop()
for(t in (p+1):T) mu[t] = rnorm(1, alpha + sum( phi * mu[(t-1):(t-p)] ), sigma_t)
for(t in 1:T) {
  reg_mean = sum( x[t,]*beta )  
}

y[t] = rnorm(1, mean = alpha + ar_mean + ma_mean + reg_mean, sd = sigma)
plot(1:T,y,type='l')

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code = '
model
{
  # Set up residuals
  for(t in 1:max(p,q)) {
    eps[t] <- z[t] - alpha
  }
  # Likelihood
  for (t in (max(p,q)+1):T) {
    z[t] ~ dnorm(alpha + ar_mean[t] + ma_mean[t] + reg_mean[t], tau)
    ma_mean[t] <- inprod(theta, eps[(t-q):(t-1)])
    ar_mean[t] <- inprod(phi, z[(t-p):(t-1)])
    reg_mean[t] <- inprod(beta, x[t,])
    eps[t] <- z[t] - alpha - ar_mean[t] - ma_mean[t] - reg_mean[t]
  }

  # Priors
  alpha ~ dnorm(0.0,0.01)
  for (i in 1:q) {
    theta[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:p) {
    phi[i] ~ dnorm(0.0,0.01)
  }
  for(i in 1:k) {
    beta[i] ~ dnorm(0.0,0.01)
  }
  tau <- 1/pow(sigma,2) # Turn precision into standard deviation
  sigma ~ dunif(0.0,10.0)
}
'

# Set up the data
model_data = list(T = T, z = y, x=x, q = 1, p = 1, k=2)

# Choose the parameters to watch
model_parameters =  c("alpha","theta","phi","beta","sigma")

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
print(model_run)


# Real example ------------------------------------------------------------

# Data wrangling and jags code to run the model on a real data set in the data directory
hadcrut = read.csv('https://raw.githubusercontent.com/andrewcparnell/tsme_course/master/data/hadcrut.csv')
head(hadcrut)
with(hadcrut,plot(Year,Anomaly,type='l'))

# Look at the ACF/PACF
acf(hadcrut$Anomaly)
pacf(hadcrut$Anomaly)

# Set up the data
real_data = with(hadcrut,
                 list(T = nrow(hadcrut),
                      z = Anomaly,
                      x = matrix(Year,ncol=1),
                      q = 1,
                      p = 1,
                      k = 1))

# This needs a longer run to get decent convergence
real_data_run = jags(data = real_data,
                     parameters.to.save = model_parameters,
                     model.file=textConnection(model_code),
                     n.chains=4,
                     n.iter=10000,
                     n.burnin=2000,
                     n.thin=8)

# Plot output
print(real_data_run) # beta small and convergence not as good

traceplot(real_data_run, mfrow=c(1,2), varname = 'beta', ask = FALSE)
hist(real_data_run$BUGSoutput$sims.list$beta, breaks=30)
par(mfrow=c(1,1))

# Create some predictions off into the future
# Using the trick first covered in the jags_ARIMA function
T_future = 20 # Number of future data points
year_future = (max(hadcrut$Year)+1):(max(hadcrut$Year)+T_future)

real_data_future = with(hadcrut,
                        list(T = nrow(hadcrut) + T_future,
                             z = c(Anomaly, rep(NA,T_future)),
                             x = matrix(c(Year,year_future),ncol=1),
                             q = 1,
                             p = 1,
                             k = 1))

# Just watch y now
model_parameters =  c("z")

# Run the model
real_data_run_future = jags(data = real_data_future,
                            parameters.to.save = model_parameters,
                            model.file=textConnection(model_code),
                            n.chains=4,
                            n.iter=10000,
                            n.burnin=2000,
                            n.thin=8)

# Print out the above
print(real_data_run_future)

# Get the future values
y_all = real_data_run_future$BUGSoutput$sims.list$z
# If you look at the above object you'll see that the first columns are all identical because they're the data
y_all_mean = apply(y_all,2,'mean')
# Also create the upper/lower 95% CI values
y_all_low = apply(y_all,2,'quantile',0.025)
y_all_high = apply(y_all,2,'quantile',0.975)
year_all = c(hadcrut$Year,year_future)

# Plot these all together
plot(year_all,
     y_all_mean,
     type='n',
     ylim=range(c(hadcrut$Anomaly,y_all_low,y_all_high)))
lines(year_all,y_all_mean,col='red')
lines(year_all,y_all_low,col='red',lty='dotted')
lines(year_all,y_all_high,col='red',lty='dotted')
with(hadcrut,lines(Year,Anomaly))

# Other tasks -------------------------------------------------------------

# Perhaps exercises, or other general remarks

# 1) It might be that a linear function of time is not a good model. Run a model which includes k=2 explanatory variables where the second one is a quadratic in year (hint: it might be a good idea to divide year by 1000 before you start to avoid numerical overflow)
# 2) A linear term in year seems like a bit of a waste of the ARIMAX model. What else might be useful to include as an explanatory variable? See if you can find worldwide yearly values for whatever you choose and include it as a covariate
# 3) (Harder) An ARIMAX prediction at time t is made up of four components: the overall mean, the ar terms, the ma terms, and the regression terms. Save these components individually in the 'parameters.to.save' argument and create a plot to show how important they are across the time series and how they behave in the future projections.

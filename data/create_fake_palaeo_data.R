# Create some fake palaeoclimate data

# Preamble
rm(list=ls())
#setwd("~/GitHub/tsme_course/")
library(R2jags)

# Load in the hadcrut data set
temp = read.csv('data/hadcrut.csv')

# Need to create a time grid from 1000 AD to 2015 AD with fill in temps from 1855-2015 AD and filled in proxy values for the entire range
# Notation: t is time, x is global mean temp, y is proxy data
t_all = 1000:2015
i_all = t_all - 999 # turn into indices for easy access later
N = length(t_all)
t_mod = temp$Year
x_mod = temp$Anomaly
x_all = x_all_hidden = rep(NA,length(t_all))
x_all[t_all %in% t_mod] = x_mod
x_all_hidden[t_all %in% t_mod] = x_mod

# Simulate old temperatures
sigma = 0.01
set.seed(123)
for(i in rev(i_all)) {
  if(is.na(x_all[i])) {
    x_all[i] = rnorm(1, x_all[i+1], sigma)
  }
}

plot(t_all, x_all)

# Now create proxy values
alpha_y = 4.2
beta_y = 3.5
sigma_y = 0.2
y_all = rnorm(length(t_all), alpha_y + beta_y * x_all, sigma_y)

plot(t_all, y_all)

df = data.frame('year'=t_all,'proxy'=round(y_all,3),'temp'=x_all_hidden)

write.csv(df,file='data/palaeo.csv',quote=FALSE,row.names=FALSE)

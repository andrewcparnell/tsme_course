rm (list = ls( )) 
########Required Packages and Libraries########
###Uncomment to install packages
#install.packages("rjags")
library(rjags)
#install.packages("R2jags")
library(R2jags)
#install.packages("coda")
library(coda)
#install.packages("runjags")
library(runjags)
#install.packages("MASS")
library(MASS)
#install.packages("fields")
library(fields)
##############################################################

#######Set your working directory#######
setwd("~")

########Read in and set up the data########
data<-read.csv("~.csv")
N<-nrow(data)
x1<-data[,1] #If x data is in column 1    
x<-x1/1000
y<-data[,2] #If x data is in column 2
plot(x,y) #Plot the data

###################################One Change Point Model#################################
CP1model="
model
{
for(i in 1:N)
{
y[i] ~ dnorm(mu[i],tau)
mu[i]<-alpha + beta[J[i]]*(x[i]-x.change)

J[i] <- 1 + step(x[i] - x.change)
}

###Priors
alpha ~ dnorm(0.0,1.0E-6) #Uninformative

beta[1]~ dnorm(0.0,1.0E-6) #Uninformative
beta[2]~ dnorm(0.0,1.0E-6) #Uninformative

x.change~dunif(1.88,2.014) #CP prior is uniform over the range of the data

tau <- 1/pow(sigmay,2) ##JAGS works with precision instead of variance
sigmay ~ dgamma(4,1)

}##End model
"

##The required data
mydata <- list(N=N,y=y,x=x) 
##Set up initial values for parameters
myinitial<-function(){list("alpha"=rnorm(1,0,3),
                           "beta"=(rnorm(2,0,3)),
                           "x.change"=runif(1,1.88,2.014),
                           "sigmay"=runif(1,0,10) )}
##Parameters to look at
mypars <- c("alpha","beta","x.change","sigmay")  

#Set up (this may be changed)
myiterations <- 20000
myburnin <-2000
mythin <-3

##Run the model
run.mcmc1 <- jags(data=mydata, inits=myinitial, parameters.to.save=mypars, model.file=textConnection(CP1model),n.chains=1, n.iter=myiterations, n.burnin=myburnin,n.thin=mythin,DIC=TRUE)

##Look at the results and summary statistics
mcmcmodel1 <- run.mcmc1$BUGSoutput$sims.list
results1.jags <- mcmc(matrix(unlist(run.mcmc1$BUGSoutput$sims.list),ncol=length(run.mcmc1$BUGSoutput$sims.list)+1))
colnames(results1.jags) <- c("alpha","beta1","beta2",names(run.mcmc1$BUGSoutput$sims.list)[-2][-1])
plot(results1.jags)
sumstat<-summary(results1.jags)

####Save the results
save(run.mcmc1,file="~run.mcmc1.RData")
save(results1.jags,file="~results1jags.RData")
stop()

###Diagnostics
acfplot(results1.jags)
geweke.plot(results1.jags)
run.mcmc1$BUGSoutput$DIC

###################################Two Change Point Model#################################
CP2model="
model
{
for(i in 1:N)
{
y[i] ~ dnorm(mu[i],tau)
mu[i]<-alpha[K[i]]+beta[J[i]+K[i]]*(x[i]-x.change[K[i]])

J[i] <- step(x[i]-x.change[1])
K[i]<-1+step(x[i]-x.change[2])
}

###Priors
alpha[1] ~ dnorm(0.0,1.0E-6)
alpha[2] ~ dnorm(0.0,1.0E-6)

beta[1]~dnorm(0.0,1.0E-6)
beta[2]<-(alpha[2]-alpha[1])/(x.change[2]-x.change[1]) #Deterministic
beta[3]~dnorm(0.0,1.0E-6)

x.change.temp[1] ~ dunif(1.88,2.014)
x.change.temp[2] ~ dunif(1.88,2.014)
x.change[1:2]<-sort(x.change.temp) #CPs need to be ordered

tau <- 1/pow(sigmay,2)
sigmay ~ dgamma(4,1)

}##End model
"

##The required data
mydata <- list(N=N,x=x,y=y)  
##Set up initial values for parameters
myinitial<-function(){list("alpha"=rnorm(2,0,3),
                           "beta"=c(rnorm(1,0,3),NA,rnorm(1,0,3)),
                           "x.change.temp"=sort(runif(2,1.88,2.014)),
                           "sigmay"=runif(1,0,10))}
##Parameters to look at
mypars <- c("alpha","beta","x.change","sigmay")  

#Set up 
myiterations <- 20000
myburnin <-2000
mythin <-3

##Run the model
run.mcmc2 <- jags(data=mydata, inits=myinitial, parameters.to.save=mypars, model.file=textConnection(CP2model),n.chains=1, n.iter=myiterations, n.burnin=myburnin,n.thin=mythin,DIC=TRUE)

##Look at the results and summary statistics
mcmcmodel2 <- run.mcmc2$BUGSoutput$sims.list
results2.jags <- mcmc(matrix(unlist(run.mcmc2$BUGSoutput$sims.list),ncol=length(run.mcmc2$BUGSoutput$sims.list)+4))
colnames(results2.jags) <- c("alpha1","alpha2","beta1","beta2","beta3","deviance","sigmay","xchange1","xchange2",names(run.mcmc2$BUGSoutput$sims.list)[-2][-1][-2][-1][-1][-1])
plot(results2.jags)
sumstat<-summary(results2.jags)

####Save the results
save(run.mcmc2,file="~run.mcmc2.RData")
save(results2.jags,file="~results2jags.RData")
stop()

###Diagnostics
acfplot(results2.jags)
geweke.plot(results2.jags)
run.mcmc2$BUGSoutput$DIC

###################################Three Change Point Model#################################
CP3model="
model
{
for(i in 1:N)
{
y[i] ~ dnorm(mu[i],tau)
mu[i] <- alpha[K[i]+L[i]] + beta[J[i]+K[i]+L[i]]*(x[i]-x.change[K[i]+L[i]])

J[i] <- step(x[i] - x.change[1])
K[i]<-step(x[i]-x.change[2])
L[i]<-1 + step(x[i]-x.change[3])

}

##Priors
alpha[1]~dnorm(0.0,1.0E-6)
alpha[2]~dnorm(0.0,1.0E-6)
alpha[3]~dnorm(0.0,1.0E-6)

beta[1]~dnorm(0.0,1.0E-6)
beta[2]<-(alpha[2]-alpha[1])/(x.change[2]-x.change[1])
beta[3]<-(alpha[3]-alpha[2])/(x.change[3]-x.change[2])
beta[4]~dnorm(0.0,1.0E-6)

x.change.temp[1] ~ dunif(1.88,2.014)
x.change.temp[2] ~ dunif(1.88,2.014)
x.change.temp[3] ~ dunif(1.88,2.014)
x.change[1:3]<-sort(x.change.temp)

tau <- 1/pow(sigmay,2)
sigmay ~ dgamma(4,1)

}##End model
"

##The required data
mydata <- list(N=N,x=x,y=y) 
##Set up initial values for parameters
myinitial<-function(){list("alpha"=rnorm(3,0,3),
                           "beta"=c(rnorm(1,0,3),NA,NA,rnorm(1,0,3)),
                           "x.change.temp"=sort(runif(3,1.88,1.998)),
                           "sigmay"=runif(1,0,10))}

##Parameters to look at
mypars <- c("alpha","beta","x.change","sigmay")  
 
#Set up 
myiterations <- 20000
myburnin <-2000
mythin <-3

##Run the model
run.mcmc3 <- jags(data=mydata, inits=myinitial, parameters.to.save=mypars, model.file=textConnection(CP3model),n.chains=1, n.iter=myiterations, n.burnin=myburnin,n.thin=mythin,DIC=TRUE)

##Look at the results and summary statistics
mcmcmodel3 <- run.mcmc3$BUGSoutput$sims.list
results3.jags <- mcmc(matrix(unlist(run.mcmc3$BUGSoutput$sims.list),ncol=length(run.mcmc3$BUGSoutput$sims.list)+7))
colnames(results3.jags) <- c("alpha1","alpha2","alpha3","beta1","beta2","beta3","beta4","deviance","sigmay","xchange1","xchange2","xchange3",names(run.mcmc3$BUGSoutput$sims.list)[-1][-1][-5][-1][-1][-1][-1])
plot(results3.jags)
sumstat<-summary(results3.jags)

####Save the results
save(run.mcmc3,file="~run.mcmc3.RData")
save(results3.jags,file="~results3jags.RData")
stop()

###Diagnostics
acfplot(results3.jags)
geweke.plot(results3.jags)
run.mcmc3$BUGSoutput$DIC

###################################Four Change Point#################################
CP4model="
model
{
for(i in 1:N)
{
y[i] ~ dnorm(mu[i],tau)
mu[i]<-alpha[K[i]+L[i]+M[i]]+beta[J[i]+K[i]+L[i]+M[i]]*(x[i]-x.change[K[i]+L[i]+M[i]])

J[i]<-step(x[i]-x.change[1])
K[i]<-step(x[i]-x.change[2])
L[i]<-step(x[i]-x.change[3])
M[i]<-1+step(x[i]-x.change[4])
}

alpha[1]~dnorm(0.0,1.0E-6)
alpha[2]~dnorm(0.0,1.0E-6)
alpha[3]~dnorm(0.0,1.0E-6)
alpha[4]~dnorm(0.0,1.0E-6)

beta[1]~dnorm(0.0,1.0E-6)
beta[2]<-(alpha[2]-alpha[1])/(x.change[2]-x.change[1])
beta[3]<-(alpha[3]-alpha[2])/(x.change[3]-x.change[2])
beta[4]<-(alpha[4]-alpha[3])/(x.change[4]-x.change[3])
beta[5]~dnorm(0.0,1.0E-6)

x.change.temp[1] ~ dunif(1.88,2.014)
x.change.temp[2] ~ dunif(1.88,2.014)
x.change.temp[3] ~ dunif(1.88,2.014)
x.change.temp[4] ~ dunif(1.88,2.014)
x.change[1:4]<-sort(x.change.temp)

tau <- 1/pow(sigmay,2)
sigmay ~ dgamma(2,1)

}##End model
"

##The required data
mydata <- list(N=N,x=x,y=y)  
##Set up initial values for parameters
myinitial<-function(){list("alpha"=rnorm(4,0,3),
                           "beta"=c(rnorm(1,0,3),NA,NA,NA,rnorm(1,0,3)),
                           "x.change.temp"=sort(runif(4,1.88,2.01)),
                           "sigmay"=runif(1,0,10))}
##Parameters to look at
mypars <- c("alpha","beta","x.change","sigmay")  

#Set up 
myiterations <- 20000
myburnin <-2000
mythin <-3

##Run the model
run.mcmc4 <- jags(data=mydata, inits=myinitial, parameters.to.save=mypars, model.file=textConnection(CP4model),n.chains=1, n.iter=myiterations, n.burnin=myburnin,n.thin=mythin,DIC=TRUE)

##Look at the results and summary statistics
mcmcmodel4 <- run.mcmc4$BUGSoutput$sims.list
results4.jags <- mcmc(matrix(unlist(run.mcmc4$BUGSoutput$sims.list),ncol=length(run.mcmc4$BUGSoutput$sims.list)+10))
colnames(results4.jags) <- c("alpha1","alpha2","alpha3","alpha4","beta1","beta2","beta3","beta4","beta5","deviance","sigmay","xchange1","xchange2","xchange3","xchange4",names(run.mcmc4$BUGSoutput$sims.list)[-2][-1][-2][-1][-1])
plot(results4.jags)
sumstat<-summary(results4.jags)

####Save the results
save(run.mcmc3,file="~run.mcmc3.RData")
save(results3.jags,file="~results3jags.RData")
stop()

###Diagnostics
acfplot(results4.jags)
geweke.plot(results4.jags)
run.mcmc4$BUGSoutput$DIC


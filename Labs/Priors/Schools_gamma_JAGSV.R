# JAGS model: uniform prior on sigma
model 
{
for (j in 1:J){ # J = the number of schools
  y[j] ~ dnorm (theta[j], tau.y[j]) # data model: the likelihood
  theta[j] <- mu.theta + eta[j]
  tau.y[j] <- pow(sigma.y[j], -2)
}

for(j in 1:J){
    eta[j]~dnorm(0,tau.theta)
}
#priors

mu.theta ~ dnorm(0,1.0E-6)
#uniform====
#tau.theta <- pow(sigma.theta, -2)
#sigma.theta ~ dunif(0,upper.sigma)  #no miscalibration from 100 to 1000
tau.theta ~ dgamma(.001,.001)
sigma.theta <-1/sqrt(tau.theta)
}

    
    



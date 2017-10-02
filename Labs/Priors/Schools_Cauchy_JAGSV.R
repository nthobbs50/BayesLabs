model 
{
for (j in 1:J){ # J = the number of schools
  y[j] ~ dnorm (theta[j], tau.y[j]) # data model: the likelihood
  theta[j] <- mu.theta + xi*eta[j]
  tau.y[j] <- pow(sigma.y[j], -2)
}
xi ~ dnorm (0, tau.xi)
tau.xi <- pow(prior.scale, -2)
for (j in 1:J){
  eta[j] ~ dnorm (0, tau.eta) # hierarchical model for theta
}
sigma.theta <- abs(xi) / sqrt(tau.eta)
tau.eta ~ dgamma(.5, 5)
#tau.eta ~ dgamma (.5, .5) # chi^2 with 1 d.f.
#sigma.theta <- abs(xi)/sqrt(tau.eta) # cauchy = normal/sqrt(chi^2)
mu.theta ~ dnorm (0.0, 1.0E-6) # noninformative prior on mu
}
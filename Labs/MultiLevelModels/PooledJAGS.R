
model{
#priors
alpha ~ dnorm(0,.0001)
beta ~ dnorm(0,.0001)
sigma ~ dunif(0,100)
tau.reg <- 1/sigma^2
#likelihood
#note that data have been log-transformed on the R side. 
 for(i in 1:length(y.emission)){
    log_mu[i] <- alpha + beta * y.n.input[i]
    y.emission[i] ~ dnorm(log_mu[i], tau.reg)
 }

}
    


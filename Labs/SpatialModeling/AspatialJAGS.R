
    model{
    for(i in 1:y.p){
      beta[i] ~ dnorm(0,.00001)
    }
    sigma ~ dunif(0,5)
    tau <- 1/sigma^2
    #mu = X %*% beta
    #prec.mat = inverse(sigma^2*y.I) #y.I = identity matrix
    #y ~ dmnorm(mu,prec.mat)
    for(i in 1:length(y)){
      mu[i] = beta[1] + beta[2]*X[i,2]
      y[i] ~ dnorm(mu[i], tau)
      y.new[i] ~ dnorm(mu[i], tau)
    }
    p.sd=step(sd(y.new) - sd(y))
    p.mean = step(mean(y.new[]) - mean(y[]))
    epsilon = y - mu
    
}
    

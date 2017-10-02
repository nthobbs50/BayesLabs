
    model{
    for(i in 1:y.p){
    beta[i] ~ dnorm(0,.00001)
    }
    sigma ~ dunif(0,5)  #unstructured error
    #tau ~ dgamma(2,5)
    tau ~ dgamma(.01,.10)
    sigma_sq_struc <- 1 / tau 
    phi ~ dunif(3/y.max.prior.range, 3/y.min.prior.range)
    mu = X %*% beta
    #Define exponential variance covariance matrix
    for(i in 1:length(y)){
       for(j in 1:length(y)){
        Cov.mat[i,j] <- sigma_sq_struc * exp(-D[i,j] * phi)
      }
    }
    prec.mat = inverse(Cov.mat + sigma^2*y.I)
    y[] ~ dmnorm(mu[],prec.mat) 
    epsilon = y - mu
    }
    

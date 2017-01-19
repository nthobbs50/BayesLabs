
model{
  z <- X %*% beta # the regression model in matrix form, returns a vector of length
  for(i in 1:n)   { 
    y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(z[i])
    #calculate predicitve density for use in WAIC
    pd[i] <- dpois(y[i],lambda[i]) #note when the lhs of the <- is not data, dpois() returns a probability
    #calculate the log predicitve density for use in WAIC
    log_pd[i] <- log(dpois(y[i],lambda[i]))
    #simulate new data sets for posterior predictive loss
    y.new[i]~dpois(lambda[i])
  }
  
  
  # PRIORS
  # p = number of coefficients, including intercept
  for(i in 1:p) {  
    beta[i] ~ dnorm(0, 0.01)
  }
}


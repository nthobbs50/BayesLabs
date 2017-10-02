
model {
  
  # LIKELIHOOD
  # n= number of states
  # y = number of birds in each state
  
  for(i in 1:n)   { 
    y[i] ~ dpois(lambda[i])
    z[i] <- beta
    lambda[i] <- exp(z[i])
    #calculate predicitve density for use in WAIC
    pd[i] <- dpois(y[i],lambda[i]) #note when the lhs of the <- is not data, dpois() returns a probability
    #calculate the log predicitve density for use in WAIC
    log_pd[i] <- log(dpois(y[i],lambda[i]))
    #simulate new data sets for posterior predictive loss
    y.new[i]~dpois(lambda[i])
  }
  
  
  # PRIOR

    beta ~ dnorm(0, 0.01)

  
}

    


model{
  
  ### priors 
  tau ~ dgamma(0.001, 0.001)
  alpha ~ dnorm(0, 1e-4)

  # derived quantities
  sd2 <- 1 / tau
  
  ### likelihood

  # because we are taking the log of the product alpha*x_ct, we need to restrict values to be       
  # positive. by looping through all values and choosing the maximum between a very small number positive   
  # number and each alpha*x_c[i], we restrict all values to be non-zero and positive

   for (i in 1:n){
     mu[i] <- max(1e-10, x_ct[i] * alpha)
   }
  
  # remember, y_lnwt is already log-transformed
  for (i in 1:n){
    y_lnwt[i] ~ dnorm(log(mu[i]), tau)
  }
  
  
  }# end model
    


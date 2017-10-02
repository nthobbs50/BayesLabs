
 model{
  ## priors
  # vague priors
  tau.p ~ dgamma(0.001, 0.001)
  sigmasq.p <- 1 / tau.p
  phi ~ dunif(1e-10, 5)
  for(i in 1:p){
     beta[i]~dnorm(0,0.0001)
   }
   
  # informed priors
   tau.c  ~ dgamma(17.1841127, 0.2118985)
   alpha.c ~ dgamma(2505.777, 9639.005)
  

  ## create covariance matrix
  # account for spatial autocorrelation using a covariance matrix specified by
  # distance matrix (D), structured variance (sigma.p.str^2), and Matern variable (phi)
   
   for(i in 1:n){
     for(j in 1:n){
       Cov.mat[i,j]<-sigmasq.p^2*exp(-D[i,j]*phi)
     }
   }
   
  # invert covariance matrix to obtain precision matrix (use in multivariate normal)
    Prec.mat <- inverse(Cov.mat)
   
  # mean number of seeds governed by beta priors and covariates
    mu <- X%*%beta
  
  # restrict mu to be positive  (or you could exponentiate)
   for (i in 1:n){
     mu.pos[i]<- max(0.00000001, mu[i])
   }
   
  # multivariate lognormal (all positive values of seed number) with precision matrix
   lnz ~ dmnorm(log(mu.pos), Prec.mat)
   
  # likelihood (data conditioned on log(alpha.c*z[i]) which is the same as
  # log(alpha.c) + log(z[i])
    for (i in 1:n){
      y[i] ~ dnorm(log(alpha.c) + lnz[i], tau.c)
      y.new[i] ~ dnorm(log(alpha.c) + lnz[i], tau.c)
      
      # squared residuals for model fit
      sq.y[i] <- (y[i] - (log(alpha.c) + lnz[i])) ^ 2
      sq.ynew[i] <- (y.new[i] - (log(alpha.c) + lnz[i])) ^ 2
    }
      # posterior predictive check
   
      # sum of squared residuals (from above) for model fit
      sumsq.y <- sum(sq.y[])
      sumsq.ynew <- sum(sq.ynew[])
      p.sumsq <- step(sumsq.ynew - sumsq.y)
      
 }


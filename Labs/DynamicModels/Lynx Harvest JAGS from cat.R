
model{
#priors==============
sigma.p ~ dunif(0,5)
tau.p <- 1/sigma.p^2
lambda ~ dunif(0,10)
p ~ dbeta(y.a,y.b)  #Get parameters a and b from mean and sd using moment matching to make this prior informative.  These are calcuated on R side and read in as data.

#Informative priors on initial conditions based on first year's observation of family groups
fg[1] ~ dpois(y[1])
N[1]~dlnorm(log(y[1]/p),tau.p)


#process model===============
for(t in 2:(y.endyr + 1)){  # the last year is a forecast with known harvest data
	mu[t] <- log(max(.0001,lambda*(N[t-1]-y.H[t-1])))
	N[t] ~ dlnorm(mu[t], tau.p)
	fg[t]<-N[t] * p
	}#end of process model
		
#data model===============
for(t in 2:y.endyr){   
		y[t] ~ dpois(p*N[t])  
	    	}  #end of data model
	
#simulate new data for posterior predicitve check
for(t in 1:y.endyr){
	    y.rep[t] ~ dpois(p*N[t])
	    #accumlate test statistics for posterior predictive check
	    sq[t] <- (y[t]-p*N[t])^2
	    sq.rep[t] <-(y.rep[t] - p*N[t])^2
}
#calculate Bayesian P value
fit <- sum(sq[])
fit.new <- sum(sq.rep[])
pvalue <- step(fit.new - fit)

##forecast effects of different harvest regeimes on next year's number of family grops
	for(i in 1:length(h)){
		#mu.hat is the forecast 1 year beyond y.endyr +1, i.e., 2011
		mu.hat[i] <- log(max(.001,lambda*(N[y.endyr+1]-h[i]))) 
		N.hat[i] ~ dlnorm(mu.hat[i], tau.p)	#Nhat forecasts 2 years out
		fg.hat[i] <- N.hat[i] * p
    }
	
	
} #end of model



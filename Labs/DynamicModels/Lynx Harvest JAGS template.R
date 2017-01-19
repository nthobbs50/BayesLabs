# file: "/Users/Tom/Documents/Ecological Modeling Course/_A_Master_Lab_Exercises/Dynamic models_lynx_problem/Lynx Harvest JAGS.R"

#This code does estimation of poplation state and population growth rate for lynx in Sweden. Modified by Tom Hobbs on May 25, 2015

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
#Include code for process model here. (3 lines)
for(t in 2:(y.endyr + 1)){  # the last year is a *forecast* of N at time y.endyr + 1 with known harvest data
		
		}#end of process model
		
#data model===============
#Include code for data model here (1 line)
for(t in 2:y.endyr){   
	
	    	}  #end of data model
	
#simulate new data for posterior predicitve check
for(t in 1:y.endyr){
	  
}
#calculate Bayesian P value

##forecast effects of different harvest regeimes on next year's number of family grops
	    }
	
	
} #end of model.  
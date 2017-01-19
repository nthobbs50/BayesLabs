# file: "/Users/Tom/Documents/Ecological Modeling Course/_A_Master_Lab_Exercises/Dynamic models_lynx_problem/Lynx Harvest JAGS.R"
rm(list=ls())
library(rjags)
setwd("/Users/Tom/Documents/Ecological Modeling Course/_A_Master_Lab_Exercises/Dynamic models_lynx_problem/")
y=read.csv("Lynx data.csv")

#Function to get beta shape parameters from moments
shape_from_stats <- function(mu = mu.global, sigma = sigma.global){
		 a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
		 b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
		shape_ps <- c(a,b)
		return(shape_ps)
}



#get shape parameters for population multiplier, 1/p
shapes=shape_from_stats(.163,.012)
#check prior on p using simulated data from beta distribution
x = seq(0,1,.001)
p=dbeta(x,shapes[1],shapes[2])
plot(x,p,typ="l",xlim=c(.1,.3))


#visually estimate some data for initial conditions
endyr = nrow(y)
n=numeric(endyr+1)
mu=numeric(endyr+1) #use this for family groups
lambda=1.02
sigma.p=.00001
n[1] = y$census[1]

for(t in 2: (endyr+1)){
	n[t] <- lambda*(y$census[t-1] - .16 * y$harvest[t-1])  #use this for family groups
	}
plot(y$census,ylim=c(0,100))
lines(n)


# Levels of  Harvest to evaluate relative to goals.
h=c(0, 10, 25, 50, 75)

#Data for JAGS
data = list(
	y.endyr = endyr,
	y.a=shapes[1], 
	y.b=shapes[2],
	y.H=y$harvest,
	y=y$census,
	h=h
)




inits = list(
	list(
	lambda = 1.2,
	sigma.p = .01,
	N=n
	),
	list(
	lambda = 1.01,
	sigma.p = .2,
	N=n*1.2),
	list(
	lambda = .95,
	sigma.p = .5,
	N=n*.5
	))


model = "Lynx Harvest JAGS.R"

n.update=10000
n.iter=50000
n.adapt=5000
n.thin=1

jm = jags.model(model,data=data,inits=inits, n.adapt=n.adapt, n.chains=length(inits))
update(jm, n.iter=n.update)

z = coda.samples(jm,variable.names=c("lambda","sigma.p","p"), n.iter=n.iter, thin=n.thin)

zj=jags.samples(jm,variable.names=c("N","N.hat", "fg", "fg.hat", "pvalue", "fit", "fit.new"), n.iter=n.iter, thin=n.thin)
summary(z)

#check convergence

#look at Bayesian P value


#Do goodness of fit plot.
par(mfrow=c(1,1))
plot(zj$fit.new,zj$fit, xlab = "Discrepancy observed", ylab= "Discrepancy simulated", cex=.05, xlim=c(0,3000), ylim=c(0,3000))
abline(0,1)
p=summary(zj$pvalue,mean)$stat
text(500,2500, paste("Bayesian P value = ", as.character(signif(p,2))))

#Plot quantiles of the true, unobserved state vs observations and forecasts.



#calculate probability of meeting goals






	
	


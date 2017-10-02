#preliminaries
rm(list=ls())
library(rjags)
library(reshape)
library(ggplot2)
set.seed(5)
setwd("/Users/Tom/Documents/Ecological Modeling Course/_A_Master_Lab_Exercises/Multi-level models NO2/")
y=read.csv("NO_2 emission all data for exercise.csv")
y.n.sites = length(unique(y$group))
qplot(n.input, emission, data=y, color =  group)

group_from_index = function(group, group.index, output ){
  #group is a vector of group names or numbers
  #group.index is vector of indices to group names or numbers.  It is a sequence of integers 1 to length(group)
  #output is a matrix or dataframe of output with number of rows = length(group). Each row contains statistics, etc for each group.
  a = unique(as.vector(group)) 
  b = unique(group.index)
  group.key=as.data.frame(t(rbind(a,b))) #columns containing indices paired with group name or number
  names(group.key)= c(names(as.data.frame(group)), names(as.data.frame(group.index))) 
  link.column.name = names(group.key)[2] #name of column for merging output data with groups
  output2 = cbind(seq(1,nrow(output)),output) #give the output data sequential index the same as 
  colnames(output2)[1]=link.column.name
  group.data=as.data.frame(merge(group.key, output2, by = link.column.name )) #merge the output with the groups
  return(group.data)
}
#data for all models except last one
data = list(
  y.emission = log(y$emission),
  y.n.input = log(y$n.input) - mean(log(y$n.input)), #center the data to speed convergence and aid in interpretation. Can recover 0 intercept if needed.
  j = y$group.index,  #use j to index groups
  k= y$fert.index, #use k to index fertilizer types
  y.n.sites = length(unique(y$group))
)

#Revise code to includ model checks

####Pooled model
sink("Pooled")
cat("
model{
#priors
alpha ~ dnorm(0,.0001)
beta ~ dnorm(0,.0001)
sigma ~ dunif(0,100)
tau.reg <- 1/sigma^2
#likelihood
 for(i in 1:length(y.emission)){
    mu[i] <- alpha + beta * y.n.input[i]
    y.emission[i] ~ dnorm(mu[i], tau.reg)
    #simulated data for posterior predictive checks
    
    #create vectors for discrepancy statistic
   
 }
#Bayesian P values
    
",fill=TRUE)
sink()

inits = list(
  list(
    alpha = 0,
    beta = .5,
    sigma = 50
  ),
  list(
    alpha = 1,
    beta = 1.5,
    sigma = 10
  )
)

n.update=10000
n.iter=10000
n.update = 3000
jm.pooled = jags.model("Pooled", data=data, n.adapt = 3000, inits=inits, n.chains=length(inits))
update(jm.pooled, n.iter = n.update)
#zc.pooled = coda.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma"), n.iter=n.iter)
zj.pooled = jags.samples(jm.pooled, variable.names = c("alpha", "beta", "sigma", "p.sd", "p.mean", "p.discrep", "y.emission.sim"), n.iter=n.iter)

hist(data$y.emission, breaks=20, freq=FALSE, main="Simulated and real data", xlab=expression(paste("log(", N0[2], " emission)")), cex.lab=1.2) #note that this is the log transformed data
lines(density(zj.pooled$y.emission.sim), col="red")
legend(-10,.2,"simulated", col="red", lty="solid")
plot(density(zj.pooled$beta))
summary(zc.pooled)
gelman.diag(zc.pooled)

#Bayesian p values here:


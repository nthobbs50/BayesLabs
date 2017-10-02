rm(list=ls())

pdf_flag=FALSE

#This is fashioned from the example in Gelman 2006, Prior distributions for variance parameters in
#hierarchical models
#schools=read.table( "/Users/Tom/Documents/SESYNC/Bayes course August 2015/Gelman and Hill code and data/ARM_Data/schools/schools.dat", header=T)
#write.csv(schools, file="schools.csv")

setwd("/Users/Tom/Documents/Sweden Visiting Professor/Oct 2015 visit/Workshop/")
schools=read.csv("Priors/schools.csv")



J <- nrow (schools)
y <- schools$estimate
sigma.y <- schools$sd
prior.scale <- 25
data<-list(J=J,y=y,sigma.y=sigma.y,upper.sigma=100)

inits=list(
  list(mu.theta=0,sigma.theta=1),
  list( mu.theta=0,sigma.theta=10),
  list( mu.theta=0,sigma.theta=100))

library(rjags)
n.iter=10000
n.update=10000
n.adapt=3000

#Run this code with uniform prior on sigma.theta with all 8 schools
library(rjags)
jm=jags.model("Priors/Schools_uniform_JAGSV.R",  data=data,n.chains=3, n.adapt=n.adapt)
update(jm, n.iter=n.update)
zcu=coda.samples(jm,variable.names=c("theta", "mu.theta", "sigma.theta"),n.iter=n.iter)
zju=jags.samples(jm, variable.names=c("sigma.theta"), n.iter=n.iter)
plot(zcu)
summary(zcu)
gelman.diag(zcu)
if(pdf_flag) pdf(file="Uniform_8_schools.pdf", height=5, width=5)
hist(zju$sigma.theta,freq=FALSE,xlim=c(0,30), main=expression(paste("MCMC ouptut, uniform prior on ",sigma[theta])), xlab=expression(sigma[theta]), cex.lab=1.5,cex.main=1.5,ylim=c(0,.15),breaks=50, sub="8 schools")
x=seq(0,40)
y=dunif(x,0,100)
lines(x,y,col="red",lwd=2)
legend(8,.10,c("uniform prior"),col=c("red"), lty="solid", bty="n", lwd=2)
if(pdf_flag) dev.off()
p20u=ecdf(zju$sigma.theta)(20)


#8 schools with gamma(.001,.001) on precision
data<-list(J=J,y=y,sigma.y=sigma.y)

inits=list(
  list(mu.theta=0,tau.theta=1),
  list( mu.theta=0,tau.theta=1/10^2),
  list( mu.theta=0,tau.theta=1/100^2))

library(rjags)
n.iter=10000
n.update=100000
n.adapt=3000

library(rjags)
jm=jags.model("Priors/Schools_gamma_JAGSV.R",  data=data,n.chains=3, n.adapt=n.adapt)
update(jm, n.iter=n.update*2)
zcg=coda.samples(jm,variable.names=c("theta", "mu.theta", "sigma.theta"),n.iter=n.iter)
zjg=jags.samples(jm, variable.names=c("sigma.theta", "tau.theta"), n.iter=n.iter)
#plot(zcg)
summary(zcg)
gelman.diag(zcg)
par(mfrow=c(1,1))
if(pdf_flag) pdf(file="Gamma_8_schools.pdf", height=5, width=5)
hist(zjg$sigma.theta,freq=FALSE, main=expression(paste("MCMC ouptut, gamma prior on ",tau)), xlab=expression(sigma[theta]), cex.lab=1.5,cex.main=1.5, xlim=c(0,30),breaks=50, sub="8 schools")
if(pdf_flag) dev.off()

#now use three schools and uniform prior
schools=schools[1:3,]
y <- schools$estimate
sigma.y <- schools$sd
J=length(y)
data<-list(J=J,y=y,sigma.y=sigma.y,upper.sigma=1000)
inits=list(
  list(mu.theta=0,sigma.theta=1),
  list( mu.theta=0,sigma.theta=10),
  list( mu.theta=0,sigma.theta=100))

library(rjags)


library(rjags)
jm=jags.model("Priors/Schools_uniform_JAGSV.R",  data=data,n.chains=3, n.adapt=n.adapt)
update(jm, n.iter=n.update)
zcu=coda.samples(jm,variable.names=c("theta", "mu.theta", "sigma.theta"),n.iter=n.iter)
zju=jags.samples(jm, variable.names=c("sigma.theta"), n.iter=n.iter)
#plot(zc)
#summary(zcu)
#gelman.diag(zcu)
if(pdf_flag) pdf(file="Uniform_3_schools.pdf", height=5, width=5)
hist(zju$sigma.theta,freq=FALSE,xlim=c(0,200), main=expression(paste("MCMC ouptut, uniform prior on ",sigma[theta])), xlab=expression(sigma[theta]), cex.lab=1.5,cex.main=1.5, breaks=200, sub="3 schools")
x=seq(0,200)
y=dunif(x,0,1000)
lines(x,y,col="red",lwd=2)
legend(40,.010,c("uniform prior"),col=c("red"), lty="solid", bty="n", lwd=2)
if(pdf_flag) dev.off()
p20u=ecdf(zju$sigma.theta)(20)



#set to duplicate Gelman 2006 figure
n.iter=1000
n.update=1000
n.adapt=3000

#Now change the priors to use half-cauchy
data<-list(J=J,y=y,sigma.y=sigma.y, prior.scale=25)
jm=jags.model("Priors/Schools_Cauchy_JAGSV.R", data=data,n.chains=3, n.adapt=n.adapt)
update(jm, n.iter=n.update)
zcc=coda.samples(jm,variable.names=c("theta", "mu.theta", "sigma.theta"),n.iter=n.iter*10)
zjc=jags.samples(jm, variable.names=c("sigma.theta"), n.iter=n.iter*10)
if(pdf_flag) pdf(file="Cauchy_3_schools.pdf", height=5, width=5)
hist(zjc$sigma.theta,freq=FALSE,xlim=c(0,200), main=expression(paste("MCMC ouptut, half-Cauchy prior on ",sigma[theta])), xlab=expression(sigma[theta]), cex.lab=1.3,cex.main=1.3, breaks=50,sub="3 schools")
x=seq(0,200)
y=2*dcauchy(x,0,25)
lines(x,y,typ="l",xlim=c(0,200),col="red")
legend(50,.02, c("proper half-Cauchy prior"), col="red", lty="solid", bty="n")
if(pdf_flag)   dev.off()



####
####  Load Libraries
####

library(geoR)
library(maps)
library(mvtnorm)
library(rjags)
library(spBayes)
library(rgdal)
library(maptools)
library(dplyr)
setwd("/Users/Tom/Documents/Teaching/BayesWorkshop/Labs/SpatialModeling/")
gpclibPermit()

####
####  Change Projection to UTM
####

ut.st=map("state",regions="utah",fill=T,plot=F)
ut.coords=cbind(ut.st$x,ut.st$y)
plot(ut.coords,type="b")
ut.coords=rbind(ut.coords,ut.coords[1,])
ut.latlon=SpatialPoints(ut.coords,proj4string=CRS("+proj=longlat +datum=WGS84"))
str(ut.latlon)
ut.utm=spTransform(ut.latlon,CRS("+proj=utm +zone=12  +datum=WGS84"))
plot(ut.utm)

ut.range=ut.utm@bbox
xg=seq(ut.range[1],ut.range[3],25)
yg=seq(ut.range[2],ut.range[4],35)
utgrid.locs=as.matrix(expand.grid(xg,yg))
plot(utgrid.locs,pch=20,asp=1,xlab="x.utms",ylab="y.utms")
points(ut.utm@coords,type="l",lwd=3)
ntot=dim(utgrid.locs)[1]

####
####  Simulate Geostatistical process epsilon
####

D=as.matrix(dist(utgrid.locs))
s2=2
phi=1.5*10^-5
#phi=1  # to simulate absence of spatial structure
Sigma=s2*exp(-D*phi)

plot(seq(0,max(D),,20),s2*exp(-seq(0,max(D),,20)*phi),type="o",ylab="cov",xlab="distance")

set.seed(13)
eps=as.vector(rmvnorm(1,matrix(0,ntot,1),Sigma,method="chol"))     # may take some time
image(matrix(eps,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)))
myImagePlot(matrix(eps,length(xg),length(yg)))
image(matrix(eps,length(xg),length(yg)),x=xg,y=yg,col=heat.colors(50))
points(ut.utm@coords,type="l",lwd=3)

idxkeep=sort(sample(1:ntot,round(0.2*ntot)))
epsmask=matrix(0,ntot,1)
epsmask[idxkeep,]=1
image(matrix(eps,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),asp=TRUE)
image(matrix(epsmask,length(xg),length(yg)),x=xg,y=yg,col=c("white","transparent"),asp=TRUE,add=TRUE)
points(ut.utm@coords,type="l",lwd=3)

####
####  Examine empirical spatial structure in simulated process
####

n=length(idxkeep)
utgrid.sm.locs=utgrid.locs[idxkeep,]
eps.gd=as.geodata(cbind(utgrid.sm.locs,eps[idxkeep]))

eps.v=variog(eps.gd,max.dist=max(D))
plot(eps.v)
lines(eps.v$uvec,s2-s2*exp(-eps.v$uvec*phi),type="o",pch=20)

####
####  Simulate (or specify) spatial covariates and then calculate data 
####

p=3
X=matrix(1,ntot,p)
X[,2]=-cos(scale(utgrid.locs[,1])-.5)*cos(scale(utgrid.locs[,2]))
X[,3]=scale(utgrid.locs[,1])+scale(utgrid.locs[,2])

beta=c(1,-2,1)
y=X%*%beta+eps

layout(matrix(1:4,2,2))
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y",asp=TRUE)
image(matrix(epsmask,length(xg),length(yg)),x=xg,y=yg,col=c("white","transparent"),asp=TRUE,add=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(eps,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="epsilon",asp=TRUE)
image(matrix(epsmask,length(xg),length(yg)),x=xg,y=yg,col=c("white","transparent"),asp=TRUE,add=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(X[,2],length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main=bquote(x[1]),asp=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(X[,3],length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main=bquote(x[2]),asp=TRUE)
points(ut.utm@coords,type="l",lwd=3)

####
####  Fit Bayesian Geostatistical Model to Data using JAGS
####
####  Note: These functions will take a while to run.
####

y.sm=y[idxkeep]
X.sm=X[idxkeep,]
D.sm=D[idxkeep,idxkeep]
min.prior.range=max(D.sm)/500
max.prior.range=max(D.sm)
data=list(y = y.sm, X = X.sm, n = n, p = p, D = D.sm,min.prior.range=min.prior.range,max.prior.range=max.prior.range)
inits <- list(beta = rep(0,p), tau = .5, phi = 2*10^-5 )
var.names <- c("beta","phi","tau")

jags.sp.m <- jags.model(file = "lab_8_jags_model.txt", data = data, inits = inits, n.chains = 1, n.adapt = 0)
jags.out <- jags.samples(model=jags.sp.m,variable.names=var.names,n.iter=1000)

zc=coda.samples(jags.sp.m, variable.names=var.names, n.iter = 1000)

par(mfrow=c(2,3))
for(i in 1:p){
  plot(jags.out$beta[i,,1],type="l",col=i,main=paste(c("beta",i-1)))
  abline(h=beta[i],lwd=3)
}
plot(jags.out$phi[1,,1],type="l",col=4,main="phi")
abline(h=phi,lwd=3)
plot(1/jags.out$tau[1,,1],type="l",col=5,main="Sigma^2")
abline(h=s2,lwd=3)

####
####  Fit Bayesian Geostatistical Model to Data using the spBayes package
####
####  Note: spBayes is much faster than JAGS for spatial modeling. 
####  These functions will take a while to run.
####

starting <- list("phi"=4/max(D), "sigma.sq"=1)
tuning <- list("phi"=.1, "sigma.sq"=.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/max.prior.range, 3/min.prior.range),"sigma.sq.IG"=c(2, 5))
cov.model <- "exponential"

m.1 <- spLM(y.sm~0+X.sm, coords=utgrid.sm.locs, starting=starting, tuning=tuning,priors=priors, cov.model=cov.model, n.samples=10000)

#spDiag(m.1)

m.1.betas <- spRecover(m.1,get.w=FALSE)

layout(matrix(1:3,3,1))
plot(m.1.betas$p.theta.samples[,1],type="l",lty=1)  # Check convergence for covariance parameters
abline(v=c(s2),col=8,lwd=2)
plot(m.1.betas$p.theta.samples[,2],type="l",lty=1)  # Check convergence for covariance parameters
matplot(m.1.betas$p.beta.recover.samples,type="l",lty=1)  # Check convergence for covariance parameters
abline(h=beta,col=8,lwd=2)

####
####  Obtain Bayesian Spatial Predictions and Prediction SD
####

m.1.pred <- spPredict(m.1, pred.covars=X, pred.coords=utgrid.locs, start=2000)
y.pred=apply(m.1.pred$p.y.predictive.samples, 1, mean)
y.pred.sd=apply(m.1.pred$p.y.predictive.samples, 1, sd)

layout(matrix(1:4,2,2))
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y truth",asp=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(y,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y data",asp=TRUE)
image(matrix(epsmask,length(xg),length(yg)),x=xg,y=yg,col=c("white","transparent"),asp=TRUE,add=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(y.pred,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y prediction",asp=TRUE)
points(ut.utm@coords,type="l",lwd=3)
image(matrix(y.pred.sd,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="prediction sd",asp=TRUE)
points(ut.utm@coords,type="l",lwd=3)

##########################################
##########################################
####
####  CHALLENGE:  Iowa Temperature 
####
##########################################
##########################################


####
####  Read in Data
####

iowa.df=read.table("iowa_temperature.txt")
names(iowa.df)=c("long","lat","temp")
head(iowa.df)

####
####  Plot Data in Default Lat/Long Projection 
####

plot(iowa.df[,1:2],cex=1+4*(iowa.df$temp-min(iowa.df$temp))/diff(range(iowa.df$temp)),type="p")
map("state",lwd=2,col=1,add=TRUE)

####
####  Convert to UTM 
####

ia.st=map("state",regions="iowa",fill=T,plot=F)
ia.coords=cbind(ia.st$x,ia.st$y)
ia.coords=rbind(ia.coords,ia.coords[1,])
ia.latlon=SpatialPoints(ia.coords,proj4string=CRS("+proj=longlat +datum=WGS84"))
ia.utm=spTransform(ia.latlon,CRS("+proj=utm +zone=15  +datum=WGS84"))
iowa.data.latlon=SpatialPoints(iowa.df[,1:2],proj4string=CRS("+proj=longlat +datum=WGS84"))
iowa.data.utm=spTransform(iowa.data.latlon,CRS("+proj=utm +zone=15  +datum=WGS84"))

ia.range=iowa.data.utm@bbox
ia.diff=apply(ia.range,1,diff)
n.x=50 
n.y=round(n.x/(ia.diff[1]/ia.diff[2]))
xg=seq(ia.range[1],ia.range[3],,n.x)
yg=seq(ia.range[2],ia.range[4],,n.y)
iagrid.locs=as.matrix(expand.grid(xg,yg))
plot(iagrid.locs,pch=20,asp=1,xlab="x.utms",ylab="y.utms")
lines(ia.utm@coords,lwd=3)
points(iowa.data.utm@coords,col=2,lwd=1)
n.grid=dim(iagrid.locs)[1]
n.data=dim(iowa.data.utm@coords)[1]
all.locs=rbind(iowa.data.utm@coords,iagrid.locs)
D.max=max(dist(all.locs))

####
####  Create Design Matrix using Locs as Covariates 
####

p=3
X=matrix(1,n.data+n.grid,p)
X[,2:3]=scale(all.locs)
y=as.vector(iowa.df$temp)

####
####  Fit Uncorrelated Error Model and View Residuals 
####

temp.resid=resid(lm(iowa.df$temp~0+X[1:n.data,]))
temp.gd=as.geodata(cbind(all.locs[1:n.data,],temp.resid))
temp.v=variog(temp.gd,max.dist=D.max/2)
plot(variog4(temp.gd,max.dist=D.max/2))
plot(temp.v)

####
####  Fit Bayesian Geostatistical Model to Data using the spBayes package
####
####  Note: These functions will take a while to run.
####

min.prior.range=60000
max.prior.range=D.max/2
starting <- list("phi"=10/D.max, "sigma.sq"=1,"tau.sq"=1)
tuning <- list("phi"=.1, "sigma.sq"=.1,"tau.sq"=.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/max.prior.range, 3/min.prior.range),"sigma.sq.IG"=c(2, 5),"tau.sq.IG"=c(2,5))
cov.model <- "exponential"

m.1 <- spLM(y~0+X[1:n.data,], coords=iowa.data.utm@coords, starting=starting, tuning=tuning,priors=priors, cov.model=cov.model, n.samples=10000)

m.1.betas <- spRecover(m.1,get.beta=TRUE,get.w=TRUE)
spDiag(m.1.betas)

layout(matrix(1:3,1,3))
plot(as.vector(m.1.betas$p.theta.samples[,1]),type="l",lty=1)  # Check convergence for covariance parameters
plot(as.vector(m.1.betas$p.theta.samples[,2]),type="l",lty=1)  # Check convergence for covariance parameters
matplot(m.1.betas$p.beta.recover.samples,type="l",lty=1)  # Check convergence for covariance parameters

####
####  Obtain Bayesian Spatial Predictions and Prediction SD
####

m.1.pred <- spPredict(m.1, pred.covars=X[-(1:n.data),], pred.coords=iagrid.locs, start=2000)
y.pred=apply(m.1.pred$p.y.predictive.samples, 1, mean)
y.pred.sd=apply(m.1.pred$p.y.predictive.samples, 1, sd)

layout(matrix(1:2,1,2))
image(matrix(y.pred,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y prediction",asp=TRUE)
points(iowa.data.utm@coords,lwd=2)
lines(ia.utm@coords,lwd=2)
image(matrix(y.pred.sd,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="prediction sd",asp=TRUE)
points(iowa.data.utm@coords,lwd=2)
lines(ia.utm@coords,lwd=2)

####
####  Create Design Matrix using long, as Covariates
####

p=5
X=matrix(1,n.data+n.grid,p)
X[,2:3]=scale(all.locs)
X[,4]=X[,2]^2
X[,5]=X[,2]^3
y=as.vector(iowa.df$temp)

####
####  Fit Uncorrelated Error Model and View Residuals
####

temp.resid=resid(lm(iowa.df$temp~0+X[1:n.data,]))
temp.gd=as.geodata(cbind(all.locs[1:n.data,],temp.resid))
temp.v=variog(temp.gd,max.dist=D.max/2)
plot(temp.v)

####
####  Fit Bayesian Geostatistical Model to Data using the spBayes package
####
####  Note: These functions will take a while to run.
####

min.prior.range=60000
max.prior.range=D.max/2
starting <- list("phi"=10/D.max, "sigma.sq"=1,"tau.sq"=1)
tuning <- list("phi"=.1, "sigma.sq"=.1,"tau.sq"=.1)
priors <- list("beta.Flat", "phi.Unif"=c(3/max.prior.range, 3/min.prior.range),"sigma.sq.IG"=c(2, 5),"tau.sq.IG"=c(2,5))
cov.model <- "exponential"

m.1 <- spLM(y~0+X[1:n.data,], coords=iowa.data.utm@coords, starting=starting, tuning=tuning,priors=priors, cov.model=cov.model, n.samples=10000)

m.1.betas <- spRecover(m.1,get.beta=TRUE,get.w=TRUE)
spDiag(m.1.betas)


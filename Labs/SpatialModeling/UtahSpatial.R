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
gpclibPermit()
setwd("/Users/Tom/Documents/Teaching/BayesWorkshop/Labs/SpatialModeling/")

par(mfrow=c(1,1))
####
####  Create spatial grid for Utah in utm
####
par(mfrow=c(1,1))
ut.st=map("state",regions="utah",fill=T,plot=F)
ut.coords=cbind(ut.st$x,ut.st$y)
plot(ut.coords,type="b")
ut.coords=rbind(ut.coords,ut.coords[1,])
ut.latlon=SpatialPoints(ut.coords,proj4string=CRS("+proj=longlat +datum=WGS84"))
str(ut.latlon)
ut.utm=spTransform(ut.latlon,CRS("+proj=utm +zone=12  +datum=WGS84"))
plot(ut.utm)

ut.range=ut.utm@bbox
xg=seq(ut.range[1],ut.range[3],,25)
yg=seq(ut.range[2],ut.range[4],,35)
utgrid.locs=as.matrix(expand.grid(xg,yg))
plot(utgrid.locs,pch=20,asp=1,xlab="x.utms",ylab="y.utms")
points(ut.utm@coords,type="l",lwd=3)
ntot=dim(utgrid.locs)[1]

#Make some spatially structured, standarized covariates and add them to the gridded data

X2=-cos(scale(utgrid.locs[,1])-.5)*cos(scale(utgrid.locs[,2]))
X3=scale(utgrid.locs[,1])+scale(utgrid.locs[,2])

dat = cbind(utgrid.locs,X2,X3)
colnames(dat) = c("easting", "northing", "X2", "X3")



#Make an unstructured response variable
p=3

X=matrix(1,ntot,p)
X[,2]=dat[,3]
X[,3]=dat[,4]
#Note that the matrix has a column of 1's
#set some coefficients
beta=c(1,-2,1)
sigma = .05
y=rnorm(ntot,X%*%beta,sigma)
dat = as.data.frame(cbind(y,dat))
head(dat)

##Look at the data
par(mfrow=c(2,2))
image(matrix(dat$y,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y",asp=TRUE, xlab="Easting", ylab="Northing")

image(matrix(dat$X2,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="x2",asp=TRUE,xlab="Easting", ylab="Northing")
#ImagePlot(matrix(dat$X2,length(xg),length(yg)))
image(matrix(dat$X3,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="x3",asp=TRUE,xlab="Easting", ylab="Northing")
image(matrix(rnorm(ntot,mean=0,sd=1),length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="unstructured y",asp=TRUE,xlab="Easting", ylab="Northing")

image.plot(rnorm(ntot,mean=0,sd=1))

#fit aspatial model
sink("AspatialJAGS.R")
cat("
    model{
    for(i in 1:y.p){
      beta[i] ~ dnorm(0,.00001)
    }
    sigma ~ dunif(0,.5)
    mu = X %*% beta
    prec.mat = inverse(sigma^2*y.I) #y.I = identity matrix
    y ~ dmnorm(mu,prec.mat)
    epsilon = y - mu
}
    ",fill=TRUE)
sink()
inits=list(
  list(
  beta=rep(1,p),
  sigma=.02
  )
)

#make index for subsetting data to increase computation speed
idxkeep=sort(sample(1:ntot,round(0.2*ntot)))
data = list(
  y=as.double(dat$y[idxkeep]),
  X = cbind(rep(1,idxkeep),dat[,4:5]), #Be sure you have a column of 1's in the matrix given to JAGS
  y.p = p,
  y.I = diag(length(idxkeep)) #The identity matrix
)
jm_a=jags.model("AspatialJAGS.R", n.adapt = 0, n.chains = 1, inits=inits, data=data)
#zc_a = coda.samples(jm_a, variable.names=c("beta","sigma"), n.iter=1000)

zj_a =jags.samples(jm_a,variable.names=c("beta","sigma","epsilon"), n.iter=1000)
summary(zj_a$beta,quantile,c(.025,.5,.975))
summary(zj_$sigma, quantile(c,.025,.5,.975))
epsilon = summary(zj_a$epsilon,mean)$stat
image(matrix(epsilon,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y",asp=TRUE, xlab="Easting", ylab="Northing")

n=length(idxkeep)
utgrid.sm.locs=utgrid.locs[idxkeep,]
epsilon.gd=as.geodata(cbind(utgrid.sm.locs,epsilon))
epsilon.v=variog(epsilon.gd,max.dist=max(D))
plot(epsilon.v)
#no spatial strcuture in the residuals

##add spatial structure to the full data set
D=as.matrix(dist(utgrid.locs))
s2=2
phi=1.5*10^-5
Sigma=s2*exp(-D*phi)

plot(seq(0,max(D),,20),s2*exp(-seq(0,max(D),,20)*phi),type="o",ylab="cov",xlab="distance")

set.seed(13)
eps=as.vector(rmvnorm(1,matrix(0,ntot,1),Sigma,method="chol"))     # may take some time
image(matrix(eps,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)))

dat$y_structured = dat$y+eps
image(matrix(dat$y_structured,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)))


#refit aspatial model with structured data
data = list(
  y=as.double(dat$y_structured[idxkeep]),
  X = X[idxkeep,],
  y.p = p,
  y.I = diag(length(idxkeep))
)
jm_a=jags.model("AspatialJAGS.R", n.adapt = 0, n.chains = 1, inits=inits, data=data)
zc_a = coda.samples(jm_a, variable.names=c("beta","sigma"), n.iter=1000)

zj_a =jags.samples(jm_a,variable.names=c("epsilon"), n.iter=1000)

epsilon = summary(zj_a$epsilon,median)$stat
n=length(idxkeep)
utgrid.sm.locs=utgrid.locs[idxkeep,]
epsilon.gd=as.geodata(cbind(utgrid.sm.locs,epsilon))
epsilon.v=variog(epsilon.gd,max.dist=max(D/4))
plot(epsilon.v)
lines(epsilon.v$uvec,s2-s2*exp(-epsilon.v$uvec*phi),type="o",pch=20)
legend(100000,.5,c("Fit","Generating"),pch=c(1,20))

##Fit model accounting for spatial depdnence
sink("SpatialJAGS.R")
cat("
    model{
    for(i in 1:y.p){
    beta[i] ~ dnorm(0,.00001)
    }
    sigma ~ dunif(0,.5)  #unstructured error
    tau ~ dgamma(2,5)
    #tau ~ dgamma(.01,.10)
    sigma_sq_struc <- 1 / tau 
    phi ~ dunif(3/y.max.prior.range, 3/y.min.prior.range)
    mu = X %*% beta
    #Define exponential variance covariance matrix
    for(i in 1:length(y)){
       for(j in 1:length(y)){
        Cov.mat[i,j] <- sigma_sq_struc * exp(-D[i,j] * phi)
      }
    }
    prec.mat = inverse(Cov.mat + sigma^2*y.I)
    y[] ~ dmnorm(mu[],prec.mat) 
    epsilon = y - mu
    }
    ",fill=TRUE)
sink()
inits=list(
  list(
    beta=rep(1,p),
    sigma=.02,
    phi = 2*10^-5,
    tau = .5 
  )
)
#increase number of data points
idxkeep=sort(sample(1:ntot,round(0.4*ntot)))
D.sm = D[idxkeep,idxkeep]
data = list(
  y=as.double(dat$y_structured[idxkeep]),
  X = X[idxkeep,],
  y.p = p,
  y.I = diag(length(idxkeep)),
  D = D.sm,
  y.min.prior.range=max(D.sm)/500,
  y.max.prior.range=max(D.sm)
)

n.iter=15000
jm_s=jags.model("SpatialJAGS.R", n.adapt = 1000, n.chains = 1, inits=inits, data=data)
update(jm_s, n.iter=n.iter)
#zc_s = coda.samples(jm_s, variable.names=c("beta","sigma","phi"), n.iter=n.iter)

zj_s =jags.samples(jm_s,variable.names=c("epsilon","beta","sigma","phi"), n.iter=n.iter)
epsilon_s = summary(zj_s$epsilon,mean)$stat
image(matrix(epsilon_s,length(xg),length(yg)),x=xg,y=yg,col=rev(rainbow(100,start=0,end=.7)),main="y",asp=TRUE, xlab="Easting", ylab="Northing")

n=length(idxkeep)
utgrid.sm.locs=utgrid.locs[idxkeep,]
epsilon.gd=as.geodata(cbind(utgrid.sm.locs,epsilon_s))
epsilon.v=variog(epsilon.gd,max.dist=max(D))
plot(epsilon.v)


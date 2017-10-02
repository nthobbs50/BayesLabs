z=rgamma(10000,2,5)
z2=1/sqrt(z)
hist(z2,breaks=500,xlim=c(1,6), xlab=expression(sigma[struc]),freq=FALSE, main=expression(paste("Prior on ", sigma[struc])))

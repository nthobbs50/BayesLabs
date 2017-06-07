
rm(list=ls())
alpha=38
gamma=1.7
c=5
sigma=2
set.seed(4)
par(mfrow=c(1,1))
L=sort(runif(50,min=12, max = 100))
mu = alpha*(L - c) / (alpha/gamma + (L-c))
plot(L,mu, typ="l")
y=rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
plot(L,y)
lines(L,mu)
model=nls(y ~ alpha*(L - c) / (alpha/gamma + (L-c)), start=list(alpha=50,gamma=4,c=2))
s=summary(model)
p=coef(model)
alpha=p[1]
gamma=p[2]
c=p[3]
y.hat=alpha*(L - c) / (alpha/gamma + (L-c))
lines(L,y.hat,col="red")
legend(40,18, c("generating", "nls fit"), lty=c("solid", "solid"), col=c("black", "red"), bty="n")
	





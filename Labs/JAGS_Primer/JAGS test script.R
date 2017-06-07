library(rjags)

sink("phi_transmission")
cat("
model{
phi ~ dbeta(1,1)
for(t in 1:3){
	y[t] ~ dbin((1-phi)^t, n[t])
	}
}
",fill=TRUE)
sink()


n=c(28,29,40)
y=c(26,25,31)


data=list(
y=y,
n=n
)
jm=jags.model("phi_transmission",n.adapt=1000, n.chains=4, data=data)
update(jm, n.iter=3000)
zm=coda.samples(jm,variable.names=c("phi"), n.iter=5000)
summary(zm)


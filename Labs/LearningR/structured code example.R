rm(list=ls())
#A function for the light limitation of trees
#============
g=function(alpha, gamma,c, L){
	mu =alpha*(L - c) / (alpha/gamma + (L-c))
return(mu)
}
#=============
#get_data() simulates data for a Michaelis–Menten model 
#It uses the function g() for the model
#=============
get_data=function(alpha, gamma,c, sigma){
	set.seed(4)
	par(mfrow=c(1,1))
	L=sort(runif(50,min=12, max = 100))
	mu = g(alpha=alpha, gamma=gamma, c=c, L=L)
	y=rgamma(length(mu), mu^2/sigma^2, mu/sigma^2)
	gen.y= mu
	return(list(y=y, x=L,gen.y=gen.y))
}
#=============
#Function to do non-linear least squared fit 
#x are the independent variables
#y are the dependent variables
#Guesses is a named list of initial values for the parameters to be estimated. 
#The order of parameters must be the same as the arguments to g( )
#================
fit = function(x,y, guesses){
	model=nls(y ~ g(alpha=alpha,gamma=gamma,c=c, L=x), start=guesses)
	s=summary(model)
	p=coef(model)
	y.hat=g(alpha=p[1],gamma=p[2],c=p[3], L=x)
	return(y.hat)
}

#================
###Function for plotting generating data and nls() fit
#data is the x and y data and gen.y is the generating line
#================
plot_fit = function(data,fit){
	plot(data$x, data$y, xlab="Light level (lumens)", ylab="Tree growth (cm / yr)")
	lines(data$x, data$gen.y)
	lines(data$x,fit,col="red")
	legend(40,18, c("generating", "nls fit"), lty=c("solid", "solid"), col=c("black", "red"), bty="n")
}

#=================
## Maestro function
#==================
run_functions=function(){
	tree.data=get_data(alpha=38, gamma=1.7, c=5, sigma=2)
	nls.fit=fit(x=tree.data$x, y=tree.data$y, guesses=list(alpha=50,gamma=4, c=2))
	plot_fit(data=tree.data, fit=nls.fit)
}
#================

run_functions()


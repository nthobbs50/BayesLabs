library(rjags)
rm(list=ls())

##Summarize original data to deal with missing values--not done by students but instructive
setwd("/Users/Tom/Documents/Ecological Modeling Course/_A_Master_Lab_Exercises/Swiss breeding birds occupancy model/occupancy lab materials/")

obs=read.csv("wtmatrix.csv",header=TRUE,sep=",",na.strings=c("NA"))
obs=obs[,c("y.1","y.2","y.3","elev","forest")]
y<-as.matrix(obs[,c("y.1","y.2","y.3")])
n<-apply(!is.na(y),1,sum)  #nice way to sum non missing, n is number of visits
#M<-nrow(y)
y<-apply(y,1,sum,na.rm=TRUE) #y is number of times birds were observed

obs=as.data.frame(cbind(n,y,obs$elev,obs$forest))
names(obs)=c("number_visits", "number_detections", "elev", "forest")
write.csv(obs,file="Swiss BB data.csv")


######Student work starts here
obs=read.csv(file="Swiss BB data.csv")
#create vectors of standardized covariates
elev<-as.vector(scale(obs[,"elev"],center=TRUE))
forest<-as.vector(scale(obs[,"forest"],center=TRUE))
elev2<-elev*elev
M<-nrow(obs)

#JAGS code
sink("model.txt")
cat("
      model { 
      	
      	
        
       #Write code for likelihood and true state of patch here.  Also simulate new data for use in posterior predictive checks.
            for(i in 1:length(y)){ 
           		
           		
         } 
         
        #Find probablility of occupancy at the average elevation and forest cover
        
        
        #Do posterior predictive checks using mean and standard deviation as test statistics.
                      
        
       #Estimate optimum elevation at mean forest cover for  scaled values
       
       
       #back transform to get optimum elevation in meters
      
       #
       #Do derived quanties to allow plotting probability of occupancy as a function of elevation at the averaged forest cover.
       
       for(j in 1:length(elev.x)){
       	logit(psi.elev[j]) <- b0 + b1*(elev.x[j]-mu.elev)/sd.elev + b2*((elev.x[j]-mu.elev)/sd.elev)^2
		}
     
      }#end of model

      

",fill=TRUE)
sink()



#make some elevations for predicting occupancy
elev.x=seq(500,2500,100)
data <- list ( y=obs$number_detections,n=obs$number_visits,forest=forest,elev=elev,elev2=elev2, mu.elev=mean(obs$elev), sd.elev=sd(obs$elev),elev.x=elev.x)

#Notice that z *must* be initialized at 1.  If not the code will not run.
inits= list(
list( z=rep(1,M),b0=rnorm(1),p=runif(1),b1=rnorm(1),b2=rnorm(1),b3=rnorm(1)),
list ( z=rep(1,M),b0=rnorm(1),p=runif(1),b1=rnorm(1),b2=rnorm(1),b3=rnorm(1)) 
)
parameters <- c("p","b0","b1","b2","b3") # use this for variable list in coda object


###Compile JAGS model, do updates, obtain coda obect for parameters, and JAGS object for plotting posterior distribution of optimum elevation, Bayesian p value, and probability of occupancy as a function of elevation

#check convergence


#output derived quantities and do plots


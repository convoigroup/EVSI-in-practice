## Run baseline model
source("Colorectal Cancer Model.R")

# Update sample size for Menzies method
n.sim.menz<-2500

# Define likelihood function for Menzies method
likeli<-function(ya,lambda,g.){
  y<-ya[,1]
  a<-ya[,2]
  p<-1-exp(-lambda*a^g.)
  neg.ll<-prod(p^y*(1-p)^(1-y))
  return(as.numeric(neg.ll))
}

## EVSI estimation - Menzies
mu.X<-array(NA,dim=c(n.sim,length(n.)))
EVSI.men<-array(NA,dim=length(n.))
for(i in 1:length(n.)){
  for(j in 1:n.sim.menz){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Resimulate ages
    ages<-sample(ages.long,n.[i],replace=TRUE)
    mu<-1-exp(-lambda_1*ages^g)
    dat<-rbinom(n.[i],1,mu)
    #Set to a multiple function to use mapply instead of nested loop
    app<-function(lambda,g){
      r<-likeli(ya=cbind(dat,ages),lambda,g)
      return(r)
    }
    #The likelihood of data conditional on each parameter set.
    like<-mapply(app,lambda=df.psa.params.informed[1:n.sim.menz,1],g=df.psa.params.informed[1:n.sim.menz,2])
    weights<-like/sum(like)
    mu.X[j,i]<-weights%*%save[1:n.sim.menz]
  }
  EVSI.men[i]<-mean(pmax(0,mu.X[,i]),na.rm=TRUE)-max(0,mean(mu.X[,i],na.rm=TRUE))
}

### EVSI
EVSI.men
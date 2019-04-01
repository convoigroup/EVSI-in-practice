## Run baseline model
source("Colorectal Cancer Model.R")

## Define likelihood function to maximise
likeli.s<-function(params,y,a){
  lambda<-params[1]/1e4
  g.<-params[2]
  p<-1-exp(-lambda*a^g.)
  neg.ll<-(prod(p^y*(1-p)^(1-y)))
  return(as.numeric(neg.ll))
}

## EVSI estimation - Strong et al.
EVSI.strong<-array(NA,dim=length(n.))

i<-1
for(N in n.){
  sum.stats<-array(NA,dim=c(n.sim,3))
  for(j in 1:n.sim){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Simulate enrollement age for patients - resimulate each time to represent a random study of this design.
    ages<-sample(ages.long,N,replace=TRUE)
    mu<-1-exp(-lambda_1*ages^g)
    dat<-rbinom(N,1,mu)
    # Maximise likelihood to get summary statistics
    sum.stats[j,1:2]<- hjkb(c(lambda_1*1e4,g),likeli.s,y=dat,a=ages,
                            lower=c(0.00000012,1.22),upper=c(59.06,5.7),control=list(tol=1e-5,maximize=TRUE))$par*c(1e-4,1)
    sum.stats[j,3]<-pweibull(50,shape=sum.stats[j,2],
                             scale=sum.stats[j,1]^(-1/sum.stats[j,2]))
  }
  # Calculate regression curve
  strong.mod<-gam(INB~te(sum.stats[,1],sum.stats[,2],sum.stats[,3]))
  mean(pmax(fitted(strong.mod),0))-max(0,mean(fitted(strong.mod)))
  EVSI.strong[i]<-mean(pmax(fitted(strong.mod),0))-max(0,mean(fitted(strong.mod)))
  i<-i+1
}

### EVSI
EVSI.strong
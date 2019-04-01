set.seed(1234) # It is a good practice to set the seed so results can be replicated

# Run baseline cost-effectiveness analysis
source("Chemotherapy Model.R")

# Reduce PSA size for Menzies method
Size.Outer<-20000

pi2<-theta[1:Size.Outer,"pi1"]*theta[1:Size.Outer,"rho"]
recover.amb<--log(1-theta[1:Size.Outer,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[1:Size.Outer,"lambda.2.3.fix"])

# All model parameters
params<-cbind(theta[1:Size.Outer,"pi1"],
              pi2,
              theta[1:Size.Outer,"gamma"],
              theta[1:Size.Outer,"gamma2"],
              recover.amb,
              recover.hosp)

# Calculate EVPPI from full data and reduce size
y<-gam(INB ~ te(pi1,rho,gamma,gamma2,lambda.2.3.fix,lambda.1.3.fix,k=3),
       data=theta)$fitted.values[1:Size.Outer]

## EVSI estimate Menzies

# Likelihood function
likelihood<-function(D,theta){
  n<-150
  PSA.sim<-dim(theta)[1]
  colnames(D)<-c("SE1","SE2","amb","hosp","N.die","N.hosp")
  colnames(theta)<-c("pi1","pi2","gamma","gamma2","amb","hosp")
  l<-dbinom(D[,"SE1"],n,theta[,"pi1"],log = TRUE)+
    dbinom(D[,"SE2"],n,theta[,"pi2"],log=TRUE)+
    dbinom(D[,"N.hosp"],D[,"SE1"]+D[,"SE2"],theta[,"gamma"],log=TRUE)+
    dbinom(D[,"N.die"],D[,"N.hosp"],theta[,"gamma2"],log=TRUE)+
    log(theta[,"amb"])*(D[,"SE1"]+D[,"SE2"]-D[,"N.hosp"])-theta[,"amb"]*D[,"amb"]+
    log(theta[,"hosp"])*(D[,"N.hosp"]-D[,"N.die"])-theta[,"hosp"]*D[,"hosp"]
  return(exp(l))
}


prepost.M<-array(dim=Size.Outer)
for(i in 1:Size.Outer){
  # Generate the data
  n<-150
  X.SE1<-rbinom(1,n,params[i,1])
  X.SE2<-rbinom(1,n,pi2[i])
  X.N.hosp<-rbinom(1,X.SE1+X.SE2,params[i,3])
  X.N.die<-rbinom(1,X.N.hosp,params[i,4])
  N.amb<-X.SE1+X.SE2-X.N.hosp
  T.re.amb<-rexp(N.amb,recover.amb[i])
  X.amb<-sum(T.re.amb)
  N.hosp<-X.N.hosp-X.N.die
  T.re.hosp<-rexp(N.hosp,recover.hosp[i])
  X.hosp<-sum(T.re.hosp)
  # Create dataframe with data
  D<-as.data.frame(cbind(X.SE1,X.SE2,X.amb,X.hosp,X.N.die,X.N.hosp))  
  # Calculate likelihood
  ll<-likelihood(D,params)
  w<-ll/sum(ll)
  # Reweight fitted values
  prepost.M[i]<-w%*%y
}

# Calculate EVSI
evsi.Menzies<-mean(pmax(0,prepost.M))-max(0,mean(prepost.M))

## EVSI
evsi.Menzies
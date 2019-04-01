set.seed(1234) # It is a good practice to set the seed so results can be replicated

# Run baseline cost-effectiveness analysis
source("Chemotherapy Model.R")

# Extract the key parameters
pi1<-theta[,"pi1"]
pi2<-theta[,"pi1"]*theta[,"rho"]
gamma<-theta[,"gamma"]
gamma2<-theta[,"gamma2"]
recover.amb<--log(1-theta[,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[,"lambda.2.3.fix"])

### EVSI estimation Strong et al.
X.SE1<-X.SE2<-X.N.die<-X.N.hosp<-X.amb<-X.hosp<-array(dim=Size.Prior)
for(i in 1:Size.Prior){
  # Generate the data
  n<-150
  X.SE1[i]<-rbinom(1,n,pi1[i])
  X.SE2[i]<-rbinom(1,n,pi2[i])
  
  X.N.hosp[i]<-rbinom(1,X.SE1[i]+X.SE2[i],gamma[i])
  X.N.die[i]<-rbinom(1,X.N.hosp[i],gamma2[i])
  
  N.amb<-X.SE1+X.SE2-X.N.hosp
  T.re.amb<-rexp(N.amb[i],recover.amb[i])
  X.amb[i]<-sum(T.re.amb)
  N.hosp<-X.N.hosp-X.N.die
  T.re.hosp<-rexp(N.hosp[i],recover.hosp[i])
  X.hosp[i]<-sum(T.re.hosp)
}

# Create a data frame of summary statistics
dat<-as.data.frame(cbind(INB,X.amb,X.SE1,X.SE2,X.N.hosp,X.hosp,X.N.die))

# Find regression between summart statistics and incremental net benefit
prepost.s<-gam(INB~te(X.amb, X.SE1, X.SE2, X.N.hosp, X.hosp, X.N.die, k=3),data=dat)

EVSI.strong<-mean(pmax(0,prepost.s$fitted.values))-max(0,mean(prepost.s$fitted.values))

### EVSI
EVSI.strong

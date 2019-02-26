#Required packages
library(earth)
library(R2jags)
library(R2OpenBUGS)
library(BCEA)
library(mgcv)
# Performs the baseline cost-effectiveness analysis
# Estimates EVPPI using the Heath et al method - object names evppi.half

#Set Working Directory to Source File Location
source("Chemotherapy Model.R")
# Loads function used for the Jalal methods
source('predict_ga.R', encoding = 'WINDOWS-1252')

# Performs the baseline cost-effectiveness analysis for a specified willingness to pay.
wtp<-30000
NB<-e*wtp-c

# Matrix of parameters of interest from baseline model.
extra.lines<-(Size.Prior+1):dim(prior.model$BUGSoutput$sims.matrix)[1]
theta<-as.data.frame(prior.model$BUGSoutput$sims.matrix[-extra.lines,c("pi1","rho","gamma","gamma2","lambda.2.3.fix","lambda.1.3.fix","SE[1]","SE[2]")])
colnames(theta)<-c("pi1","rho","gamma","gamma2","lambda.2.3.fix","lambda.1.3.fix","SE1","SE2")    

## Number of simulations
n.sim        <- nrow(NB)
## Number of strategies
n.strategies <- ncol(NB)

## Find the incremental net benefit for each treatment
d.star<-which.max(apply(NB,2,mean))
INB<-NB-NB[,d.star]

### EVPI
evpi <- mean(apply(INB,1,max))
evpi

### EVPPI estimation
y<-INB[,-d.star]
lmm1<-gam(y~s(gamma)+s(gamma2)+s(lambda.2.3.fix)+s(lambda.1.3.fix)+s(SE1)+s(SE2),data=theta)

evppi <- mean(pmax(0,lmm1$fitted.values))
evppi
save<-lmm1$fitted.values
plot(y,y-save)

#Number of Side Effects
pi1<-theta[,"pi1"]
pi2<-theta[,"pi1"]*theta[,"rho"]

#Treatment of Side Effects
gamma<-theta[,"gamma"]
gamma2<-theta[,"gamma2"]
recover.amb<--log(1-theta[,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[,"lambda.2.3.fix"])


#Model
model.dat<-function(){
  for(i in 1:Size.Outer){
    # Side effects analysis
    num.se[i] ~ dbin(pi[1,i], num.pat)     # sampling distribution
    pi[1,i] ~ dbeta(1, 1)               # prior distribution
    rho[i] ~ dnorm(m.rho, tau.rho)
    pi[2,i] <- rho[i] * pi[1,i]
    X.SE1[i] ~ dbin(pi[1,i],n)
    X.SE2[i] ~ dbin(pi[2,i],n)
    
    for (t in 1:2) {
      SE[t,i] ~ dbin(pi[t,i], N)         # Expected number of patients with side effects
    }
    
    num.amb[i] ~ dbin(gamma[i], num.se[i])      # sampling distribution
    gamma[i] ~ dbeta(1, 1)                # prior distribution
    num.death[i] ~ dbin(gamma2[i],num.se[i]-num.amb[i])
    gamma2[i] ~ dbeta(1,4)
    X.N.hosp[i] ~ dbin(gamma[i],X.SE1[i]+X.SE2[i])
    X.N.die[i] ~ dbin(gamma2[i],X.N.hosp[i])
    
    lambda.1.3.fix[i] ~ dbeta(p1.1.3,p2.1.3)
    lambda.2.3.fix[i] ~ dbeta(p1.2.3,p2.2.3)
    
    recover.amb[i]<--log(1-lambda.1.3.fix[i])
    recover.hosp[i]<--log(1-lambda.2.3.fix[i])
    for(j in 1:N.amb[i]){
      T.re.amb[i,j] ~ dexp(recover.amb[i])}
    for(j in 1:N.hosp[i]){
      T.re.hosp[i,j] ~ dexp(recover.hosp[i])}
     
    
    #You go to hospital 
    lambda.1.2[i]<-gamma[i]/TH
    #OR you recover with high prob
    lambda.1.3[i]<-(1-lambda.1.2[i])*lambda.1.3.fix[i]
    #OR stay in the same state
    lambda.1.1[i]<-(1-lambda.1.3.fix[i])*(1-lambda.1.2[i])
    
    #Either you die
    lambda.2.4[i]<-gamma2[i]/TH
    #OR you recover with high prob
    lambda.2.3[i]<-(1-lambda.2.4[i])*lambda.2.3.fix[i]
    #OR you stay in the same state
    lambda.2.2[i]<-(1-lambda.2.3.fix[i])*(1-lambda.2.4[i])
  }
  
  # Costs
  # These are sampled direct from distributions, so no prior is needed.
  c.amb ~ dlnorm(m.amb, tau.amb)     # Cost of ambulatory care 
  c.hosp ~ dlnorm(m.hosp, tau.hosp)  # Cost of hospitalization
  c.death ~ dlnorm(m.death,tau.death)#Cost of death
  
  # Effects
  e.chemo ~ dbeta(p1.chemo,p2.chemo)
  e.amb ~ dbeta(p1.amb,p2.amb)
  e.hosp ~ dbeta(p1.hosp,p2.hosp)
  

  
} 

filein <- file.path(tempdir(),fileext="datmodel.txt")
write.model(model.dat,filein)

#####################Jalal and Alarid-Escudero#####################
###Estimating n0
Size.Outer<-1000
Size.Inner<-10000
start<-Sys.time()

# Set the number of chains in iterations
n.chains <- 3     # Number of chains
n.burnin <- 1000  # Number of burn in iterations
n.iter <- ceiling(Size.Inner/n.chains) + n.burnin # Number of iterations per chain

X.SE1<-X.SE2<-X.N.hosp<-X.N.die<-N.amb<-N.hosp<-array(NA,dim=Size.Outer)
n<-30
T.re.amb<-T.re.hosp<-array(NA,dim=c(Size.Outer,2*n))
for(i in 1:Size.Outer){

  X.SE1[i]<-rbinom(1,n,pi1[i])
  X.SE2[i]<-rbinom(1,n,pi2[i])
  
  X.N.hosp[i]<-rbinom(1,X.SE1[i]+X.SE2[i],gamma[i])
  X.N.die[i]<-rbinom(1,X.N.hosp[i],gamma2[i])
  
  N.amb[i]<-X.SE1[i]+X.SE2[i]-X.N.hosp[i]
  if(N.amb[i]>0){
    T.re.amb[i,1:N.amb[i]]<-rexp(N.amb[i],recover.amb[i])
    }
  N.hosp[i]<-X.N.hosp[i]-X.N.die[i]
  if(N.hosp[i]>0){
    T.re.hosp[i,1:N.hosp[i]]<-rexp(N.hosp[i],recover.hosp[i])}
}
T.re.amb<-T.re.amb[,1:max(N.amb)]
T.re.hosp<-T.re.hosp[,1:max(N.hosp)]
  
data.dat <- list("X.SE1","X.SE2","n",
                   "X.N.hosp","X.N.die",
                   "N.amb","T.re.amb","N.hosp","T.re.hosp","Size.Outer")
#num.se<-rep(num.se,Size.Outer)
#num.amb<-rep(num.amb,Size.Outer)
#num.death<-rep(num.death,Size.Outer)

data.full<-append(data,data.dat)
#parameters.to.save<-parameters.to.save[-c(1,2)]
#parameters.to.save<-c("pi",parameters.to.save)

  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data.full,
    inits = NULL,
    parameters.to.save = parameters.to.save,
    model.file = filein, 
    n.chains = n.chains, 
    n.iter = n.iter, 
    n.thin = 1, 
    n.burnin = n.burnin) 
    
  gamma.mean<-apply(bugs.data$BUGSoutput$sims.list$gamma,2,mean)
  gamma2.mean<-apply(bugs.data$BUGSoutput$sims.list$gamma2,2,mean)
  lambda.2.3.fix.mean<-apply(bugs.data$BUGSoutput$sims.list$lambda.2.3.fix,2,mean)
  lambda.1.3.fix.mean<-apply(bugs.data$BUGSoutput$sims.list$lambda.1.3.fix,2,mean)
  SE.mean<-apply(bugs.data$BUGSoutput$sims.list$SE,2,mean)
  SE1.mean<-apply(bugs.data$BUGSoutput$sims.list$SE[,1,],2,mean)
  SE2.mean<-apply(bugs.data$BUGSoutput$sims.list$SE[,2,],2,mean)

##Estimating n0 separately for each element of the meta-model
n<-30
var.theta<-apply(theta,2,var)
#pi1.n0<-n*(var.theta[1]/var(pi1.mean)-1)
#rho.n0<-n*(var.theta[2]/var(rho.mean)-1)
gamma.n0<-n*(var.theta[3]/var(gamma.mean)-1)
gamma2.n0<-n*(var.theta[4]/var(gamma2.mean)-1)
lambda2.n0<-n*(var.theta[5]/var(lambda.2.3.fix.mean)-1)
lambda1.n0<-n*(var.theta[6]/var(lambda.1.3.fix.mean)-1)
SE1.n0<-n*(var.theta[7]/var(SE1.mean)-1)
SE2.n0<-n*(var.theta[8]/var(SE2.mean)-1)

n0<-c(gamma.n0,gamma2.n0,lambda2.n0,lambda1.n0,SE1.n0,SE2.n0)
n<-rep(150,6)
llpred<-predict.ga(lmm1,n0=n0,n=n)
evsi.Jal <- mean(pmax(0,llpred))-max(mean(llpred),0)
evsi.Jal
end<-Sys.time()

end-start
  
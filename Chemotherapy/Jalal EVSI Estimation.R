set.seed(1234) # It is a good practice to set the seed so results can be replicated

# Run baseline cost-effectiveness analysis
source("Chemotherapy Model.R")

# Function to run Jalal et al method
source('predict_ga.R', encoding = 'WINDOWS-1252')

### EVPPI estimation
lmm1<-gam(INB~s(pi1)+s(rho)+s(gamma)+s(gamma2)+s(lambda.2.3.fix)+s(lambda.1.3.fix),data=theta)

# EVPPI
mean(pmax(0, lmm1$fitted.values)) - max(0, mean(lmm1$fitted.values))

# Bayesian Model
model.dat<-function(){
  X.SE1 ~ dbin(pi[1],n)
  X.SE2 ~ dbin(pi[2],n)
  X.N.hosp ~ dbinom(gamma,X.SE1+X.SE2)
  X.N.die ~ dbin(gamma2,X.N.hosp)
  
  #DATA - Recovery
  recover.amb<--log(1-lambda.1.3.fix)
  recover.hosp<--log(1-lambda.2.3.fix)
  for(i in 1:N.amb){
    T.re.amb[i] ~ dexp(recover.amb)}
  for(i in 1:N.hosp){
    T.re.hosp[i] ~ dexp(recover.hosp)}
  
  
  # Side effects analysis
  num.se ~ dbin(pi[1], num.pat)     # sampling distribution
  pi[1] ~ dbeta(1, 1)               # prior distribution
  # rho = Relative risk of hematological side effects given new treatment compared to SoC
  # pi[2] = Probability of hematological side effects given new treatment
  rho ~ dnorm(m.rho, tau.rho)
  pi[2] <- rho * pi[1]
  
  # Treatment of side effects analysis -  Markov Model
  # gamma = probability that side effects are treated by ambulatory care
  # 1- gamma = probability that side effects are treated in hospital
  
  num.amb ~ dbin(gamma, num.se)      # sampling distribution
  gamma ~ dbeta(1, 1)                # prior distribution
  
  num.death ~ dbin(gamma2,num.se-num.amb)
  gamma2 ~ dbeta(1,4)
  
  #Recover
  lambda.1.3.fix ~ dbeta(p1.1.3,p2.1.3)
  lambda.2.3.fix ~ dbeta(p1.2.3,p2.2.3)
  
  #You go to hospital 
  lambda.1.2<-gamma/TH
  #OR you recover with high prob
  lambda.1.3<-(1-lambda.1.2)*lambda.1.3.fix
  #OR stay in the same state
  lambda.1.1 <-(1-lambda.1.3.fix)*(1-lambda.1.2)
  
  #Either you die
  lambda.2.4<-gamma2/TH
  #OR you recover with high prob
  lambda.2.3<-(1-lambda.2.4)*lambda.2.3.fix
  #OR you stay in the same state
  lambda.2.2<-(1-lambda.2.3.fix)*(1-lambda.2.4)
  
  # Costs
  # These are sampled direct from distributions, so no prior is needed.
  c.amb ~ dlnorm(m.amb, tau.amb)     # Cost of ambulatory care 
  c.hosp ~ dlnorm(m.hosp, tau.hosp)  # Cost of hospitalization
  c.death ~ dlnorm(m.death,tau.death)#Cost of death
  
  # Effects
  e.chemo ~ dbeta(p1.chemo,p2.chemo)
  e.amb ~ dbeta(p1.amb,p2.amb)
  e.hosp ~ dbeta(p1.hosp,p2.hosp)
  
  # Predictive distributions on the clinical outcomes
  for (t in 1:2) {
    SE[t] ~ dbin(pi[t], N)         # Expected number of patients with side effects
  }
  
  #VoI
  pi1<-pi[1]
  pi2<-pi[2]
  
} 

filein <- file.path(tempdir(),fileext="datmodel.txt")
write.model(model.dat,filein)

### EVSI Estimation from Jalal et al.
###Estimating n0
Size.Outer<-1000
Size.Inner<-10000

# Set the number of chains in iterations
n.chains <- 3     # Number of chains
n.burnin <- 1000  # Number of burn in iterations
n.iter <- ceiling(Size.Inner/n.chains) + n.burnin # Number of iterations per chain

# Parameters to generate data
pi2<-theta[,"pi1"]*theta[,"rho"]
recover.amb<--log(1-theta[,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[,"lambda.2.3.fix"])

pi1.mean<-rho.mean<-gamma.mean<-gamma2.mean<-lambda.2.3.fix.mean<-lambda.1.3.fix.mean<-SE1.mean<-SE2.mean<-array(NA,dim=Size.Outer)
for(i in 1:Size.Outer){
  # Generate Data
  # Change sample size as the same n0 can be used for each sample size and it is faster to run 
  # an economic model with n<-30
  n<-30
  X.SE1<-rbinom(1,n,theta[i,"pi1"])
  X.SE2<-rbinom(1,n,pi2[i])
  
  X.N.hosp<-rbinom(1,X.SE1+X.SE2,theta[i,"gamma"])
  X.N.die<-rbinom(1,X.N.hosp,theta[i,"gamma2"])
  
  N.amb<-X.SE1+X.SE2-X.N.hosp
  T.re.amb<-rexp(N.amb,recover.amb[i])
  N.hosp<-X.N.hosp-X.N.die
  T.re.hosp<-rexp(N.hosp,recover.hosp[i])
  
  data.dat <- list("X.SE1","X.SE2","n",
                   "X.N.hosp","X.N.die",
                   "N.amb","T.re.amb","N.hosp","T.re.hosp")
  
  data.full<-append(data,data.dat)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data.full,
    inits = inits,
    parameters.to.save = parameters.to.save,
    model.file = filein, 
    n.chains = n.chains, 
    n.iter = n.iter, 
    n.thin = 1, 
    n.burnin = n.burnin) 
  
  pi1.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$pi1)
  rho.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$rho)
  gamma.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$gamma)
  gamma2.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$gamma2)
  lambda.2.3.fix.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$lambda.2.3.fix)
  lambda.1.3.fix.mean[i]<-mean(bugs.data$BUGSoutput$sims.list$lambda.1.3.fix)
  SE.mean<-apply(bugs.data$BUGSoutput$sims.list$SE,2,mean)
  SE1.mean[i]<-SE.mean[1]
  SE2.mean[i]<-SE.mean[2]
}

## Estimating n0 separately for each element of the meta-model
n<-30
var.theta<-apply(theta,2,var)
pi1.n0<-n*(var.theta[1]/var(pi1.mean)-1)
rho.n0<-n*(var.theta[2]/var(rho.mean)-1)
gamma.n0<-n*(var.theta[3]/var(gamma.mean)-1)
gamma2.n0<-n*(var.theta[4]/var(gamma2.mean)-1)
lambda2.n0<-n*(var.theta[5]/var(lambda.2.3.fix.mean)-1)
lambda1.n0<-n*(var.theta[6]/var(lambda.1.3.fix.mean)-1)

n0<-c(pi1.n0,rho.n0,gamma.n0,gamma2.n0,lambda2.n0,lambda1.n0)
n<-rep(150,6)

# Estimate EVSI from predict.ga function.
llpred<-predict.ga(lmm1,n0 = n0,n=n)
evsi.Jal <- mean(pmax(0,llpred))-max(mean(llpred),0)

## EVSI
evsi.Jal

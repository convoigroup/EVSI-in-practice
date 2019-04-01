set.seed(1234) # It is a good practice to set the seed so results can be replicated

# Run baseline cost-effectiveness analysis
source("Chemotherapy Model.R")

### EVPPI estimation
save<-gam(INB~te(pi1,rho,gamma,gamma2,lambda.2.3.fix,lambda.1.3.fix,k=3),data=theta)$fitted.values

### EVPPI
mean(pmax(save,0)) - max(0, mean(save))

### EVSI Estimation - Heath et al.
# Model with Data
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

# Number of Side Effects
pi1<-theta[,"pi1"]
pi2<-theta[,"pi1"]*theta[,"rho"]

# Treatment of Side Effects
gamma<-theta[,"gamma"]
gamma2<-theta[,"gamma2"]
recover.amb<--log(1-theta[,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[,"lambda.2.3.fix"])

## Generate data from proposed future study
Q<-50
Size.Tot<-500000

n.burnin <- 200  # Number of burn in iterations
n.iter <- ceiling(Size.Tot/(Q*n.chains)) + n.burnin # Number of iterations per chain
Size.Mat<-ceiling(Size.Tot/Q)

# Find quantiles for the parameters of interest.
pi1.q<-sample(quantile(pi1,probs=1:Q/(Q+1)))
pi2.q<-sample(quantile(pi2,probs=1:Q/(Q+1)))
gamma.q<-sample(quantile(gamma,probs=1:Q/(Q+1)))
gamma2.q<-sample(quantile(gamma2,probs=1:Q/(Q+1)))
recover.amb.q<-sample(quantile(recover.amb,probs=1:Q/(Q+1)))
recover.hosp.q<-sample(quantile(recover.hosp,probs=1:Q/(Q+1)))

var.X.mat<-matrix(NA,nrow=Q,ncol=1)
for(i in 1:Q){
  # Generate data with 150 patients
  n<-150
  X.SE1<-rbinom(1,n,pi1.q[i])
  X.SE2<-rbinom(1,n,pi2.q[i])
  
  X.N.hosp<-rbinom(1,X.SE1+X.SE2,gamma.q[i])
  X.N.die<-rbinom(1,X.N.hosp,gamma2.q[i])
  
  N.amb<-X.SE1+X.SE2-X.N.hosp
  T.re.amb<-rexp(N.amb,recover.amb.q[i])
  N.hosp<-X.N.hosp-X.N.die
  T.re.hosp<-rexp(N.hosp,recover.hosp.q[i])
  
  data.dat <- list("X.SE1","X.SE2","n",
                   "X.N.hosp","X.N.die",
                   "N.amb","T.re.amb","N.hosp","T.re.hosp")
  
  # Feed data to Bayesian model
  data.full<-append(data,data.dat)
  
  # Perform the MCMC simulation with JAGS.
  jags.data <- jags(
    data =  data.full,
    inits = inits,
    parameters.to.save = parameters.to.save,
    model.file = filein, 
    n.chains = n.chains, 
    n.iter = n.iter, 
    n.thin = 1, 
    n.burnin = n.burnin,progress.bar = "none") 
  
  # Run health economic model
  e<-matrix(NA,ncol=2,nrow=Size.Mat)
  c<-matrix(NA,ncol=2,nrow=Size.Mat)
  for(l in 1:Size.Mat){
    e[l,]<-effects(jags.data$BUGSoutput$sims.list[["pi1"]][l],
                   jags.data$BUGSoutput$sims.list[["pi2"]][l],
                   jags.data$BUGSoutput$sims.list[["SE"]][l,1],
                   jags.data$BUGSoutput$sims.list[["SE"]][l,2],
                   jags.data$BUGSoutput$sims.list[["lambda.1.1"]][l],
                   jags.data$BUGSoutput$sims.list[["lambda.1.2"]][l],
                   jags.data$BUGSoutput$sims.list[["lambda.1.3"]][l],
                   jags.data$BUGSoutput$sims.list[["lambda.2.2"]][l],
                   jags.data$BUGSoutput$sims.list[["lambda.2.3"]][l],
                   jags.data$BUGSoutput$sims.list[["lambda.2.4"]][l],
                   jags.data$BUGSoutput$sims.list[["e.chemo"]][l],
                   jags.data$BUGSoutput$sims.list[["e.amb"]][l],
                   jags.data$BUGSoutput$sims.list[["e.hosp"]][l],
                   N,TH)
    c[l,]<-costs(jags.data$BUGSoutput$sims.list[["pi1"]][l],
                 jags.data$BUGSoutput$sims.list[["pi2"]][l],
                 jags.data$BUGSoutput$sims.list[["SE"]][l,1],
                 jags.data$BUGSoutput$sims.list[["SE"]][l,2],
                 jags.data$BUGSoutput$sims.list[["lambda.1.1"]][l],
                 jags.data$BUGSoutput$sims.list[["lambda.1.2"]][l],
                 jags.data$BUGSoutput$sims.list[["lambda.1.3"]][l],
                 jags.data$BUGSoutput$sims.list[["lambda.2.2"]][l],
                 jags.data$BUGSoutput$sims.list[["lambda.2.3"]][l],
                 jags.data$BUGSoutput$sims.list[["lambda.2.4"]][l],
                 c.drug,
                 jags.data$BUGSoutput$sims.list[["c.amb"]][l],
                 jags.data$BUGSoutput$sims.list[["c.hosp"]][l],
                 jags.data$BUGSoutput$sims.list[["c.death"]][l],
                 N,TH)
  }
  
  # Calculate the net monetary benefit
  NB.X<-wtp*e-c
  # Calculate the variance of incremental net benefit
  var.X.mat[i]<-var(NB.X[,2]-NB.X[,1])
}

# Calculate EVSI by rescaling the EVPPI simulations.
prepost.MM<-(save-mean(save))/sd(save)*sqrt(var(INB)-mean(var.X.mat))+mean(save)
EVSI.calc<-mean(pmax(prepost.MM,0))-max(0,mean(prepost.MM))

### EVSI
EVSI.calc


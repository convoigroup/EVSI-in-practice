set.seed(1234) # It is a good practice to set the seed so results can be replicated

# The following R packages are needed.  
library(BCEA)        # Bayesian Cost-Effectiveness Analysis
library(R2OpenBUGS)  # Interfaces R with OpenBUGS
library(ggplot2)                        # A plotting package
library(reshape)                        # Data manipulation package
library(plyr)                           # Data manipulation package
library(R2jags)
library(mgcv)

library(foreach)
library(doParallel)
source('predict_ga.R', encoding = 'WINDOWS-1252')

model<-function(){
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

lognPar <- function(m,s) {
  s2 <- s^2
  mulog <- log(m) - 0.5 * log(1+s2/m^2)
  s2log <- log(1+(s2/m^2))
  sigmalog <- sqrt(s2log)
  list(mulog = mulog, sigmalog = sigmalog)
}

betaPar <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

############################################################
# DATA FOR THE MODEL
############################################################

# Safety Data
num.pat <- 111 # Number of patients in observed data
num.se <- 27   # Number of patients with side effects, given standard-of-care
num.amb <- 17  # Number of patient with hospital care following side effect, given stardard-or-care
num.death <- 1 #Number of Deaths

#Priors recovery
m.1.3<-0.45
v.1.3<-0.02
p1.1.3<-betaPar(m.1.3,v.1.3)$alpha
p2.1.3<-betaPar(m.1.3,v.1.3)$beta

m.2.3<-0.35
v.2.3<-0.02
p1.2.3<-betaPar(m.2.3,v.2.3)$alpha
p2.2.3<-betaPar(m.2.3,v.2.3)$beta

# Probability reduction (=1-rho) of side effect given new treatment 
# pi[2] <- rho * pi[1]
m.rho <- 0.65 # Mean 
s.rho <- 0.1 # Standard deviation
tau.rho <- 1/s.rho^2 # Precision

# Costs
# Cost of ambulatory care 
mu.amb <- 2300 # Mean
sd.amb <- 90  # Standard deviation
m.amb <- lognPar(mu.amb,sd.amb)$mulog      # Mean of log normal distribution
s.amb <- lognPar(mu.amb,sd.amb)$sigmalog   # SD of log normal distribution
tau.amb <- 1/s.amb^2

# Cost of hospitalization
mu.hosp <- 6500
sd.hosp <- 980
m.hosp <- lognPar(mu.hosp,sd.hosp)$mulog
s.hosp <- lognPar(mu.hosp,sd.hosp)$sigmalog
tau.hosp <- 1/s.hosp^2

#Cost of Death
mu.death <- 4200
sd.death <- 560
m.death <- lognPar(mu.death,sd.death)$mulog
s.death <- lognPar(mu.death,sd.death)$sigmalog
tau.death <- 1/s.death^2

# Drug costs
c.drug <- c(110, 420)

# Effects
# QALY on Chemo
mu.chemo<- 0.98
var.chemo<- 0.001
p1.chemo<-betaPar(mu.chemo,var.chemo)$alpha
p2.chemo<-betaPar(mu.chemo,var.chemo)$beta

# QALY on Amb
mu.e.amb<- 0.5
var.e.amb<- 0.02
p1.amb<-betaPar(mu.e.amb,var.e.amb)$alpha
p2.amb<-betaPar(mu.e.amb,var.e.amb)$beta

# QALY on hosp
mu.e.hosp<- 0.2
var.e.hosp<- 0.03
p1.hosp<-betaPar(mu.e.hosp,var.e.hosp)$alpha
p2.hosp<-betaPar(mu.e.hosp,var.e.hosp)$beta

# Number of patients in the population to model
N <- 1000

# Time Horizon on Side Effects
TH<-15


############################################################
# RUN THE BUGS MODEL - PRIOR
############################################################
Size.Prior<-100000

# Data needed for model. These match with variable names above
data <- list("num.pat", "num.se",
             "num.amb", "num.death",
             "p1.1.3","p2.1.3",
             "p1.2.3","p2.2.3",
             "m.rho", "tau.rho",
             "m.amb", "tau.amb",
             "m.hosp", "tau.hosp",
             "m.death", "tau.death",
             "N",
             "p1.chemo","p2.chemo",
             "p1.amb","p2.amb",
             "p1.hosp","p2.hosp",
             "TH")

# The initial values generator function 
inits <- function(){
  list(     pi = c(runif(1),NA),
            gamma = runif(1)#,
            # c.hosp = rnorm(1,Exp,5)
  )
}

n.chains <- 3     # Number of chains
n.burnin <- 1000  # Number of burn in iterations
n.iter <- ceiling(Size.Prior/n.chains) + n.burnin # Number of iterations per chain


# Choose the parameters in the model to monitor
parameters.to.save <- c("pi1", "pi2", "rho", "gamma","gamma2", "SE",
                        "lambda.1.1","lambda.1.2","lambda.1.3","lambda.2.2","lambda.2.3","lambda.2.4",
                        "lambda.2.3.fix","lambda.1.3.fix",
                        "c.amb", "c.hosp","c.death",
                        "e.chemo","e.amb","e.hosp")

filein <- file.path(tempdir(),fileext="psitemp.txt")
write.model(model,filein)

# Perform the MCMC simulation with OpenBUGS.
# Close OpenBUGS once it has finished (if debug is set to TRUE)
prior.model <- jags(
  data =  data,
  inits = inits,
  parameters.to.save = parameters.to.save,
  model.file = filein, 
  n.chains = n.chains, 
  n.iter = n.iter, 
  n.thin = 1, 
  n.burnin = n.burnin) 


####MODEL####
effects<-function(pi1,pi2,SE1,SE2,lambda.1.1,lambda.1.2,lambda.1.3,lambda.2.2,lambda.2.3,lambda.2.4,
                  e.chemo,e.amb,e.hosp,N,TH){
  #QALY of CHEMO
  q.chemo1<-(N-SE1)*e.chemo*(TH+1)
  q.chemo2<-(N-SE2)*e.chemo*(TH+1)
  
  #QALY of SIDE EFFECTS
  MM.mat<-matrix(c(lambda.1.1,lambda.1.2,lambda.1.3,0,
                   0,lambda.2.2,lambda.2.3,lambda.2.4,
                   0,0,1,0,
                   0,0,0,1),nrow=4,ncol=4,byrow=TRUE)
  init1<-c(SE1,0,0,0)
  trace1<-matrix(NA,nrow=4,ncol=TH+1)
  trace1[,1]<-init1
  init2<-c(SE2,0,0,0)
  trace2<-matrix(NA,nrow=4,ncol=TH+1)
  trace2[,1]<-init2
  for(i in 2:(TH+1)){
    trace1[,i]<-trace1[,i-1]%*%MM.mat
    trace2[,i]<-trace2[,i-1]%*%MM.mat
  }
  
  #QALY of STATES
  e.states<-c(e.amb,e.hosp,e.chemo,0)
  
  q.se1<-sum(e.states%*%trace1)
  q.se2<-sum(e.states%*%trace2)
  
  effs<-c(q.se1+q.chemo1,q.se2+q.chemo2)/(N*(TH+1))
  return(effs)
}
e<-matrix(NA,ncol=2,nrow=Size.Prior)
for(l in 1:Size.Prior){
  e[l,]<-effects(prior.model$BUGSoutput$sims.list[["pi1"]][l],
                 prior.model$BUGSoutput$sims.list[["pi2"]][l],
                 prior.model$BUGSoutput$sims.list[["SE"]][l,1],
                 prior.model$BUGSoutput$sims.list[["SE"]][l,2],
                 prior.model$BUGSoutput$sims.list[["lambda.1.1"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.1.2"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.1.3"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.2.2"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.2.3"]][l],
                 prior.model$BUGSoutput$sims.list[["lambda.2.4"]][l],
                 prior.model$BUGSoutput$sims.list[["e.chemo"]][l],
                 prior.model$BUGSoutput$sims.list[["e.amb"]][l],
                 prior.model$BUGSoutput$sims.list[["e.hosp"]][l],
                 N,TH)
}


costs<-function(pi1,pi2,SE1,SE2,lambda.1.1,lambda.1.2,lambda.1.3,lambda.2.2,lambda.2.3,lambda.2.4,
                c.drug,c.amb,c.death,c.hosp,N,TH){
  
  #QALY of SIDE EFFECTS
  MM.mat<-matrix(c(lambda.1.1,lambda.1.2,lambda.1.3,0,
                   0,lambda.2.2,lambda.2.3,lambda.2.4,
                   0,0,1,0,
                   0,0,0,1),nrow=4,ncol=4,byrow=TRUE)
  init1<-c(SE1,0,0,0)
  trace1<-matrix(NA,nrow=4,ncol=TH+1)
  trace1[,1]<-init1
  for(i in 2:(TH+1)){
    trace1[,i]<-trace1[,i-1]%*%MM.mat
  }
  
  init2<-c(SE2,0,0,0)
  trace2<-matrix(NA,nrow=4,ncol=TH+1)
  trace2[,1]<-init2
  for(i in 2:(TH+1)){
    trace2[,i]<-trace2[,i-1]%*%MM.mat
  }
  
  #QALY of STATES
  c.states<-c(c.amb,c.hosp,0,c.death)
  
  c.se1<-sum(c.states%*%trace1)/(N*(TH+1))
  c.se2<-sum(c.states%*%trace2)/(N*(TH+1))
  
  costs<-c(c.drug[1]+c.se1,c.drug[2]+c.se2)
  return(costs)
}

c<-matrix(NA,ncol=2,nrow=Size.Prior)
for(l in 1:Size.Prior){
  c[l,]<-costs(prior.model$BUGSoutput$sims.list[["pi1"]][l],
               prior.model$BUGSoutput$sims.list[["pi2"]][l],
               prior.model$BUGSoutput$sims.list[["SE"]][l,1],
               prior.model$BUGSoutput$sims.list[["SE"]][l,2],
               prior.model$BUGSoutput$sims.list[["lambda.1.1"]][l],
               prior.model$BUGSoutput$sims.list[["lambda.1.2"]][l],
               prior.model$BUGSoutput$sims.list[["lambda.1.3"]][l],
               prior.model$BUGSoutput$sims.list[["lambda.2.2"]][l],
               prior.model$BUGSoutput$sims.list[["lambda.2.3"]][l],
               prior.model$BUGSoutput$sims.list[["lambda.2.4"]][l],
               c.drug,
               prior.model$BUGSoutput$sims.list[["c.amb"]][l],
               prior.model$BUGSoutput$sims.list[["c.hosp"]][l],
               prior.model$BUGSoutput$sims.list[["c.death"]][l],
               N,TH)
}

# Performs the baseline cost-effectiveness analysis for a specified willingness to pay.
wtp<-30000
NB<-e*wtp-c
rm(e,c)

# Matrix of parameters of interest from baseline model.
extra.lines<-(Size.Prior+1):dim(prior.model$BUGSoutput$sims.matrix)[1]
theta<-as.data.frame(prior.model$BUGSoutput$sims.matrix[-extra.lines,c("pi1","rho","gamma","gamma2","lambda.2.3.fix","lambda.1.3.fix","SE[1]","SE[2]")])
colnames(theta)<-c("pi1","rho","gamma","gamma2","lambda.2.3.fix","lambda.1.3.fix","SE1","SE2")
rm(prior.model,extra.lines,l)

## Find the incremental net benefit for each treatment
INB<-NB[,2]-NB[,1]
rm(NB)
### EVPPI estimation
lmm1<-gam(INB~s(pi1)+s(rho)+s(gamma)+s(gamma2)+s(lambda.2.3.fix)+s(lambda.1.3.fix)+s(SE1)+s(SE2),data=theta)

#Model
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

#####################Jalal and Alarid-Escudero#####################
###Estimating n0
Size.Outer<-1000
Size.Inner<-10000

# Set the number of chains in iterations
n.chains <- 3     # Number of chains
n.burnin <- 1000  # Number of burn in iterations
n.iter <- ceiling(Size.Inner/n.chains) + n.burnin # Number of iterations per chain

rm(m.1.3,m.2.3,mu.amb,mu.chemo,mu.death,mu.e.hosp,mu.hosp,s.amb,s.death,
   s.hosp.s.rho,sd.amb,sd.hosp,Size.Prior,Size.Tot,v.1.3,v.2.3,var.chemo,var.e.amb,var.e.hosp,
   betaPar,model,lognPar,sd.death,s.rho,s.hosp,m.e.amb)

pi2<-theta[,"pi1"]*theta[,"rho"]
recover.amb<--log(1-theta[,"lambda.1.3.fix"])
recover.hosp<--log(1-theta[,"lambda.2.3.fix"])

no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)
uncert<-200
EVSI.Jalal.uncert<-foreach(i=1:(uncert),.combine=c,
                           .export=c("Size.Outer","Size.Inner",
                                     "n.chains","n.burnin","n.iter",
                                     "num.pat", "num.se",
                                     "num.amb", "num.death",
                                     "p1.1.3","p2.1.3",
                                     "p1.2.3","p2.2.3",
                                     "m.rho", "tau.rho",
                                     "m.amb", "tau.amb",
                                     "m.hosp", "tau.hosp",
                                     "m.death", "tau.death",
                                     "N",
                                     "p1.chemo","p2.chemo",
                                     "p1.amb","p2.amb",
                                     "p1.hosp","p2.hosp",
                                     "TH","predict.ga","lmm1",
                                     "filein","inits","parameters.to.save"),
                           .packages=c("R2jags","R2OpenBUGS","mgcv")) %dopar% {
                             pi1.mean<-rho.mean<-gamma.mean<-gamma2.mean<-lambda.2.3.fix.mean<-lambda.1.3.fix.mean<-SE1.mean<-SE2.mean<-array(NA,dim=Size.Outer)
                             for(i in 1:Size.Outer){
                               #Generate Data
                               #Change sample size as the same n0 can be used for each sample size and it is faster to run 
                               #an economic model with n<-30
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
                             rm(bugs.data)
                             
                             
                             ##Estimating n0 separately for each element of the meta-model
                             n<-30
                             var.theta<-apply(theta,2,var)
                             pi1.n0<-n*(var.theta[1]/var(pi1.mean)-1)
                             rho.n0<-n*(var.theta[2]/var(rho.mean)-1)
                             gamma.n0<-n*(var.theta[3]/var(gamma.mean)-1)
                             gamma2.n0<-n*(var.theta[4]/var(gamma2.mean)-1)
                             lambda2.n0<-n*(var.theta[5]/var(lambda.2.3.fix.mean)-1)
                             lambda1.n0<-n*(var.theta[6]/var(lambda.1.3.fix.mean)-1)
                             SE1.n0<-n*(var.theta[7]/var(SE1.mean)-1)
                             SE2.n0<-n*(var.theta[8]/var(SE2.mean)-1)
                             
                             rm(pi1.mean,rho.mean,gamma.mean,gamma2.mean,lambda.2.3.fix.mean,lambda.1.3.fix.mean,
                                SE.mean,SE1.mean,SE2.mean)
                             n0<-c(pi1.n0,rho.n0,gamma.n0,gamma2.n0,lambda2.n0,lambda1.n0,SE1.n0,SE2.n0)
                             n<-rep(150,8)
                             llpred<-predict.ga(lmm1,n0=n0,n=n)
                             evsi.Jal <- mean(pmax(0,llpred))-max(mean(llpred),0)
                             rm(llpred)
                             evsi.Jal
                           }
stopCluster(cl)

write.table(EVSI.Jalal.uncert,"/home/aheath/Chemotherapy/EVSIJalal.txt")
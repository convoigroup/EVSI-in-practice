#Required packages
library(earth)
# Performs the baseline cost-effectiveness analysis
# Estimates EVPPI using the Heath et al method - object names evppi.half
source("C:/Users/anna heath/Google Drive/PhD Files/EVSI/Method Comparison/Chemotherapy Model.R")
# Loads function used for the Jalal methods
source('C:/Users/anna heath/Google Drive/PhD Files/EVSI/Practical Paper/predict_ga.R', encoding = 'WINDOWS-1252')

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
lmm1<-gam(y~s(pi1)+s(rho)+s(gamma)+s(gamma2)+s(lambda.2.3.fix)+s(lambda.1.3.fix)+s(SE1)+s(SE2),data=theta)

evppi <- mean(pmax(0,lmm1$fitted.values))
evppi


###Heath et al#########################
Size.Tot<-500000
Size.Dist<-50
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

#Number of Side Effects
pi1<-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"pi1"]
pi2<-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"pi2"]

#Treatment of Side Effects
gamma<-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"gamma"]
gamma2<-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"gamma2"]
recover.amb<--log(1-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"lambda.1.3.fix"])
recover.hosp<--log(1-as.data.frame(prior.model$BUGSoutput$sims.matrix)[,"lambda.2.3.fix"])
##Across Q Values
Q<-50
var.X.mat<-array(NA,dim=Q)
n.burnin <- 200  # Number of burn in iterations
n.iter <- ceiling(Size.Tot/(Q*n.chains)) + n.burnin # Number of iterations per chain
Size.Mat<-ceiling(Size.Tot/Q)
#Quantiles
pi1.q<-sample(quantile(pi1,probs=1:Q/(Q+1)))
pi2.q<-sample(quantile(pi2,probs=1:Q/(Q+1)))
gamma.q<-sample(quantile(gamma,probs=1:Q/(Q+1)))
gamma2.q<-sample(quantile(gamma2,probs=1:Q/(Q+1)))
recover.amb.q<-sample(quantile(recover.amb,probs=1:Q/(Q+1)))
recover.hosp.q<-sample(quantile(recover.hosp,probs=1:Q/(Q+1)))
for(i in 1:Q){
  ##Moment Matching
  #Generate Data
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
    n.burnin = n.burnin,progress.bar = "none") 
      
  e<-matrix(NA,ncol=2,nrow=Size.Mat)
  c<-matrix(NA,ncol=2,nrow=Size.Mat)
  for(l in 1:Size.Mat){
    e[l,]<-effects(bugs.data$BUGSoutput$sims.list[["pi1"]][l],
                   bugs.data$BUGSoutput$sims.list[["pi2"]][l],
                   bugs.data$BUGSoutput$sims.list[["SE"]][l,1],
                   bugs.data$BUGSoutput$sims.list[["SE"]][l,2],
                   bugs.data$BUGSoutput$sims.list[["lambda.1.1"]][l],
                   bugs.data$BUGSoutput$sims.list[["lambda.1.2"]][l],
                   bugs.data$BUGSoutput$sims.list[["lambda.1.3"]][l],
                   bugs.data$BUGSoutput$sims.list[["lambda.2.2"]][l],
                   bugs.data$BUGSoutput$sims.list[["lambda.2.3"]][l],
                   bugs.data$BUGSoutput$sims.list[["lambda.2.4"]][l],
                   bugs.data$BUGSoutput$sims.list[["e.chemo"]][l],
                   bugs.data$BUGSoutput$sims.list[["e.amb"]][l],
                   bugs.data$BUGSoutput$sims.list[["e.hosp"]][l],
                   N,TH)
    c[l,]<-costs(bugs.data$BUGSoutput$sims.list[["pi1"]][l],
                 bugs.data$BUGSoutput$sims.list[["pi2"]][l],
                 bugs.data$BUGSoutput$sims.list[["SE"]][l,1],
                 bugs.data$BUGSoutput$sims.list[["SE"]][l,2],
                 bugs.data$BUGSoutput$sims.list[["lambda.1.1"]][l],
                 bugs.data$BUGSoutput$sims.list[["lambda.1.2"]][l],
                 bugs.data$BUGSoutput$sims.list[["lambda.1.3"]][l],
                 bugs.data$BUGSoutput$sims.list[["lambda.2.2"]][l],
                 bugs.data$BUGSoutput$sims.list[["lambda.2.3"]][l],
                 bugs.data$BUGSoutput$sims.list[["lambda.2.4"]][l],
                 c.drug,
                 bugs.data$BUGSoutput$sims.list[["c.amb"]][l],
                 bugs.data$BUGSoutput$sims.list[["c.hosp"]][l],
                 bugs.data$BUGSoutput$sims.list[["c.death"]][l],
                 N,TH)
  }
  NB.X<-wtp*e-c
      
  var.X.mat[i]<-var(NB.X[,2]-NB.X[,1])
}


prepost.MM<-(save-mean(save))/sd(save)*sqrt(var(NB[,2]-NB[,1])-mean(var.X.mat))+mean(save)
EVSI.heath<-mean(pmax(prepost.MM,0))-max(0,mean(prepost.MM))


#####################Jalal and Alarid-Escudero#####################
###Estimating n0
Size.Outer<-1000
Size.Inner<-10000

# Set the number of chains in iterations
n.chains <- 3     # Number of chains
n.burnin <- 1000  # Number of burn in iterations
n.iter <- ceiling(Size.Inner/n.chains) + n.burnin # Number of iterations per chain

pi1.mean<-rho.mean<-gamma.mean<-gamma2.mean<-lambda.2.3.fix.mean<-lambda.1.3.fix.mean<-SE1.mean<-SE2.mean<-array(NA,dim=Size.Outer)
#pred<-predict.gam(lmm1,type="lpmatrix")
#ind.means<-array(NA,dim=c(Size.Outer,dim(pred)[2]))

start<-Sys.time()
for(i in 1:Size.Outer){
  #Generate Data
  #Change sample size as the same n0 can be used for each sample size and it is faster to run 
  #an economic model with n<-30
  n<-30
  X.SE1<-rbinom(1,n,pi1[i])
  X.SE2<-rbinom(1,n,pi2[i])
  
  X.N.hosp<-rbinom(1,X.SE1+X.SE2,gamma[i])
  X.N.die<-rbinom(1,X.N.hosp,gamma2[i])
  
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
end<-Sys.time()

time.Jalal<-end-start

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

n0<-c(pi1.n0,rho.n0,gamma.n0,gamma2.n0,lambda2.n0,lambda1.n0,SE1.n0,SE2.n0)
n<-rep(150,8)
llpred<-predict.ga(lmm1,n0=n0,n=n)
evsi.Jal <- mean(pmax(0,llpred))-max(mean(llpred),0)

###############Strong et al.#########################

start<-Sys.time()
X.SE1<-X.SE2<-X.N.die<-X.N.hosp<-X.amb<-X.hosp<-array(dim=Size.Prior)
for(i in 1:Size.Prior){
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

dat<-as.data.frame(cbind(y,X.amb,X.SE1,X.SE2,X.N.hosp,X.hosp,X.N.die))
prepost.s<-gam(y~s(X.amb)+s(X.SE1)+s(X.SE2)+s(X.N.hosp)+s(X.hosp)+s(X.N.die),data=dat)
end<-Sys.time()

time.Strong<-difftime(end,start,units="auto")

EVSI.strong<-mean(pmax(0,prepost.s$fitted.values))-max(0,mean(prepost.s$fitted.values))

##############Menzies########################
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

Size.Outer<-20000
params<-cbind(pi1,pi2,gamma,gamma2,recover.amb,recover.hosp)[1:Size.Outer,]
prepost.M<-array(dim=Size.Outer)
start<-Sys.time()
for(i in 1:Size.Outer){
  n<-150
  X.SE1<-rbinom(1,n,pi1[i])
  X.SE2<-rbinom(1,n,pi2[i])
  X.N.hosp<-rbinom(1,X.SE1+X.SE2,gamma[i])
  X.N.die<-rbinom(1,X.N.hosp,gamma2[i])
  N.amb<-X.SE1+X.SE2-X.N.hosp
  T.re.amb<-rexp(N.amb,recover.amb[i])
  X.amb<-sum(T.re.amb)
  N.hosp<-X.N.hosp-X.N.die
  T.re.hosp<-rexp(N.hosp,recover.hosp[i])
  X.hosp<-sum(T.re.hosp)
  D<-as.data.frame(cbind(X.SE1,X.SE2,X.amb,X.hosp,X.N.die,X.N.hosp))  
  ll<-likelihood(D,params)
  w<-ll/sum(ll)
  prepost.M[i]<-w%*%y[1:Size.Outer]
}
end<-Sys.time()

time.Menzies<-end-start
evsi.Menzies<-mean(pmax(0,prepost.M))-max(0,mean(prepost.M))

cbind(SeA=EVSI.strong,
      M=evsi.Menzies,
      JA=evsi.Jal,
      HeA=mean(pmax(prepost.MM,0),na.rm=TRUE)-max(0,mean(prepost.MM,na.rm=TRUE)),
      Truth = evsi.true)


cbind(SeA=time.Strong,
      M=time.Menzies,
      JA=time.Jalal,
      HeA=NA,
      Truth = "42 Days")

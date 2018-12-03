#Markov Model for Chronic Pain - Adapted from Sullivan et al.
#Copyright: Anna Heath 2017
#Model has two treatment options for chronic pain and costs and effects are determined using a Markov Model
#Treatment 1: Morphine
#Treatment 2: Novel Treatment
#N: PSA sample size
#NOTE: Throughout the PSA distributions are taken as having the mean given by the parameter estimate and the standard error as 10%
#of the parameter estimate.
#t is treatment and l is the number of treatments that have been attempted

##Packages
library(BCEA)
library(R2OpenBUGS)
library(R2jags)
library(mgcv)
library(psych)

gamma.par<-function(par.est){
  alpha=100
  beta=100/par.est
  scale=par.est/100
  return(list(alpha=alpha,beta=beta,scale=scale))
}
beta.par<-function(par.est){
  alpha<-(1-par.est)/0.01-par.est
  beta<-alpha*(1/par.est-1)
  return(list(alpha=alpha,beta=beta))
}
prob.novel<-0.3
cost.novel<-6

N<-100000
#Costs
#Treatment costs are considered known and taken from the literature. The novel treatment is assumed to be
#6 times more expensive than Oxycodone
c.t1<-2.632729166666670000
c.t2<-9.2011500000*cost.novel

#The comedication cost per cycle - this is for complications associated with the pain medication
#Costs based on a previous study so inflation needs to be taken into account
PriceIndex0910<-268.6
PriceIndex1213<-289.1
Inflation<-PriceIndex1213/PriceIndex0910

c.med.t1<-rgamma(N,shape=gamma.par(2.1*Inflation)$alpha,rate=gamma.par(2.1*Inflation)$beta)
#Novel theraphy gives an improvement on Oxycodone so the cost is based on this improvement
c.med.t2<-rgamma(N,shape=gamma.par(0.04*Inflation*(1-prob.novel))$alpha,rate=gamma.par(0.04*Inflation*(1-prob.novel))$beta)

#Costs of adverse events
c.ae<-rgamma(N,shape=gamma.par(6.991009409)$alpha,rate=gamma.par(6.991009409)$beta)

#The cost of withdrawing from the theraphy is the same irrespective of the reason you withdaw
c.withdraw.ae<-rgamma(N,shape=gamma.par(106.911031273)$alpha,rate=gamma.par(106.911031273)$beta)
c.withdraw<-rgamma(N,shape=gamma.par(106.911031273)$alpha,rate=gamma.par(106.911031273)$beta)

#Cost of discontinuing treatment is based on visiting the GP
c.dist<-rgamma(N,shape=gamma.par(18.5)$alpha,rate=gamma.par(18.5)$beta)

#Cost of second round of treatment
c.l2<-9.2011500000+0.04*Inflation

#Cost of third round of treatment
c.l3<-2.632729166666670000+2.1*Inflation

####Utilities###
#No adverse effects for the first treatment
u.l1.noae<-rbeta(N,beta.par(0.695000000)$alpha,beta.par(0.695000000)$beta)
#Adverse effects
u.l1.ae<-rbeta(N,beta.par(0.583000000)$alpha,beta.par(0.583000000)$beta)
#Withdraw from treatment due to adverse effects
u.l1.withdraw.ae<-rbeta(N,beta.par(0.503000000)$alpha,beta.par(0.503000000)$beta)
#Withdraw due to other reasons
u.l1.withdraw.noae<-rbeta(N,beta.par(0.405000000)$alpha,beta.par(0.405000000)$beta)

#Multiplier to give the utilities when on the 2nd treatment options
u.l2<-0.900

#Utility for 3rd case of treatment
u.l3<-(u.l1.noae+u.l1.ae)/2

#Discontinuing treatment
u.dist<-u.l1.withdraw.noae*0.8

####Transition Probabilities####
#For the first round of treatments
#probability of adverse effects
p.ae.l1.t1<-rbeta(N,beta.par(0.436159243)$alpha,beta.par(0.436159243)$beta)
p.ae.l1.t2<-p.ae.l1.t1*(1-prob.novel)
#probility of withdrawal due to adverse effects
p.with.ae.l1.t1<-rbeta(N,beta.par(0.055739588)$alpha,beta.par(0.055739588)$beta)
p.with.ae.l1.t2<-rbeta(N,beta.par(0.022958454)$alpha,beta.par(0.022958454)$beta)
#probability of withdrawal due to other reasons
p.with.l1.t1<-rbeta(N,beta.par(0.012741455)$alpha,beta.par(0.012741455)$beta)
p.with.l1.t2<-rbeta(N,beta.par(0.001612408)$alpha,beta.par(0.001612408)$beta)
#probability of discontinuation
p.dist.l1<-rbeta(N,beta.par(0.050000000)$alpha,beta.par(0.050000000)$beta)

#For the second round of treatment that only has one treatment

#probability of adverse effects
p.ae.l2.t1<-rbeta(N,beta.par(0.463500000)$alpha,beta.par(0.463500000)$beta)
p.ae.l2.t2<-rbeta(N,beta.par(0.463500000)$alpha,beta.par(0.463500000)$beta)
#probility of withdrawal due to adverse effects
p.with.ae.l2.t1<-rbeta(N,beta.par(0.032797792)$alpha,beta.par(0.032797792)$beta)
p.with.ae.l2.t2<-rbeta(N,beta.par(0.032797792)$alpha,beta.par(0.032797792)$beta)
#probability of withdrawal due to other reasons
p.with.l2.t1<-rbeta(N,beta.par(0.002303439)$alpha,beta.par(0.002303439)$beta)
p.with.l2.t2<-p.with.l2.t1#rbeta(N,beta.par(0.002303439)$alpha,beta.par(0.002303439)$beta)
#probability of discontinuation
p.dist.l2<-rbeta(N,beta.par(0.100000000)$alpha,beta.par(0.100000000)$beta)

###Transition Matrices###
#First line of treatment l1
#For treatment 1 t1
No.AE.l1.t1<-cbind((1-p.with.ae.l1.t1-p.with.l1.t1)*(1-p.ae.l1.t1),
                (1-p.with.ae.l1.t1-p.with.l1.t1)*(p.ae.l1.t1),
                p.with.ae.l1.t1,
                p.with.l1.t1,
                rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N))

AE.l1.t1<-cbind((1-p.with.ae.l1.t1-p.with.l1.t1)*(1-p.ae.l1.t1),
                (1-p.with.ae.l1.t1-p.with.l1.t1)*(p.ae.l1.t1),
                p.with.ae.l1.t1,
                p.with.l1.t1,
                rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N))

With.AE.l1.t1<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                  (1-p.dist.l1)*(1-p.ae.l2.t1),
                  (1-p.dist.l1)*(p.ae.l2.t1),
                  rep(0,N),rep(0,N),rep(0,N),
                  p.dist.l1)
With.l1.t1<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
               (1-p.dist.l1)*(1-p.ae.l2.t1),
               (1-p.dist.l1)*(p.ae.l2.t1),
               rep(0,N),rep(0,N),rep(0,N),
               p.dist.l1)

#Second line of treatment l2
No.AE.l2.t1<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                (1-p.with.ae.l2.t1-p.with.l2.t1)*(1-p.ae.l2.t1),
                (1-p.with.ae.l2.t1-p.with.l2.t1)*(p.ae.l2.t1),
                p.with.ae.l2.t1,p.with.l2.t1,rep(0,N),rep(0,N))

AE.l2.t1<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
             (1-p.with.ae.l2.t1-p.with.l2.t1)*(1-p.ae.l2.t1),
             (1-p.with.ae.l2.t1-p.with.l2.t1)*(p.ae.l2.t1),
             p.with.ae.l2.t1,p.with.l2.t1,rep(0,N),rep(0,N))

#First line of treatment l1
#For treatment 2 t2
No.AE.l1.t2<-cbind((1-p.with.ae.l1.t2-p.with.l1.t2)*(1-p.ae.l1.t2),
                   (1-p.with.ae.l1.t2-p.with.l1.t2)*(p.ae.l1.t2),
                   p.with.ae.l1.t2,
                   p.with.l1.t2,rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                   rep(0,N),rep(0,N))

AE.l1.t2<-cbind((1-p.with.ae.l1.t2-p.with.l1.t2)*(1-p.ae.l1.t2),
                (1-p.with.ae.l1.t2-p.with.l1.t2)*(p.ae.l1.t2),
                p.with.ae.l1.t2,
                p.with.l1.t2,rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                rep(0,N),rep(0,N))

With.AE.l1.t2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                     (1-p.dist.l1)*(1-p.ae.l2.t2),
                     (1-p.dist.l1)*(p.ae.l2.t2),
                     rep(0,N),rep(0,N),rep(0,N),
                     p.dist.l1)
With.l1.t2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                  (1-p.dist.l1)*(1-p.ae.l2.t2),
                  (1-p.dist.l1)*(p.ae.l2.t2),
                  rep(0,N),rep(0,N),rep(0,N),
                  p.dist.l1)

#Second line of treatment l2
No.AE.l2.t2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                   (1-p.with.ae.l2.t2-p.with.l2.t2)*(1-p.ae.l2.t2),
                   (1-p.with.ae.l2.t2-p.with.l2.t2)*(p.ae.l2.t2),
                   p.with.ae.l2.t2,p.with.l2.t2,rep(0,N),rep(0,N))

AE.l2.t2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),
                (1-p.with.ae.l2.t2-p.with.l2.t2)*(1-p.ae.l2.t2),
                (1-p.with.ae.l2.t2-p.with.l2.t2)*(p.ae.l2.t2),
                p.with.ae.l2.t2,p.with.l2.t2,rep(0,N),rep(0,N))

With.AE.l2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(1,N))
                  #1-p.dist.l2,p.dist.l2)
With.l2<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(1,N))
               #1-p.dist.l2,p.dist.l2)

##Absorbing states
Subs.treat<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(1,N),rep(0,N))
Dist<-cbind(rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(0,N),rep(1,N))



PSA.Trans.Mat.t1<-array(NA,dim=c(10,10,N))
PSA.Trans.Mat.t1[1,,]<-(t(No.AE.l1.t1))
PSA.Trans.Mat.t1[2,,]<-(t(AE.l1.t1))
PSA.Trans.Mat.t1[3,,]<-(t(With.AE.l1.t1))
PSA.Trans.Mat.t1[4,,]<-(t(With.l1.t1))
PSA.Trans.Mat.t1[5,,]<-(t(No.AE.l2.t1))
PSA.Trans.Mat.t1[6,,]<-(t(AE.l2.t1))
PSA.Trans.Mat.t1[7,,]<-(t(With.AE.l2))
PSA.Trans.Mat.t1[8,,]<-(t(With.l2))
PSA.Trans.Mat.t1[9,,]<-(t(Subs.treat))
PSA.Trans.Mat.t1[10,,]<-(t(Dist))

PSA.Trans.Mat.t2<-PSA.Trans.Mat.t1
PSA.Trans.Mat.t2[1,,]<-(t(No.AE.l1.t2))
PSA.Trans.Mat.t2[2,,]<-(t(AE.l1.t2))
PSA.Trans.Mat.t2[3,,]<-(t(With.AE.l1.t2))
PSA.Trans.Mat.t2[4,,]<-(t(With.l1.t2))
PSA.Trans.Mat.t2[5,,]<-(t(No.AE.l2.t2))
PSA.Trans.Mat.t2[6,,]<-(t(AE.l2.t2))

c.Mat.t1<-cbind(c.t1+c.med.t1,
                c.t1+c.med.t1+c.ae,
                c.withdraw.ae,
                c.withdraw,
                c.l2,
                c.l2+c.ae,
                c.withdraw.ae,
                c.withdraw,
                c.l3,
                c.dist)

c.Mat.t2<-cbind(c.t2+c.med.t2,
                c.t2+c.med.t2+c.ae,
                c.withdraw.ae,
                c.withdraw,
                c.l2,
                c.l2+c.ae,
                c.withdraw.ae,
                c.withdraw,
                c.l3,
                c.dist)

u.Mat<-cbind(u.l1.noae,
                u.l1.ae,
                u.l1.withdraw.ae,
                u.l1.withdraw.noae,
                u.l1.noae*u.l2,
                u.l1.ae*u.l2,
                u.l1.withdraw.ae*u.l2,
                u.l1.withdraw.noae*u.l2,
               u.l3,
               u.dist)*7/365.25


Time_Horizen<-52
InitVector<-c(1,0,0,0,0,0,0,0,0,0)

###Markov Model####
Markov_Prob <- function(TransArray){

  trace_matrix <- matrix(NA, nrow=Time_Horizen, ncol=ncol(TransArray))
  trace_matrix[1,] <- InitVector
  trace_matrix[2,] <- InitVector %*% TransArray
  for (i in 3:nrow(trace_matrix)){
    trace_matrix[i,] <- trace_matrix[i-1,] %*% TransArray
  }
  #hcc = half cycle corrected
  hcc_trace_matrix <- matrix(NA, nrow=Time_Horizen, ncol=ncol(TransArray))
  hcc_trace_matrix[1,] <- 0.5*InitVector + 0.5*trace_matrix[2,]
  for (i in 2:Time_Horizen-1){
    hcc_trace_matrix[i,] <- 0.5*trace_matrix[i,] + 0.5*trace_matrix[i+1,]
  }
  hcc_trace_matrix[Time_Horizen,] <- trace_matrix[Time_Horizen,]
  return(hcc_trace_matrix)
}

Prob.Array.1<-array(NA,c(52,10,N))
for(i in 1:N){
  Prob.Array.1[,,i]<-Markov_Prob(PSA.Trans.Mat.t1[,,i])}

Prob.Array.2<-array(NA,c(52,10,N))
for(i in 1:N){
  Prob.Array.2[,,i]<-Markov_Prob(PSA.Trans.Mat.t2[,,i])}

costs.t1<-array(NA,N)
for(i in 1:N){
  costs.t1[i]<-  sum(Prob.Array.1[,,i]%*%c.Mat.t1[i,])}

costs.t2<-array(NA,N)
for(i in 1:N){
  costs.t2[i]<- sum(Prob.Array.2[,,i]%*%c.Mat.t2[i,])
}


effects.t1<-array(NA,N)
for(i in 1:N){
  effects.t1[i]<-  sum(Prob.Array.1[,,i]%*%u.Mat[i,])}

effects.t2<-array(NA,N)
for(i in 1:N){
  effects.t2[i]<- sum(Prob.Array.2[,,i]%*%u.Mat[i,])
}
####EVPPI####
#This gives the standard cost-effectiveness analysis for the Chronic Pain model.
#It also has the EVPPI calculations to determine where to focus analysis.

discount.15<-sum(1/(1+0.035)^(0:15))
m<-bcea(discount.15*cbind(effects.t1,effects.t2),discount.15*cbind(costs.t1,costs.t2),ref=2,wtp=c(0,20000))
var.pr<-var(m$ib[which(m$k==20000),])

pars<-cbind(c.ae,c.dist,c.med.t1,c.med.t2,c.withdraw,c.withdraw.ae,p.ae.l1.t1,p.ae.l1.t2,p.dist.l1,
            p.dist.l2,p.with.ae.l1.t1,p.with.ae.l1.t2,p.with.l1.t1,p.with.l1.t2,p.ae.l2.t1,
            p.ae.l2.t2,p.with.ae.l2.t1,p.with.ae.l2.t2,p.with.l2.t1,p.with.l2.t2,
            u.dist,u.l1.ae,u.l1.noae,u.l1.withdraw.ae,u.l1.withdraw.noae,u.l3)

#Utility
system.time(save.full<-evppi(c("u.l1.noae","u.l1.withdraw.noae"),pars,m,method="gam"))
save<-20000*save.full$fitted.effects-save.full$fitted.costs
save<-save[,1]

betaPar <- function(m,s){
  a <- m*( (m*(1-m)/s^2) -1 )
  b <- (1-m)*( (m*(1-m)/s^2) -1 )
  list(a=a,b=b)
}

###################################################
#There are several approximation methods for the ##
#EVSI - this section calculates the EVSI          #
###################################################
####HEATH et al method#############################
sig.X.noae<-0.300
sig.X.with<-0.310
#Utilities
s.u.l1.noae<-beta.par(0.695000000)$alpha
r.u.l1.noae<-beta.par(0.695000000)$beta
#Withdraw due to other reasons
s.u.l1.withdraw.noae<-beta.par(0.405000000)$alpha
r.u.l1.withdraw.noae<-beta.par(0.405000000)$beta

##Model to update focal parameters
model<-function(){
  for(i in 1:n.model){
    X2[i]~dbeta(s.X2,r.X2)
    X4[i]~dbeta(s.X4,r.X4)
  }
  
  s.X2<- u.l1.noae*( (u.l1.noae*(1-u.l1.noae)/sig.X.noae^2) -1 )
  r.X2<- (1-u.l1.noae)*((u.l1.noae*(1-u.l1.noae)/sig.X.noae^2) -1 )
  s.X4<- u.l1.withdraw.noae*( (u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/sig.X.with^2) -1 )
  r.X4<- (1-u.l1.withdraw.noae)*((u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/sig.X.with^2) -1 )
  
  #Utilities
  u.l1.noae~dbeta(s.u.l1.noae,r.u.l1.noae)
  #Withdraw due to other reasons
  u.l1.withdraw.noae~dbeta(s.u.l1.withdraw.noae,r.u.l1.withdraw.noae)
}
filein <- file.path("~/",fileext="psitemp.txt")
write.model(model,filein)


#Setting up initial conditions 
N<-10000
n.chains<-2
n.burnin <- 1000  # Number of burn in iterations
n.thin<-1
n.iter <- ceiling(N*n.thin/n.chains) # Number of iterations per chain

# Choose the parameters in the model to monitor
parameters.to.save <-cbind("u.l1.noae","u.l1.withdraw.noae")  
datanames<-c("r.u.l1.noae","r.u.l1.withdraw.noae","s.u.l1.noae","s.u.l1.withdraw.noae",
             "X2","X4","n.model","sig.X.noae","sig.X.with")

Q<-50
n<-trunc(seq(10,150,length.out=Q))
phi2<-quantile(u.l1.noae,prob=(1:Q)/(Q+1))
phi4<-quantile(u.l1.withdraw.noae,prob=(1:Q)/(Q+1))
while(abs(cor(n,phi2))>0.001){phi2<-sample(quantile(u.l1.noae,prob=(1:Q)/(Q+1)),replace=F)}
while(abs(cor(n,phi4))>0.001){phi4<-sample(quantile(u.l1.withdraw.noae,prob=(1:Q)/(Q+1)),replace=F)}
discount.15<-sum(1/(1+0.035)^(0:15))
Var.X.prob<-array(NA,Q)
start<-Sys.time()
for(j in 1:Q){
  n.model<-n[j]
  #Data Simulation
  #X1<-rbeta(n.model,betaPar(phi1[j],sig.X)$a,betaPar(phi1[j],sig.X)$b)
  X2<-rbeta(n.model,betaPar(phi2[j],sig.X.noae)$a,betaPar(phi2[j],sig.X.noae)$b)
  #X3<-rbeta(n.model,betaPar(phi3[j],sig.X)$a,betaPar(phi3[j],sig.X)$b)
  X4<-rbeta(n.model,betaPar(phi4[j],sig.X.with)$a,betaPar(phi4[j],sig.X.with)$b)
  
  #Missingness
  missingness<-rbinom(n.model,1,0.687)#0.757
  #X1[which(X1*missingness==0)]<-NA
  X2[which(X2*missingness==0)]<-NA
  #X3[which(X3*missingness==0)]<-NA
  X4[which(X4*missingness==0)]<-NA
  
  
  data<-list(r.u.l1.noae,r.u.l1.withdraw.noae,s.u.l1.noae,s.u.l1.withdraw.noae,
             X2,X4,n.model,sig.X.noae,sig.X.with)
  names(data)<-datanames
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  while(class(tryCatch({
    Model.JAGS<-jags.model(filein, data =  data,
                           inits=list(u.l1.noae=rep(0.8,1),u.l1.withdraw.noae=rep(0.6,1)),
                           n.chains=n.chains,quiet=TRUE)
    update(Model.JAGS,n.iter=n.burnin,progress.bar="none")
    coda.save<-coda.samples(Model.JAGS,parameters.to.save,n.iter=n.iter,thin=n.thin,progress.bar="none")
  },error=function(e) 1))!="mcmc.list"){
    X2<-rbeta(n.model,betaPar(phi2[j],sig.X.noae)$a,betaPar(phi2[j],sig.X.noae)$b)
    X2[which(X2*missingness==0)]<-NA
    data<-list(r.u.l1.noae,r.u.l1.withdraw.noae,s.u.l1.noae,s.u.l1.withdraw.noae,
               X2,X4,n.model,sig.X.noae,sig.X.with)
    names(data)<-datanames
  }
  
  coda.full<-rbind(coda.save[[1]],coda.save[[2]])
  #Utility for 3rd case of treatment
  u.l3<-(coda.full[,1]+u.l1.ae)/2
  
  #Discontinuing treatment
  u.dist<-coda.full[,2]*0.8
  
  u.Mat<-cbind(coda.full[,1],
               u.l1.ae[1:N],
               u.l1.withdraw.ae[1:N],
               coda.full[,2],
               coda.full[,1]*u.l2,
               u.l1.ae[1:N]*u.l2,
               u.l1.withdraw.ae*u.l2,
               coda.full[,2]*u.l2,
               u.l3[1:N],
               u.dist[1:N])*7/365.25
  
  effects.t1.X<-array(NA,N)
  for(i in 1:N){
    effects.t1.X[i]<-sum(Prob.Array.1[,,i]%*%u.Mat[i,])}
  
  effects.t2.X<-array(NA,N)
  for(i in 1:N){
    effects.t2.X[i]<-sum(Prob.Array.2[,,i]%*%u.Mat[i,])}
  
  Var.X.prob[j]<- var(discount.15*(20000*(effects.t2.X-effects.t1.X)-(costs.t2[1:N]-costs.t1[1:N])))
}
end<-Sys.time()
var.calc.time<-difftime(end,start,units="secs")

model.ab<-function(){
  beta~dnorm(Nmax/2,shape.Nmax)%_%T(0,)
  for(i in 1:N){
    y[i]~dnorm(mu[i],tau)
    mu[i]<-var.PI*(x[i]/(x[i]+beta))
  }
  sigma~dt(sigma.mu,sigma.tau,3)%_%T(0,)#dunif(0.01,50)
  tau<-1/sigma^2
}


data<-list(sigma.mu=sd(var(m$ib[which(m$k==20000),])-Var.X.prob)/2,
           sigma.tau=1/(sd(var(m$ib[which(m$k==20000),])-Var.X.prob)),
           N=length(n),
           shape.Nmax=0.0005/max(n),
           var.PI=var(save),
           Nmax=max(n),
           y=as.vector(t(var(m$ib[which(m$k==20000),])-Var.X.prob)),
           x=as.vector(rep(n,1)))

n.chains<-3
n.burnin <- 1000  # Number of burn in iterations
n.thin<-5
n.iter <- ceiling(10000*n.thin/n.chains) + n.burnin # Number of iterations per chain

inits<-function(){list(params=c(0,20))}
# Choose the parameters in the model to monitor
parameters.to.save <- c("beta","sigma","mu")
library(R2jags)
library(R2OpenBUGS)
filein <- file.path("~/",fileext="abmodel.txt")
write.model(model.ab,filein)

# Perform the MCMC simulation with OpenBUGS.
# Close OpenBUGS once it has finished (if debug is set to TRUE)
bugs.a.b<- jags(
  data =  data,
  inits = NULL,
  parameters.to.save = parameters.to.save,
  model.file = filein, 
  n.chains = n.chains, 
  n.iter = n.iter, 
  n.thin = n.thin, 
  n.burnin = n.burnin) 

n.<-c(10,25,50,100,150)
med<-median(bugs.a.b$BUGSoutput$sims.list$beta)


EVSI.heath<-array(NA,dim=5)
var.pred<-var(save)*(n./(n.+med))
for(l in 1:5){
    samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
    EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
}


####GAM Fitting - STRONG et al method####
EVSI.gam<-array()
var.fit<-array()
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)
N<-100000
start<-Sys.time()
for(i in 1:5){
  n.loop<-n.[i]
  X.gam<-array(NA,dim=c(N,2))
  X.gam.with<-array(NA,dim=c(N,2))
  for(j in 1:N){
    rand<-rbinom(n.loop,1,0.687)
    samp<-rbeta(sum(rand),betaPar(u.l1.noae[j],sig.X.noae)$a,betaPar(u.l1.noae[j],sig.X.noae)$b)
    samp.with<-rbeta(sum(rand),betaPar(u.l1.ae[j],sig.X.with)$a,betaPar(u.l1.ae[j],sig.X.with)$b)
    X.gam[j,]<-c(psych::geometric.mean(samp),psych::geometric.mean(1-samp))
    X.gam.with[j,]<-c(psych::geometric.mean(samp.with),psych::geometric.mean(1-samp.with))
    }
gam.fit<-gam(m$ib[which(m$k==20000),1:N]~te(X.gam,X.gam.with))
EVSI.gam[i]<-mean(pmax(0,gam.fit$fitted))-max(0,mean(gam.fit$fitted))
var.fit[i]<-var(gam.fit$fitted)}
end<-Sys.time()

####MENZIES####
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)
N<-5000

sim_data<-function(u.noae,u.ae,n,sig.X.noae=0.300,sig.X.with=0.310){
  rand<-rbinom(n,1,0.687)
  samp<-rbeta(sum(rand),betaPar(u.noae,sig.X.noae)$a,betaPar(u.noae,sig.X.noae)$b)
  samp.with<-rbeta(sum(rand),betaPar(u.ae,sig.X.with)$a,betaPar(u.ae,sig.X.with)$b)
  r<-cbind(samp,samp.with)
  return(r)
}
func.dbeta<-function(u.l1.noae,u.l1.ae,X,sig.X.noae=0.300,sig.X.with=0.310){
  noae.par<-betaPar(u.l1.noae,sig.X.noae)
  ae.par<-betaPar(u.l1.ae,sig.X.with)
  d<-exp(sum(dbeta(X[,1],noae.par$a,noae.par$b,log=T))+
           sum(dbeta(X[,2],ae.par$a,ae.par$b,log=T)))
  return(d)
}

start<-Sys.time()
mu.X<-array(NA,dim=N)
EVSI.men<-array(NA,dim=length(n.))
for(i in 1:length(n.)){
  for(j in 1:N){
dat.M<-sim_data(u.l1.noae[j],u.l1.ae[j],n.[i])
app<-function(u.l1.noae,u.l1.ae){
  r<-func.dbeta(u.l1.noae,u.l1.ae,X=dat.M)
  return(r)
}
likeli<-mapply(app,u.l1.noae=u.l1.noae[1:N],u.l1.ae=u.l1.ae[1:N])
weights<-likeli/sum(likeli)
mu.X[j]<-weights%*%save[1:N]
}
EVSI.men[i]<-mean(pmax(0,mu.X),na.rm=TRUE)-max(0,mean(mu.X,na.rm=TRUE))
}

####JALAL ET AL####
source('~/predict_ga.R', encoding = 'WINDOWS-1252')
library(R2jags)
library(R2OpenBUGS)
INB<-20000*m$delta.e-m$delta.c
mod<-gam(INB~s(u.l1.noae)+s(u.l1.withdraw.noae))
###Estimating n0
Size.Outer<-1000
Size.Inner<-10000

parameters.to.save <-cbind("u.l1.noae","u.l1.withdraw.noae")
#Utilities
s.u.l1.noae<-beta.par(0.695000000)$alpha
r.u.l1.noae<-beta.par(0.695000000)$beta
#Withdraw due to other reasons
s.u.l1.withdraw.noae<-beta.par(0.405000000)$alpha
r.u.l1.withdraw.noae<-beta.par(0.405000000)$beta
#Model
model<-function(){
  for(i in 1:n.model){
    X2[i]~dbeta(s.X2,r.X2)
    X4[i]~dbeta(s.X4,r.X4)
  }
  
  s.X2<- u.l1.noae*( (u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae)) -1 )
  r.X2<- (1-u.l1.noae)*((u.l1.noae*(1-u.l1.noae)/(sig.X.noae*sig.X.noae)) -1 )
  s.X4<- u.l1.withdraw.noae*( (u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with)) -1 )
  r.X4<- (1-u.l1.withdraw.noae)*((u.l1.withdraw.noae*(1-u.l1.withdraw.noae)/(sig.X.with*sig.X.with)) -1 )
  
  #Utilities
  u.l1.noae~dbeta(s.u.l1.noae,r.u.l1.noae)
  #Withdraw due to other reasons
  u.l1.withdraw.noae~dbeta(s.u.l1.withdraw.noae,r.u.l1.withdraw.noae)
}
filein <- file.path(getwd(),fileext="psitemp.txt")
write.model(model,filein)

Mean.X<-array(NA,dim=c(Size.Outer,2))

start<-Sys.time()
for(j in 1:Size.Outer){
  #CHANGE PROB 0.687
  missingness<-rbinom(40,1,0.687)
  X2<-rbeta(40,betaPar(u.l1.noae[j],sig.X.noae)$a,betaPar(u.l1.noae[j],sig.X.noae)$b)
  X4<-rbeta(40,betaPar(u.l1.withdraw.noae[j],sig.X.with)$a,betaPar(u.l1.withdraw.noae[j],sig.X.with)$b)
  
  X2[which(X2*missingness==0)]<-NA
  X4[which(X4*missingness==0)]<-NA
  
  data<-list(r.u.l1.noae,r.u.l1.withdraw.noae,s.u.l1.noae,s.u.l1.withdraw.noae,
             X2,X4,n.model=40,sig.X.noae,sig.X.with)
  names(data)<-c("r.u.l1.noae","r.u.l1.withdraw.noae","s.u.l1.noae","s.u.l1.withdraw.noae",
                 "X2","X4","n.model","sig.X.noae","sig.X.with")
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  stopping<-1
  while(!(class(
    tryCatch({
      if(stopping>=21){Mean.X[j]<-c(NA,NA)}
      else{
        Model.JAGS<-jags.model(filein, data =  data,
                               inits=list(u.l1.noae=rep(0.8,1),u.l1.withdraw.noae=rep(0.6,1)),
                               n.chains=1,quiet=TRUE)
        update(Model.JAGS,n.iter=200,progress.bar="none")
        coda.save<-coda.samples(Model.JAGS,parameters.to.save,n.iter=Size.Inner,thin=1,progress.bar="none")
      }
    }
    ,error=function(e){
      print(stopping)
    }
    )
  )%in%c("mcmc.list","logical"))){
    X2<-rbeta(40,betaPar(u.l1.noae[j],sig.X.noae)$a,betaPar(u.l1.noae[j],sig.X.noae)$b)
    X2[which(X2*missingness==0)]<-NA
    data<-list(r.u.l1.noae,r.u.l1.withdraw.noae,s.u.l1.noae,s.u.l1.withdraw.noae,
               X2,X4,40,sig.X.noae,sig.X.with)
    names(data)<-c("r.u.l1.noae","r.u.l1.withdraw.noae","s.u.l1.noae","s.u.l1.withdraw.noae",
                   "X2","X4","n.model","sig.X.noae","sig.X.with")
    stopping<-stopping+1
  }
  
  if(stopping<=20){
    coda.full<-coda.save[[1]]
    Mean.X[j,]<-apply(coda.save[[1]],2,mean)
  }
}

end<-Sys.time()

time.Jalal<-end-start

start<-Sys.time()
ae.n0<-40*(var(u.l1.noae)/var(Mean.X[,1],na.rm=TRUE)-1)
noae.n0<-40*(var(u.l1.withdraw.noae)/var(Mean.X[,2],na.rm=TRUE)-1)

evsi.Jal<-array(NA,dim=length(n.))
n0<-c(ae.n0,noae.n0)
for(i in 1:length(n.)){
n<-rep(n.[i],2)
llpred<-predict.ga(mod,n0=n0,n=n)
evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)}
end<-Sys.time()

difftime(end,start,unit="mins")



EVSI.truth<-c(450,643,737,803,835)
EVSI.heath<-c(467,651,750,812,836)
EVSI.jalal<-c(441,621,722,792,819)
EVSI.menzies<-c(467,657,767,816,799)
EVSI.strong<-c(455,650,733,782,802)


which((EVSI.truth-EVSI.heath)/EVSI.truth>=0.05)
which((EVSI.truth-EVSI.jalal)/EVSI.truth>=0.05)
which((EVSI.truth-EVSI.menzies)/EVSI.truth>=0.05)
which((EVSI.truth-EVSI.strong)/EVSI.truth>=0.05)

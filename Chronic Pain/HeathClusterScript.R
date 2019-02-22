#Markov Model for Chronic Pain - Adapted from Sullivan et al.
#Copyright: Anna Heath 2017
#Model has two treatment options for chronic pain and costs and effects are determined using a Markov Model
#Treatment 1: Morphine
#Treatment 2: Novel Treatment
#N: PSA sample size
#NOTE: Throughout the PSA distributions are taken as having the mean given by the parameter estimate and the standard error as 10%
#of the parameter estimate.
#t is treatment and l is the number of treatments that have been attempted

set.seed(1234)
##Packages
library(BCEA)
library(R2OpenBUGS)
library(R2jags)
library(mgcv)
library(psych)
library(foreach)
library(doParallel)

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

##MENZIES N
#N<-5000
##ALL N
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

rm(PSA.Trans.Mat.t1,PSA.Trans.Mat.t2,"AE.l1.t1","AE.l1.t2","AE.l2.t1","AE.l2.t2", "No.AE.l1.t1","No.AE.l1.t2","No.AE.l2.t1","No.AE.l2.t2","p.ae.l1.t1","p.ae.l1.t2","p.ae.l2.t1",
   "p.ae.l2.t2","p.dist.l1","p.dist.l2","p.with.ae.l1.t1","p.with.ae.l1.t2","p.with.ae.l2.t1","p.with.ae.l2.t2",
   "p.with.l1.t1","p.with.l1.t2","p.with.l2.t1","p.with.l2.t2","With.AE.l1.t1",
   "With.AE.l1.t2","With.AE.l2","With.l1.t1","With.l1.t2","With.l2",Dist,Subs.treat,c.ae,c.dist,
   c.med.t1,c.med.t2,c.withdraw,c.withdraw.ae)

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
rm(u.Mat,c.Mat.t1,c.Mat.t2)
####EVPPI####
#This gives the standard cost-effectiveness analysis for the Chronic Pain model.
#It also has the EVPPI calculations to determine where to focus analysis.

discount.15<-sum(1/(1+0.035)^(0:15))
INB<-discount.15*((effects.t1-effects.t2)*20000-(costs.t1-costs.t2))
save<-gam(INB~te(u.l1.noae,u.l1.withdraw.noae))$fitted.values
rm(effects.t1,effects.t2)

betaPar <- function(m,s){
  a <- m*( (m*(1-m)/s^2) -1 )
  b <- (1-m)*( (m*(1-m)/s^2) -1 )
  list(a=a,b=b)
}

###Reduce Size
N<-10000
costs.t2<-costs.t2[1:N]
costs.t1<-costs.t1[1:N]
u.l1.ae<-u.l1.ae[1:N]
u.l1.withdraw.ae<-u.l1.withdraw.ae[1:N]
u.l3<-u.l3[1:N]
u.dist<-u.dist[1:N]
Prob.Array.2<-Prob.Array.2[,,1:N]
Prob.Array.1<-Prob.Array.1[,,1:N]

####HEATH ET AL####
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)
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
rm(u.l1.withdraw.noae,u.l1.noae)

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


no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)
uncert<-200
EVSI.Heath.uncert<-foreach(i=1:(uncert),.combine=rbind,
                             .export=c("phi2","phi4",
                                       "sig.X.noae",
                                       "sig.X.with","betaPar","INB",
                                       "n.","N","save",
                                       "filein","r.u.l1.noae","r.u.l1.withdraw.noae",
                                       "s.u.l1.noae","s.u.l1.withdraw.noae",
                                       "parameters.to.save"),
                            .packages=c("R2jags","R2OpenBUGS")) %dopar% {
                              discount.15<-sum(1/(1+0.035)^(0:15))
                              Var.X.prob<-array(NA,Q)
                              for(j in 1:Q){
                                n.model<-n[j]
                                X2<-rbeta(n.model,betaPar(phi2[j],sig.X.noae)$a,betaPar(phi2[j],sig.X.noae)$b)
                                X4<-rbeta(n.model,betaPar(phi4[j],sig.X.with)$a,betaPar(phi4[j],sig.X.with)$b)
                                
                                #Missingness
                                missingness<-rbinom(n.model,1,0.687)
                                X2[which(X2*missingness==0)]<-NA
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
                                u.l3<-(coda.full[,1]+u.l1.ae[1:N])/2
                                
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
                              rm(effects.t2.X,effects.t1.X)
                              model.ab<-function(){
                                beta~dnorm(Nmax/2,shape.Nmax)%_%T(0,)
                                for(i in 1:N){
                                  y[i]~dnorm(mu[i],tau)
                                  mu[i]<-var.PI*(x[i]/(x[i]+beta))
                                }
                                sigma~dt(sigma.mu,sigma.tau,3)%_%T(0,)#dunif(0.01,50)
                                tau<-1/sigma^2
                              }
                              
                              
                              data.ab<-list(sigma.mu=sd(var(INB)-Var.X.prob)/2,
                                         sigma.tau=1/(sd(var(INB)-Var.X.prob)),
                                         N=length(n),
                                         shape.Nmax=0.0005/max(n),
                                         var.PI=var(save),
                                         Nmax=max(n),
                                         y=as.vector(t(var(INB)-Var.X.prob)),
                                         x=as.vector(rep(n,1)))
                              rm(Var.X.prob)
                              
                              n.chains.ab<-3
                              n.burnin.ab <- 1000  # Number of burn in iterations
                              n.thin.ab<-5
                              n.iter.ab <- ceiling(10000*n.thin/n.chains) + n.burnin # Number of iterations per chain
                              
                              # Choose the parameters in the model to monitor
                              parameters.to.save.ab <- c("beta","sigma","mu")
                              filein.ab <- file.path("~/",fileext="abmodel.txt")
                              write.model(model.ab,filein.ab)
                              
                              # Perform the MCMC simulation with OpenBUGS.
                              # Close OpenBUGS once it has finished (if debug is set to TRUE)
                              bugs.a.b<- jags(
                                data =  data.ab,
                                inits = NULL,
                                parameters.to.save = parameters.to.save.ab,
                                model.file = filein.ab, 
                                n.chains = n.chains.ab, 
                                n.iter = n.iter.ab, 
                                n.thin = n.thin.ab, 
                                n.burnin = n.burnin.ab) 
                              
                              n.<-c(10,25,50,100,150)
                              med<-median(bugs.a.b$BUGSoutput$sims.list$beta)
                              rm(bugs.a.b,data.ab,parameters.to.save.ab,filein.ab,n.chains.ab,
                                 n.iter.ab,n.thin.ab,n.burnin.ab)
                                EVSI.heath<-array(NA,dim=5)
                              var.pred<-var(save)*(n./(n.+med))
                              for(l in 1:5){
                                samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
                                EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
                              }
                          rm(samp.pre)
                              
EVSI.heath
}
stopCluster(cl)
write.table(EVSI.Heath.uncert,"/home/aheath/Chronic_Pain/EVSIHeath.txt")
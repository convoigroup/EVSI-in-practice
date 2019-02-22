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
N<-5000
#N<-100000

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
INB<-discount.15*((effects.t1-effects.t2)*20000-(costs.t1-costs.t2))
save<-gam(INB~te(u.l1.noae,u.l1.withdraw.noae))$fitted.values


rm("AE.l1.t1","AE.l1.t2","AE.l2.t1","AE.l2.t2","beta.par","c.ae","c.dist","c.l2","c.l3","c.Mat.t1",
   "c.Mat.t2","c.med.t1","c.med.t2","c.t1","c.t2","c.withdraw","c.withdraw.ae","cost.novel","costs.t1",
   "costs.t2","Dist","effects.t1","effects.t2","gamma.par","i","Inflation","InitVector","Markov_Prob",
   "No.AE.l1.t1","No.AE.l1.t2","No.AE.l2.t1","No.AE.l2.t2","p.ae.l1.t1","p.ae.l1.t2","p.ae.l2.t1",
   "p.ae.l2.t2","p.dist.l1","p.dist.l2","p.with.ae.l1.t1","p.with.ae.l1.t2","p.with.ae.l2.t1","p.with.ae.l2.t2",
   "p.with.l1.t1","p.with.l1.t2","p.with.l2.t1","p.with.l2.t2","PriceIndex0910","PriceIndex1213","Prob.Array.1",
   "Prob.Array.2","prob.novel","PSA.Trans.Mat.t1","PSA.Trans.Mat.t2","Subs.treat","Time_Horizen","u.dist",
   "u.l1.withdraw.ae","u.l1.ae","u.l2","u.l3","u.Mat","With.AE.l1.t1",
   "With.AE.l1.t2","With.AE.l2","With.l1.t1","With.l1.t2","With.l2","discount.15")

betaPar <- function(m,s){
  a <- m*( (m*(1-m)/s^2) -1 )
  b <- (1-m)*( (m*(1-m)/s^2) -1 )
  list(a=a,b=b)
}

####MENZIES####
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)


sim_data<-function(u.noae,u.withdraw,n,sig.X.noae=0.300,sig.X.with=0.310){
  rand<-rbinom(n,1,0.687)
  samp<-rbeta(sum(rand),max(0,betaPar(u.noae,sig.X.noae)$a),max(0,betaPar(u.noae,sig.X.noae)$b))
  samp.with<-rbeta(sum(rand),max(0,betaPar(u.withdraw,sig.X.with)$a),max(0,betaPar(u.withdraw,sig.X.with)$b))
  r<-cbind(samp,samp.with)
  return(r)
}
func.dbeta<-function(u.l1.noae,u.l1.withdraw.noae,X,sig.X.noae=0.300,sig.X.with=0.310){
  noae.par<-betaPar(u.l1.noae,sig.X.noae)
  withdraw.par<-betaPar(u.l1.withdraw.noae,sig.X.with)
  d<-exp(sum(dbeta(X[,1],noae.par$a,noae.par$b,log=T))+
           sum(dbeta(X[,2],withdraw.par$a,withdraw.par$b,log=T)))
  return(d)
}

no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoParallel(cl)
uncert<-32
EVSI.Menzies.uncert<-foreach(i=1:(uncert),.combine=rbind,
                             .export=c("u.l1.withdraw.noae","u.l1.noae",
                                       "save","sig.X.noae",
                                       "sig.X.with","betaPar",
                                       "sim_data","func.dbeta","n.","N")) %dopar% {
                                         EVSI.men<-array(NA,dim=length(n.))
                                         for(i in 1:length(n.)){
                                           mu.X<-array(NA,dim=N)
                                           for(j in 1:N){
                                             dat.M<-sim_data(u.l1.noae[j],u.l1.withdraw.noae[j],n.[i])
                                             app<-function(u.l1.noae,u.l1.withdraw.noae){
                                               r<-func.dbeta(u.l1.noae,u.l1.withdraw.noae,X=dat.M)
                                               return(r)
                                             }
                                             likeli<-mapply(app,u.l1.noae=u.l1.noae[1:N],u.l1.withdraw.noae=u.l1.withdraw.noae[1:N])
                                             weights<-likeli/sum(likeli,na.rm=TRUE)
                                             weights[is.na(weights)]<-0
                                             mu.X[j]<-weights%*%save[1:N]
                                             rm(weights)
                                           }
                                           EVSI.men[i]<-mean(pmax(0,mu.X),na.rm=TRUE)-max(0,mean(mu.X,na.rm=TRUE))
                                           rm(mu.X)
                                         }
                                         EVSI.men
                                       }
stopCluster(cl)

write.table(EVSI.Menzies.uncert,"/home/aheath/Chronic_Pain/EVSIMenzies_1.txt")
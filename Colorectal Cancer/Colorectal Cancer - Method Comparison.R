#### 08.1 Load packages ####
library(ggplot2)
library(dampack)
library(psych)

###Need Access to Fernando CRC Model GitHub
setwd("~/GitHub/CRCmodR")

#### 08.2 Load inputs ####
source("report/01_model-inputs.R")
#### 08.2.1 Screening parameters ####
v.params.scr <- gen_params_scr()
# list2env(as.list(v.params.scr),globalenv()) # Add parameters to global environment
### Screening frequency
p.Screen <- numeric(length = length(v.ages))
names(p.Screen) <- v.ages
p.Screen[which(v.ages %in% c(50, 60, 70, 80))] <- 1 # Only three
### Surveillance frequency
## Low risk (LR)
p.Surv.LR <- numeric(length = length(v.ages))
names(p.Surv.LR) <- v.ages
p.Surv.LR[which(v.ages %in% seq(50, 85, by = 5))] <- 1
## High risk (HR)
p.Surv.HR <- numeric(length = length(v.ages))
names(p.Surv.HR) <- v.ages
p.Surv.HR[which(v.ages %in% seq(50, 85, by = 3))] <- 1



#### 08.3 Load functions ####
source("R/00_general_functions.R")
source("R/02_crc-nhm_functions.R")
source("R/04_calibrate-det_functions.R")
source("R/08_crc-screening_functions.R")

n.sim <- 5000
### Sample from informed priors
df.prior.informed <- sample.prior.informed(n.sim)
### Load PSA dataset
# df.psa.params <- gen_psa()
# df.psa.params.informed <- df.psa.params[1:n.sim, ]
### Substitute posterior of calibrated parameters with MAP estimates
df.psa.params <- gen_psa()
df.psa.params <- df.psa.params[1:n.sim, ]
df.psa.params.informed <- df.psa.params
df.psa.params.informed[, 1:length(v.params0)] <- df.prior.informed
### Initialize matrices for Costs and QALYs
m.c.informed <- m.e.informed <-matrix(NaN, nrow = n.sim, ncol = 2)
colnames(m.c.informed) <- colnames(m.e.informed) <- v.names.strategies
### Run PSA
for(i in 1:n.sim){ # i <- 1
  out.cea <- cea_crc_screening(v.params.scr = df.psa.params.informed[i, ])
  ### Fill Costs and QALYs matrices
  m.c.informed[i, ] <- out.cea$Costs
  m.e.informed[i, ] <- out.cea$QALYs
  #cat('\r', paste(round(i/n.sim * 100), "% done", sep = " "))
}

library(BCEA)
wtp.vec<-seq(0,120000,by=500)
m<-bcea(m.e.informed,m.c.informed,ref=2,interventions = v.names.strategies,wtp=wtp.vec)
plot(m)

wtp<-38000
which(m$k==wtp)

evppi<-list()
for(i in 1:dim(df.psa.params.informed)[2]){
  evi<-evppi(i,df.psa.params.informed,m)
  evppi[[i]]<-evi$evppi
}

imp<-unlist(lapply(evppi,function(x){x[which(m$k==wtp)]}))

names(imp)<-names(df.psa.params.informed)
imp

###SET UP
library(R2OpenBUGS)
library(R2jags)
library(rjags)
library(mgcv)

#INB for specific willingness-to-pay
INB<-m$ib[which(m$k==wtp),]

#EVPPI
evi.key<-evppi(c("lambda_1","g"),df.psa.params.informed,m)
evi.key$evppi[which(evi.key$k==wtp)]
#Fitted Values for specific WTP
save<-wtp*evi.key$fitted.effects[,1]-evi.key$fitted.costs[,1]

#Sample sizes
n.<-c(5,40,100,200,500,750,1000,1500)

#Bayesian Model for Parameters
model<-function(){
  g~dlnorm(1.039825,31.98693)
  lambda_1~dlnorm(-11.97115,2.90804)
  for(i in 1:N){
    a[i]~dbern(mu[i])
    mu[i]<-1-exp(-lambda_1*g*age[i]^(g-1))
  }
}
filein.HE <- file.path("psitemp.txt")
write.model(model,filein.HE)

##Heath et al.
Q<-50
#Quantiles for Model Parameters
lambda_1.q<-sample(quantile(df.psa.params.informed[,1],probs=(1:Q)/(1+Q)))
g.q<-sample(quantile(df.psa.params.informed[,2],probs=(1:Q)/(1+Q)))

#Sample sizes for Heath et al fitting
N<-round(seq(sqrt(40),sqrt(1500),length.out=Q)^2)
#Ages of patients in the study
ages<-sample(50:90,max(N),replace=TRUE)


###ONLY RUN ONCE TO GET NESTED SIMULATION RESULTS
var.q<-array(NA,dim=c(Q,3))
start<-Sys.time()
n.sim<-5000
for(q in 1:Q){
  #Simulate Data
  lambda_1<-lambda_1.q[q]
  g<-g.q[q]
  mu<-1-exp(-lambda_1*g*ages^(g-1))
  dat<-rbinom(N[q],1,mu)
  
  #Data in Bayesian model
  data.jags<-list(N=N[q],a=dat,age=ages[1:N[q]])
  
  #Bayesian Updating
  Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                       n.chains=1,quiet=TRUE)
  update(Model.JAGS,n.iter=200,progress.bar="none")
  coda.save<-coda.samples(Model.JAGS,c("lambda_1","g"),n.iter=n.sim,thin=1,progress.bar="none")
  
  psa.posterior<-df.psa.params.informed[1:n.sim,]
  psa.posterior[,c(1,2)]<-as.matrix(coda.save)[,c(2,1)]
  
  m.c.post<-m.e.post<-array(NA,dim=c(n.sim))
  for(i in 1:n.sim){
    out.cea <- cea_crc_screening(v.params.scr = psa.posterior[i, ])
    ### Fill Costs and QALYs matrices
    m.c.post[i] <- out.cea$Costs[2]-out.cea$Costs[1]
    m.e.post[i] <- out.cea$QALYs[2]-out.cea$QALYs[1]
    }
  var.q[q,]<-c(var(m.e.post),cov(m.e.post,m.c.post),var(m.c.post))
print(q)
}

var.INB<-wtp^2*var.q[,1]+var.q[,3]-2*wtp*var.q[,2]

#Model fitting to estimate EVSI across sample size
model.ab<-function(){
  beta~dnorm(Nmax/2,shape.Nmax)%_%T(0,)
  for(i in 1:N){
    y[i]~dnorm(mu[i],tau)
    mu[i]<-var.PI*(x[i]/(x[i]+beta))
  }
  sigma~dt(sigma.mu,sigma.tau,3)%_%T(0,)#dunif(0.01,50)
  tau<-1/sigma^2
}

data<-list(sigma.mu=sd(var(m$ib[which(m$k==wtp),])-var.INB)/2,
           sigma.tau=1/(sd(var(m$ib[which(m$k==wtp),])-var.INB)),
           N=length(N),
           shape.Nmax=0.0005/max(N),
           var.PI=var(save),
           Nmax=max(N),
           y=as.vector(t(var(m$ib[which(m$k==wtp),])-var.INB)),
           x=as.vector(N))

#Write model file
filein <- file.path("~/",fileext="abmodel.txt")
write.model(model.ab,filein)

# Perform the MCMC simulation with OpenBUGS.
# Close OpenBUGS once it has finished (if debug is set to TRUE)
bugs.a.b<- jags(
  data =  data,
  inits = NULL,
  parameters.to.save = c("beta","sigma","mu"),
  model.file = filein, 
  n.chains = 3, 
  n.iter = 11000, 
  n.thin = 3, 
  n.burnin = 1000) 

#Estimate the EVSI across sample size with the Heath method.
#Doesn't estimate uncertainty although this is possible - how to compare with other methods?
med<-median(bugs.a.b$BUGSoutput$sims.list$beta)
EVSI.heath<-array(NA,dim=length(n.))
var.pred<-var(save)*(n./(n.+med))
for(l in 1:length(n.)){
  samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
  EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
}

end<-Sys.time()
time.Heath<-end-start

###STRONG ET AL...
EVSI.strong<-array(NA,dim=length(n.))

#Estimate across sample size
start<-Sys.time()
i<-1
for(N in n.){
  sum.stats<-array(NA,dim=c(n.sim,2))
  for(j in 1:n.sim){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Simulate enrollement age for patients - resimulate each time to represent a random study of this design.
    ages<-sample(50:90,N,replace=TRUE)
    mu<-1-exp(-lambda_1*g*ages^(g-1))
    dat<-rbinom(N,1,mu)
    #Sufficient statistic come from likelihood and factorisation theorem.
    sum.stats[j,]<-c(sum(dat),sum(1-dat))
  }
  strong.mod<-gam(INB~te(sum.stats[,1],sum.stats[,2]))
  EVSI.strong[i]<-mean(pmax(fitted(strong.mod),0))-max(0,mean(INB))
  i<-i+1
}
end<-Sys.time()

time.Strong<-end-start

###JALAL ET AL...
start<-Sys.time()
mean.param<-array(NA,dim=c(n.sim,2))
for(q in 1:n.sim){
  lambda_1<-df.psa.params.informed[q,1]
  g<-df.psa.params.informed[q,2]
  #Resimulate ages each time to represent full uncertainty in trial design.
  ages<-sample(50:90,40,replace=TRUE)
  mu<-1-exp(-lambda_1*g*ages^(g-1))
  dat<-rbinom(40,1,mu)
  
  data.jags<-list(N=40,a=dat,age=ages)
  
  Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                         n.chains=1,quiet=TRUE)
  update(Model.JAGS,n.iter=200,progress.bar="none")
  coda.save<-coda.samples(Model.JAGS,c("lambda_1","g"),n.iter=n.sim,thin=1,progress.bar="none")
  
  mean.param[q,]<-apply(as.matrix(coda.save),2,mean)[c(2,1)]
}

start<-Sys.time()
source('predict_ga.R', encoding = 'WINDOWS-1252')
#Might work better with te...can't get the function to work??
mod<-gam(INB~s(df.psa.params.informed[1:n.sim,1])+s(df.psa.params.informed[1:n.sim,2]))

lambda_1.n0<-40*(var(df.psa.params.informed[,1])/var(mean.param[,1],na.rm=TRUE)-1)
g.n0<-40*(var(df.psa.params.informed[,2])/var(mean.param[,2],na.rm=TRUE)-1)

evsi.Jal<-array(NA,dim=length(n.))
n0<-c(lambda_1.n0,g.n0)
for(i in 1:length(n.)){
  n<-rep(n.[i],2)
  llpred<-predict.ga(mod,n0=n0,n=n)
  evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)
}
end<-Sys.time()
time.Jalal<-"31 minutes"
####MENZIES...
likeli<-function(ya,lambda,g){
  y<-ya[,1]
  a<-ya[,2]
  term1<-y*log(1-exp(-lambda*g*a^(g-1)))
  term2<-(1-y)*lambda*g*a^(g-1)
  #Likelihood must be on likelihood scale NOT log scale.
  return(exp(sum(term1-term2)))
}

start<-Sys.time()
mu.X<-array(NA,dim=n.sim)
EVSI.men<-array(NA,dim=length(n.))
n.sim.menz<-2500
for(i in 1:length(n.)){
  for(j in 1:n.sim.menz){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Resimulate ages
    ages<-sample(50:90,n.[i],replace=TRUE)
    mu<-1-exp(-lambda_1*g*ages^(g-1))
    dat<-rbinom(n.[i],1,mu)
    #Set to a multiple function to use mapply instead of nested loop
    app<-function(lambda,g){
      r<-likeli(ya=cbind(dat,ages),lambda,g)
      return(r)
    }
    #The likelihood of data conditional on each parameter set.
    like<-mapply(app,lambda=df.psa.params.informed[1:n.sim.menz,1],g=df.psa.params.informed[1:n.sim.menz,2])
    weights<-like/sum(like)
    mu.X[j]<-weights%*%save[1:n.sim.menz]
  }
  EVSI.men[i]<-mean(pmax(0,mu.X),na.rm=TRUE)-max(0,mean(mu.X,na.rm=TRUE))
  print(i)
}
end<-Sys.time()

time.Menzies<-end-start
EVSI.full<-cbind(evsi.Jal,EVSI.heath,EVSI.strong,EVSI.men)
c(time.Jalal,time.Heath,time.Strong,time.Menzies)
library(colorspace)
cols<-rainbow_hcl(4)

plot(log(n.),evsi.Jal,ylim=c(min(EVSI.full),max(EVSI.full)),pch=19,col=cols[2],xaxt="n",
     xlab="Sample Size (Log-scale)",ylab="EVSI",cex.axis=0.7)
axis(1,at=log(n.),labels=n.,cex.axis=0.7)
points(log(n.),EVSI.heath,pch=2,col=cols[1])
points(log(n.),EVSI.strong,pch=3,col=cols[4])
points(log(n.),EVSI.men,pch=4,col=cols[3])
abline(h=evi.key$evppi[which(evi.key$k==wtp)],lwd=2)
legend("topleft",c("Heath","Jalal","Menzies","Strong","EVPPI"),col=c(cols,"black"),pch=c(2,19,4,3,NA),lwd=c(rep(NA,4),2))

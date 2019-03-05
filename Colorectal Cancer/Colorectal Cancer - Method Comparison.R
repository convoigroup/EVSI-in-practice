#### 08.1 Load packages ####
library(ggplot2)
library(dampack)
library(psych)
library(foreach)

beta_params <- function(mean, sigma){
  alpha <- ((1-mean) / sigma^2 - 1 / mean) * mean^2 
  beta  <- alpha*(1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}


#### 08.2 Load inputs ####
setwd("C:/Users/anna heath/Documents/GitHub/CRCmodR")
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
n.sim<-50000
source("R/00_general_functions.R")
source("R/02_crc-nhm_functions.R")
source("R/04_calibrate-det_functions.R")
source("R/08_crc-screening_functions.R")
source("C:/Users/anna heath/Documents/GitHub/EVSI_Comparison/Colorectal Cancer/09_EVSI_Functions.R")
set.seed(12)
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
out.prior<-foreach(i=1:n.sim,.combine=rbind,.export=ls(globalenv()),.packages=c("msm")) %do% {
  row.pick<-df.psa.params.informed[i, ]
  out.cea <- cea_crc_screening(v.params.scr = row.pick)
  ### Fill Costs and QALYs matrices
  c(out.cea$QALYs,out.cea$Costs)
}
m.c.informed<-out.prior[,3:4]
m.e.informed<-out.prior[,1:2]

n.sim<-5000
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


wtp<-75000

evi<-list()
for(i in 1:dim(df.psa.params.informed)[2]){
  evi[[i]]<-evppi(i,df.psa.params.informed,m)
}

cbind(names(df.psa.params.informed),unlist(lapply(evi,function(x){x$evppi[which(x$k==wtp)]})))

###SET UP
library(R2OpenBUGS)
library(R2jags)
library(rjags)
library(mgcv)
library(dfoptim)
library(demography)

#INB for specific willingness-to-pay
INB<-m$ib[which(m$k==wtp),]

#EVPPI
evi.key<-evppi(c("lambda_1","g","prev.adeno"),df.psa.params.informed[1:5000,],m)
plot(evi.key)
evi.key$evppi[which(evi.key$k==wtp)]

#Fitted Values for specific WTP
save<-wtp*evi.key$fitted.effects[,1]-evi.key$fitted.costs[,1]

write.csv(cbind(INB,df.psa.params.informed[1:5000,c("lambda_1","g","prev.adeno")]),
         "C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/CRCModel_Data.csv")

var(save)/var(INB)
#Sample sizes
n.<-c(5,40,100,200,500,750,1000,1500)

#Age Distribution from Canadian demographic data 2011.
tables<-hmd.mx2("CAN",username="anna.heath@sickkids.ca",password="1546527368")
probs.age<-tables$pop$total[,"2011"]/sum(tables$pop$total[,"2011"])
ages.vec<-0:110
ages<-rmultinom(15000,1,probs.age)
ages.long<-array(NA,dim=15000)
for(i in 1:15000){
ages.long[i]<-ages.vec[which(ages[,i]==1)]}
ages.long<-ages.long[which((ages.long>=25)&(ages.long<=90))]

#Bayesian Model for Parameters
model<-function(){
  #Priors come from initial analysis - normal distribution on a chosen scale
  params~ dmnorm(Mu[1:2],InvCov[1:2,1:2] )
  lambda_1<-params[1]^pow
  g<-params[2]^pow
  for(i in 1:N){
    a[i]~dbern(mu[i])
    mu[i]<-1-exp(-lambda_1*age[i]^g)
  }
}
filein.HE <- file.path("C:/Users/anna heath/Documents/psitemp.txt")
write.model(model,filein.HE)

##Heath et al.
print("Heath")
Q<-50
#Quantiles for Model Parameters
#Data in Bayesian model
ages<-sample(ages.long,1500,replace=TRUE)
data.jags<-list(N=1500,age=ages,Mu=lamb_g_mean,InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)
#Bayesian Updating
Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                       n.chains=1,quiet=TRUE)
update(Model.JAGS,n.iter=200,progress.bar="none")
coda.save<-coda.samples(Model.JAGS,c("a"),n.iter=n.sim*5,thin=5,progress.bar="none")
prior.pred<-as.data.frame(coda.save[[1]])
centres<-kmeans(prior.pred,centers=Q,iter.max=20)$cluster
clusters<-array(NA,dim=c(Q,1500))
for(l in 1:Q){
  clusters[l,]<-apply(prior.pred[which(centres==l),],2,quantile,probs=0.5,type=1)
}

library(foreach)
library(doParallel)
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)

#lambda_1.q<-sample(quantile(df.psa.params.informed[,1],probs=(1:Q)/(1+Q)))
#g.q<-sample(quantile(df.psa.params.informed[,2],probs=(1:Q)/(1+Q)))

#Sample sizes for Heath et al fitting
N<-round((seq(sqrt(5),sqrt(1500),length.out=Q))^2)

###ONLY RUN ONCE TO GET NESTED SIMULATION RESULTS
var.q<-array(NA,dim=c(Q,3))
start<-Sys.time()
for(q in 1:Q){
  #Simulate Data
  #lambda_1<-lambda_1.q[q]
  #g<-g.q[q]
  #Ages of patients in the study
  #ages<-sample(ages.long,N[q],replace=TRUE)
  
  #Probability of adenoma in one year.
  #mu<-1-exp(-lambda_1*ages^g)
  #dat<-rbinom(N[q],1,mu)
  
  #Data in Bayesian model
  samps<-sample(1:1500,N[q])
  data.jags<-list(N=N[q],a=clusters[q,samps],age=ages[samps],Mu=lamb_g_mean,InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)
  
  #Bayesian Updating
  Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                       n.chains=1,quiet=TRUE)
  update(Model.JAGS,n.iter=200,progress.bar="none")
  coda.save<-coda.samples(Model.JAGS,c("lambda_1","g"),n.iter=n.sim*5,thin=5,progress.bar="none")
  
  psa.posterior<-df.psa.params.informed[1:n.sim,]
  psa.posterior[,c(1,2)]<-as.matrix(coda.save)[,c(2,1)]
  psa.posterior[,"prev.adeno"] <- pweibull(50,shape=psa.posterior[,2],
                                          scale=psa.posterior[,1]^(-1/psa.posterior[,2]))
  
  out.post<-foreach(i=1:(n.sim),.combine=rbind,.export=ls(globalenv()),.packages=c("msm")) %dopar% {
    row.pick<-psa.posterior[i,]
    out.cea <- cea_crc_screening(v.params.scr = row.pick)
    ### Fill Costs and QALYs matrices
    c(out.cea$QALYs[2]-out.cea$QALYs[1],out.cea$Costs[2]-out.cea$Costs[1])
  }
  var.q[q,]<-  unique(as.numeric(var(out.post)))#c(var(m.e.post,na.rm=TRUE),cov(m.e.post,m.c.post,use="complete.obs"),var(m.c.post,na.rm=TRUE))
  print(q)
  }
stopCluster()
#write.csv(var.q,"C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/var_heath_fer.csv")

#Read variance from file 
#var.q<-read.csv("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/var_heath_fer.csv")[,2:4]

var.INB<-wtp^2*var.q[,1]+var.q[,3]-2*wtp*var.q[,2]
wtp^2*var(m$delta.e)+var(m$delta.c)-2*wtp*cov(m$delta.c,m$delta.e)

var(save)/var(INB)

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

data<-list(sigma.mu=sd(var(INB)-var.INB)/2,
           sigma.tau=1/(sd(var(INB)-var.INB)),
           N=length(N),
           shape.Nmax=0.0005/max(N),
           var.PI=var(save),
           Nmax=max(N),
           y=as.vector(t(var(INB)-var.INB)),
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
med<-quantile(bugs.a.b$BUGSoutput$sims.list$beta,probs=0.5)
EVSI.heath<-array(NA,dim=length(n.))
var.pred<-var(save)*(n./(n.+med))
for(l in 1:length(n.)){
  samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
  EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
}

end<-Sys.time()
time.Heath<-end-start

###STRONG ET AL...
print("Strong")
EVSI.strong<-array(NA,dim=length(n.))
likeli.s<-function(params,y,a){
  lambda<-params[1]/1e4
  g.<-params[2]
  p<-1-exp(-lambda*a^g.)
  neg.ll<-(prod(p^y*(1-p)^(1-y)))
  return(as.numeric(neg.ll))
}
trans<-as.data.frame(t(df.psa.params.informed[,1:2]))
colnames(trans)<-paste("X",1:n.sim,sep="")

#Estimate across sample size
start<-Sys.time()
i<-1
for(N in n.){
  sum.stats<-array(NA,dim=c(n.sim,3))
  for(j in 1:n.sim){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Simulate enrollement age for patients - resimulate each time to represent a random study of this design.
    ages<-sample(ages.long,N,replace=TRUE)
    mu<-1-exp(-lambda_1*ages^g)
    dat<-rbinom(N,1,mu)
    sum.stats[j,1:2]<- hjkb(c(lambda_1*1e4,g),likeli.s,y=dat,a=ages,
           lower=c(0.00000012,1.22),upper=c(59.06,5.7),control=list(tol=1e-5,maximize=TRUE))$par*c(1e-4,1)
    sum.stats[j,3]<-pweibull(50,shape=sum.stats[j,2],
                             scale=sum.stats[j,1]^(-1/sum.stats[j,2]))
  }
  strong.mod<-gam(INB~te(sum.stats[,1],sum.stats[,2],sum.stats[,3]))
  mean(pmax(fitted(strong.mod),0))-max(0,mean(fitted(strong.mod)))
  EVSI.strong[i]<-mean(pmax(fitted(strong.mod),0))-max(0,mean(fitted(strong.mod)))
  i<-i+1
}
end<-Sys.time()

time.Strong<-end-start
write.csv(sum.stats,"C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/sum_stats_strong_fer.csv")

###JALAL ET AL
print("Jalal")
start<-Sys.time()
mean.param<-array(NA,dim=c(n.sim,3))
for(q in 1:n.sim){
  lambda_1<-df.psa.params.informed[q,1]
  g<-df.psa.params.informed[q,2]
  #Resimulate ages each time to represent full uncertainty in trial design.
  ages<-sample(ages.long,40,replace=TRUE)
  mu<-1-exp(-lambda_1*ages^g)
  dat<-rbinom(40,1,mu)
  
  data.jags<-list(N=40,a=dat,age=ages,Mu=lamb_g_mean,InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)
  
  Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                         n.chains=1,quiet=TRUE)
  update(Model.JAGS,n.iter=200,progress.bar="none")
  coda.save<-coda.samples(Model.JAGS,c("lambda_1","g"),n.iter=n.sim,thin=1,progress.bar="none")
  
  mean.param[q,1:2]<-apply(as.matrix(coda.save),2,mean)[c(2,1)]
  mean.param[q,3]<-mean(pweibull(50,shape=as.matrix(coda.save)[,1],
                           scale=as.matrix(coda.save)[,2]^(-1/as.matrix(coda.save)[,1])))
}


#write.csv(mean.param,"C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/mean_jalal_fer.csv")
#mean.param<-read.csv("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/mean_jalal_fer.csv")[,2:3]
source('C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/predict_ga.R', encoding = 'WINDOWS-1252')
#Might work better with te...can't get the function to work??
mod<-gam(INB[1:n.sim]~s(df.psa.params.informed[1:n.sim,1])+s(df.psa.params.informed[1:n.sim,2])+
           s(df.psa.params.informed[1:n.sim,"prev.adeno"]))

mean.param<-read.csv("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/mean_jalal_fer.csv")[,2:4]
plot(mod$fitted.values,mod$fitted.values-INB[1:n.sim])
lambda_1.n0<-40*(var(df.psa.params.informed[,1])/var(mean.param[,1],na.rm=TRUE)-1)
g.n0<-40*(var(df.psa.params.informed[,2])/var(mean.param[,2],na.rm=TRUE)-1)
prev.n0<-40*(var(df.psa.params.informed[,"prev.adeno"])/var(mean.param[,3],na.rm=TRUE)-1)

evsi.Jal<-array(NA,dim=length(n.))
n0<-c(lambda_1.n0,g.n0,prev.n0)
for(i in 1:length(n.)){
  n<-rep(n.[i],3)
  llpred<-predict.ga(mod,n0=n0,n=n)
  evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)
}
end<-Sys.time()
time.Jalal<-end-start

####MENZIES
print("Menzies")
likeli<-function(ya,lambda,g.){
  y<-ya[,1]
  a<-ya[,2]
  p<-1-exp(-lambda*a^g.)
  neg.ll<-prod(p^y*(1-p)^(1-y))
  return(as.numeric(neg.ll))
}

start<-Sys.time()
mu.X<-array(NA,dim=c(n.sim,length(n.)))
EVSI.men<-array(NA,dim=length(n.))
n.sim.menz<-2500
for(i in 1:length(n.)){
  for(j in 1:n.sim.menz){
    lambda_1<-df.psa.params.informed[j,1]
    g<-df.psa.params.informed[j,2]
    #Resimulate ages
    ages<-sample(ages.long,n.[i],replace=TRUE)
    mu<-1-exp(-lambda_1*ages^g)
    dat<-rbinom(n.[i],1,mu)
    #Set to a multiple function to use mapply instead of nested loop
    app<-function(lambda,g){
      r<-likeli(ya=cbind(dat,ages),lambda,g)
      return(r)
    }
    #The likelihood of data conditional on each parameter set.
    like<-mapply(app,lambda=df.psa.params.informed[1:n.sim.menz,1],g=df.psa.params.informed[1:n.sim.menz,2])
    weights<-like/sum(like)
    mu.X[j,i]<-weights%*%save[1:n.sim.menz]
  }
  EVSI.men[i]<-mean(pmax(0,mu.X[,i]),na.rm=TRUE)-max(0,mean(mu.X[,i],na.rm=TRUE))
  print(i)
}
end<-Sys.time()

time.Menzies<-end-start
write.csv(mu.X,"C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/rescaled_Menzies_fer.csv")

EVSI.full<-cbind(evsi.Jal,EVSI.heath,EVSI.strong,EVSI.men)
time.full<-c(time.Jalal,time.Heath,time.Strong,time.Menzies)

library(colorspace)
cols<-rainbow_hcl(4)
n.<-c(5,40,100,200,500,750,1000,1500)

plot(log(n.),evsi.Jal,ylim=c(min(EVSI.full,na.rm=TRUE),max(EVSI.full,na.rm=TRUE)),pch=19,col=cols[2],xaxt="n",
     xlab="Sample Size (Log-scale)",ylab="EVSI",cex.axis=0.7)
axis(1,at=log(n.),labels=n.,cex.axis=0.7)
points(log(n.),EVSI.heath,pch=2,col=cols[1])
points(log(n.),EVSI.strong,pch=3,col=cols[4])
points(log(n.),EVSI.men,pch=4,col=cols[3])
abline(h=evi.key$evppi[which(evi.key$k==wtp)],lwd=2)
legend("topleft",c("Heath","Jalal","Menzies","Strong","EVPPI"),col=c(cols,"black"),pch=c(2,19,4,3,NA),lwd=c(rep(NA,4),2))

rm(coda.save,data.jags,df.prior.informed,evi,evppi,lambda_asr,llpred,m.c.informed,m.e.informed,mod,Model.JAGS,mu.asr,out.cea,
   psa.posterior,a.,age.init,ages,d.c,d.e,dat,end,filein.HE,g,g.n0,g.q,gap,i,imp,j,lambda_1,lambda_1.n0,lambda_1.q,lambda_7,
   lambda_8,like,m.c.post,m.e.post,mu,mu.X,N,n,n.,n.s,n.s.scr,n.sim.menz,n.t,n0,p.Screen,p.Surv.HR,p.Surv.LR,prev.precl.early,
   prev.precl.late,Q,start,v.ages,v.ages.crc,v.ages.prev,v.n,v.name.states,v.name.states.scr,v.names.params,v.names.strategies,
   v.params.scr,v.params0,wtp.vec,app,cea_crc_screening,crc_nhm,crc_nhm_hist_tp,crc_nhm_micsim,crc_nhm_tp,crc_screening,
   crc_screening_tp,data_summary,gen_bounds,gen_params_init,gen_params_scr,gen_psa,l_likelihood_all,l_likelihood_all_corr,
   l_likelihood_out_all,l_likelihood_out_subset,l_likelihood_subset,l_post,l_prior_informed,l_prior_uniform,likeli,likeli.s,
   likelihood,model,multiplot,plot_nhm_det_out,plot_nhm_det_out_vs_targets,plot_nhm_micsim_out_or_targets,predict.ga,
   Predict.matrix.tensor.smooth.ga,prior,Probs,sample.prior.informed,sample.prior.unif,samplev,Predict.smooth.ga,weights,
   var.INB,q,n.sim,mean.param,var.q,sum.stats,df.psa.params,strong.mod,coefs,cov.rate,data,fit,fit2,lamb_g_cov,
   m.vcov,prior.knowledge,PSA.coefs,PSA.coefs.scale,rates,tables,a,ages.long,ages.vec,b,cor,filein,g.fit,gamm,k,l,
   lamb,lamb_g_cov,lamb_g_mean,lamb_g_rescale,lambda_1.fit,lastFuncGrad,lastFuncParam,mean.rate,means,med,n.sim.co,params,precs,
   probs.age,rate,rate.se,rates.samp,samp.pre,sd.rate,stuff,var.pred,grad,hmd.mx2,model.ab,bugs.a.b,save,INB,age.fit,cl,clusters,
   out.post,trans,prior.pred,centres,no_cores,prev.n0,samps)

#save.image("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/Colorectal Cancer EVSI Results.RData")
#load("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/Colorectal Cancer EVSI Results.RData")

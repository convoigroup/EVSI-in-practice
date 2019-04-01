## Run baseline model
source("Colorectal Cancer Model.R")
source("predict_ga.R", encoding = 'WINDOWS-1252')

## EVSI estimation - Jalal et al.

# Estimate N_0
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

# Calculate linear meta-model
mod<-gam(INB[1:n.sim]~s(df.psa.params.informed[1:n.sim,1])+s(df.psa.params.informed[1:n.sim,2])+
           s(df.psa.params.informed[1:n.sim,"prev.adeno"]))

# Calculate N_0
lambda_1.n0<-40*(var(df.psa.params.informed[,1])/var(mean.param[,1],na.rm=TRUE)-1)
g.n0<-40*(var(df.psa.params.informed[,2])/var(mean.param[,2],na.rm=TRUE)-1)
prev.n0<-40*(var(df.psa.params.informed[,"prev.adeno"])/var(mean.param[,3],na.rm=TRUE)-1)

# Calculate EVSI
evsi.Jal<-array(NA,dim=length(n.))
n0<-c(lambda_1.n0,g.n0,prev.n0)
for(i in 1:length(n.)){
  n<-rep(n.[i],3)
  llpred<-predict.ga(mod,n0=n0,n=n)
  evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)
}

### EVSI
evsi.Jal

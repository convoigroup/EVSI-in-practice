## Run Health Economic Model
source("Chronic Pain Model.R")

# Data set up
sig.X.noae<-0.300
sig.X.with<-0.310

# Sample size of future data
n.<-c(10,25,50,100,150)

# Bayesian Model
n.chains<-2
n.burnin <- 1000  # Number of burn in iterations
n.thin<-1
n.iter <- ceiling(N*n.thin/n.chains) # Number of iterations per chain

# Choose the parameters in the model to monitor
parameters.to.save <-cbind("u.l1.noae","u.l1.withdraw.noae")  
datanames<-c("r.u.l1.noae","r.u.l1.withdraw.noae","s.u.l1.noae","s.u.l1.withdraw.noae",
             "X2","X4","n.model","sig.X.noae","sig.X.with")

# Quantile Specification
Q<-50
# Sample size for Heath et al. method
n<-trunc(seq(10,150,length.out=Q))
phi2<-quantile(u.l1.noae,prob=(1:Q)/(Q+1))
phi4<-quantile(u.l1.withdraw.noae,prob=(1:Q)/(Q+1))
while(abs(cor(n,phi2))>0.001){phi2<-sample(quantile(u.l1.noae,prob=(1:Q)/(Q+1)),replace=F)}
while(abs(cor(n,phi4))>0.001){phi4<-sample(quantile(u.l1.withdraw.noae,prob=(1:Q)/(Q+1)),replace=F)}

# Reduced Bayesian Model for data update
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

## EVSI Estimation Heath et al.
Var.X.prob<-array(NA,Q)
for(j in 1:Q){
  # Data Generation
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
  
  # Perform the MCMC simulation with JAGS
  # tryCatch needed to skip data that causes issues with JAGS.
  while(class(tryCatch({
    # Try Gibbs sampler
    Model.JAGS<-jags.model(filein, data =  data,
                           inits=list(u.l1.noae=rep(0.8,1),u.l1.withdraw.noae=rep(0.6,1)),
                           n.chains=n.chains,quiet=TRUE)
    update(Model.JAGS,n.iter=n.burnin,progress.bar="none")
    coda.save<-coda.samples(Model.JAGS,parameters.to.save,n.iter=n.iter,thin=n.thin,progress.bar="none")
  },error=function(e) 1))!="mcmc.list"){
    # If error - generate different data
    X2<-rbeta(n.model,betaPar(phi2[j],sig.X.noae)$a,betaPar(phi2[j],sig.X.noae)$b)
    X2[which(X2*missingness==0)]<-NA
    data<-list(r.u.l1.noae,r.u.l1.withdraw.noae,s.u.l1.noae,s.u.l1.withdraw.noae,
               X2,X4,n.model,sig.X.noae,sig.X.with)
    names(data)<-datanames
  }
  
  # Rerun health economic model
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

# Fit Bayesian model to estimate EVSI across sample size
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

n.chains.ab<-3
n.burnin.ab <- 1000  # Number of burn in iterations
n.thin.ab<-5
n.iter.ab <- ceiling(10000*n.thin/n.chains) + n.burnin # Number of iterations per chain

# Choose the parameters in the model to monitor
parameters.to.save.ab <- c("beta","sigma","mu")
filein.ab <- file.path("~/",fileext="abmodel.txt")
write.model(model.ab,filein.ab)

# Perform the MCMC simulation with JAGS
bugs.a.b<- jags(
  data =  data.ab,
  inits = NULL,
  parameters.to.save = parameters.to.save.ab,
  model.file = filein.ab, 
  n.chains = n.chains.ab, 
  n.iter = n.iter.ab, 
  n.thin = n.thin.ab, 
  n.burnin = n.burnin.ab) 


# Summarize results
med<-median(bugs.a.b$BUGSoutput$sims.list$beta)

# Calculate EVSI across sample size
EVSI.heath<-array(NA,dim=5)
var.pred<-var(save)*(n./(n.+med))
for(l in 1:5){
  samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
  EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
}

### EVSI results
EVSI.heath

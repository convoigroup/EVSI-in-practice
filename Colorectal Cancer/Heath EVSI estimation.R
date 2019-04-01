## Run baseline model
source("Colorectal Cancer Model.R")

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

## EVSI calculation - Heath et al
Q<-50

# Generate data for the Bayesian model using clustering

# Sample age of future data points
ages<-sample(ages.long,1500,replace=TRUE)
# Forward sample from data generating process
data.jags<-list(N=1500,age=ages,Mu=lamb_g_mean,InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)

Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                       n.chains=1,quiet=TRUE)
update(Model.JAGS,n.iter=200,progress.bar="none")
coda.save<-coda.samples(Model.JAGS,c("a"),n.iter=n.sim*5,thin=5,progress.bar="none")

# Predictive distribution of data
prior.pred<-as.data.frame(coda.save[[1]])
centres<-kmeans(prior.pred,centers=Q,iter.max=20)$cluster
clusters<-array(NA,dim=c(Q,1500))
for(l in 1:Q){
  # Find the centre of the 50 data clusters.
  clusters[l,]<-apply(prior.pred[which(centres==l),],2,quantile,probs=0.5,type=1)
}

#Sample sizes for Heath et al fitting
N<-round((seq(sqrt(5),sqrt(1500),length.out=Q))^2)

### Nested simulation results
var.q<-array(NA,dim=c(Q,3))
for(q in 1:Q){
  #Data in Bayesian model
  # Randomly simulate which patients to include in data
  samps<-sample(1:1500,N[q])
  data.jags<-list(N=N[q],a=clusters[q,samps],
                  age=ages[samps],Mu=lamb_g_mean,
                  InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)
  
  #Bayesian Updating
  Model.JAGS<-jags.model(filein.HE, data =  data.jags,
                         n.chains=1,quiet=TRUE)
  update(Model.JAGS,n.iter=200,progress.bar="none")
  coda.save<-coda.samples(Model.JAGS,c("lambda_1","g"),n.iter=n.sim*5,thin=5,progress.bar="none")
  
  psa.posterior<-df.psa.params.informed[1:n.sim,]
  psa.posterior[,c(1,2)]<-as.matrix(coda.save)[,c(2,1)]
  psa.posterior[,"prev.adeno"] <- pweibull(50,shape=psa.posterior[,2],
                                           scale=psa.posterior[,1]^(-1/psa.posterior[,2]))
  
  # Rerun health economic model with updated data
  out.post <-matrix(NaN, nrow = n.sim, ncol = 2)
  for(i in 1:n.sim){
    row.pick<-psa.posterior[i,]
    out.cea <- cea_crc_screening(v.params.scr = row.pick)
    out.post[i,] <- c(out.cea$QALYs[2]-out.cea$QALYs[1],out.cea$Costs[2]-out.cea$Costs[1])
    ### Fill Costs and QALYs matrices
    m.c.informed[i, ] <- out.cea$Costs
    m.e.informed[i, ] <- out.cea$QALYs
  }
  var.q[q,]<-  unique(as.numeric(var(out.post)))
}

var.INB<-wtp^2*var.q[,1]+var.q[,3]-2*wtp*var.q[,2]

# Model fitting to estimate EVSI across sample size
model.ab<-function(){
  beta~dnorm(Nmax/2,shape.Nmax)%_%T(0,)
  for(i in 1:N){
    y[i]~dnorm(mu[i],tau)
    mu[i]<-var.PI*(x[i]/(x[i]+beta))
  }
  sigma~dt(sigma.mu,sigma.tau,3)%_%T(0,)
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
med<-quantile(bugs.a.b$BUGSoutput$sims.list$beta,probs=0.5)
EVSI.heath<-array(NA,dim=length(n.))
var.pred<-var(save)*(n./(n.+med))
for(l in 1:length(n.)){
  samp.pre<-(save-mean(save))/sd(save)*sqrt(var.pred[l])+mean(save)
  EVSI.heath[l]<-mean(pmax(samp.pre,0))-max(mean(samp.pre),0)
}

### EVSI
EVSI.heath
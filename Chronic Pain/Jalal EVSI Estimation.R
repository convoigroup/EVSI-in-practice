## Run Health Economic Model
source("Chronic Pain Model.R")

# Load function to use Jalal et al. method
source("predict_ga.R", encoding = 'WINDOWS-1252')

# Linear Meta-Model
mod<-gam(INB~s(u.l1.noae)+s(u.l1.withdraw.noae))

# Data generation
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)

### EVSI estimation Jalal et al.
Size.Outer<-1000
Size.Inner<-10000

# Reduced Model for key parameters
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

# Estimate N_0
Mean.X<-array(NA,dim=c(Size.Outer,2))
for(j in 1:Size.Outer){
  # Generating missingness
  missingness<-rbinom(40,1,0.687)
  # Generate data
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
        # Try Gibbs sampler
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
    # If error - generate different data
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

# Calculate N_0 from nested MC simulation
ae.n0<-40*(var(u.l1.noae)/var(Mean.X[,1],na.rm=TRUE)-1)
noae.n0<-40*(var(u.l1.withdraw.noae)/var(Mean.X[,2],na.rm=TRUE)-1)

# Calculate EVSI across sample size
evsi.Jal<-array(NA,dim=length(n.))
n0<-c(ae.n0,noae.n0)
for(i in 1:length(n.)){
  n<-rep(n.[i],2)
  llpred<-predict.ga(mod,n0=n0,n=n)
  evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)
  }

### EVSI
evsi.Jal

## Run Health Economic Model
source("Chronic Pain Model.R")

## Reduce Menzies sample size
N <- 5000

# Data set up
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)

# Function to simulate data
sim_data<-function(u.noae,u.withdraw,n,sig.X.noae=0.300,sig.X.with=0.310){
  rand<-rbinom(n,1,0.687)
  samp<-rbeta(sum(rand),betaPar(u.noae,sig.X.noae)$a,betaPar(u.noae,sig.X.noae)$b)
  samp.with<-rbeta(sum(rand),betaPar(u.withdraw,sig.X.with)$a,betaPar(u.withdraw,sig.X.with)$b)
  r<-cbind(samp,samp.with)
  return(r)
}

# Likelihood function
func.dbeta<-function(u.l1.noae,u.l1.withdraw.noae,X,sig.X.noae=0.300,sig.X.with=0.310){
  noae.par<-betaPar(u.l1.noae,sig.X.noae)
  withdraw.par<-betaPar(u.l1.withdraw.noae,sig.X.with)
  d<-exp(sum(dbeta(X[,1],noae.par$a,noae.par$b,log=T))+
           sum(dbeta(X[,2],withdraw.par$a,withdraw.par$b,log=T)))
  return(d)
}

# EVSI calculation
EVSI.men<-array(NA,dim=length(n.))
# Loop over sample sizes
for(i in 1:length(n.)){
  mu.X<-array(NA,dim=N)
  for(j in 1:N){
    # Simulate data
    dat.M<-sim_data(u.l1.noae[j],u.l1.withdraw.noae[j],n.[i])
    # Create a single input function
    app<-function(u.l1.noae,u.l1.withdraw.noae){
      r<-func.dbeta(u.l1.noae,u.l1.withdraw.noae,X=dat.M)
      return(r)
    }
    # Calculate likelihood
    likeli<-mapply(app,u.l1.noae=u.l1.noae[1:N],u.l1.withdraw.noae=u.l1.withdraw.noae[1:N])
    weights<-likeli/sum(likeli)
    mu.X[j]<-weights%*%save[1:N]
  }
  EVSI.men[i]<-mean(pmax(0,mu.X),na.rm=TRUE)-max(0,mean(mu.X,na.rm=TRUE))
}

### EVSI
EVSI.men

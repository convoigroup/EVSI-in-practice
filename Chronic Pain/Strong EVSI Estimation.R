## Run Health Economic Model
source("Chronic Pain Model.R")

## Data set up.
sig.X.noae<-0.300
sig.X.with<-0.310
n.<-c(10,25,50,100,150)

### EVSI Estimation Strong et al.
EVSI.gam<-array(NA,dim=5)
# Loop over sample size
for(i in 1:5){
  n.loop<-n.[i]
  X.gam<-array(NA,dim=c(N,2))
  X.gam.with<-array(NA,dim=c(N,2))
  for(j in 1:N){
    # Generate data
    rand<-rbinom(n.loop,1,0.687)
    samp<-rbeta(sum(rand),betaPar(u.l1.noae[j],sig.X.noae)$a,betaPar(u.l1.noae[j],sig.X.noae)$b)
    samp.with<-rbeta(sum(rand),betaPar(u.l1.withdraw.noae[j],sig.X.with)$a,betaPar(u.l1.withdraw.noae[j],sig.X.with)$b)
    # Calculate summary statistics - geometric mean
    X.gam[j,]<-c(psych::geometric.mean(samp),psych::geometric.mean(1-samp))
    X.gam.with[j,]<-c(psych::geometric.mean(samp.with),psych::geometric.mean(1-samp.with))
  }
  gam.fit<-gam(INB[1:N]~te(X.gam[,1], X.gam[,2], X.gam.with[,1], X.gam.with[,2])) 
  EVSI.gam[i]<-mean(pmax(0,gam.fit$fitted))-max(0,mean(gam.fit$fitted))
  }

### EVSI
EVSI.gam

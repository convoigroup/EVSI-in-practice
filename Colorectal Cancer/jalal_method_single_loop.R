library(R2jags)
library(mgcv)
setwd("~/Git/EVSI_Comparison/Colorectal Cancer/")
df.psa.params.informed <- read.csv(file = "CRCModel_Data.csv")
#df.psa.params.informed <- read.csv("C:/Users/anna heath/OneDrive/OneDrive - SickKids/EVSI Methods/CRCModel_Data.csv")
INB <- df.psa.params.informed$INB
df.psa.params.informed <- df.psa.params.informed[,-c(1,2)]
prior.lambda_1 <- df.psa.params.informed$lambda_1
prior.g <- df.psa.params.informed$g
prior.prev<-df.psa.params.informed$prev.adeno

n.sim<-length(INB)
source("09_EVSI_functions.R")

#Bayesian Model for Parameters
model<-function(){
 for(q in 1:n.sim){
  params[q,1:2]~ dmnorm(Mu[1:2],InvCov[1:2,1:2] )
  lambda_1[q]<-params[q,1]^pow
  g[q]<-params[q,2]^pow
  for(i in 1:N){
    a[q,i]~dbern(mu[q,i])
    mu[q,i]<-1-exp(-lambda_1[q]*age[q,i]^g[q])
  }
 }

}

n.sim <- length(INB)

###JALAL ET AL...
start<-Sys.time()
mean.param<-array(NA,dim=c(n.sim,3))
nn <- 40
lambda_1<-prior.lambda_1 %*% matrix(1,nrow = 1, ncol = nn)
g<-prior.g %*% matrix(1,nrow = 1, ncol = nn)

#Resimulate ages each time to represent full uncertainty in trial design.
# ? I think the code only uses the first age in each simulation to compute mu
ages<-sample(50:90,n.sim,replace=TRUE) %*% matrix(1,nrow = 1, ncol = nn)
# to allow for mu to vary by age use these two lines
#ages<-sample(50:90,nn*n.sim,replace=TRUE)
#dim(ages) <- c(n.sim, nn)

mu<-1-exp(-lambda_1*ages^g)
dat <- matrix(0, nrow = n.sim, ncol = nn)
for(i in 1:nn){
for(q in 1:n.sim){
  dat[q,i]<-rbinom(1,1,mu[q,i])
}
}
data.jags<-list(N=nn,a=dat,age=ages, n.sim = n.sim,Mu=lamb_g_mean,InvCov=solve(lamb_g_cov),pow=lamb_g_rescale)
model.sim<-jags(model.file = model, data =  data.jags,
                   n.chains=1 ,n.iter=5200, n.burnin=200,
                progress.bar="text", parameters.to.save = c("lambda_1", "g"))
end<-Sys.time()
time.n0<-end-start

post.lambda_1 <- model.sim$BUGSoutput$sims.list$lambda_1
post.g <- model.sim$BUGSoutput$sims.list$g
mean.post.lambda_1 <- colMeans(post.lambda_1)
mean.post.g <- colMeans(post.g)
post.prev<-pweibull(50,shape=post.g,scale=post.lambda_1^(-1/post.g))
mean.post.prev<-colMeans(post.prev)


start<-Sys.time()
source('~/Git/EVSI_Comparison/Chemotherapy/predict_ga.R', encoding = 'WINDOWS-1252')
#Might work better with te...can't get the function to work?? 
mod<-gam(INB~s(prior.lambda_1)+s(prior.g)+s(prior.prev))#+ti(prior.lambda_1,prior.g,prior.prev))

lambda_1.n0<-nn*(var(prior.lambda_1)/var(mean.post.lambda_1)-1)
g.n0<-nn*(var(prior.g)/var(mean.post.g)-1)
prev.n0<-nn*(var(prior.prev)/var(mean.post.prev)-1)
n.<-c(5,40,100,200,500,750,1000,1500)

evsi.Jal<-array(NA,dim=length(n.))
n0<-c(lambda_1.n0,g.n0,prev.n0)
n0
for(i in 1:length(n.)){
  n<-rep(n.[i],3)
  llpred<-predict.ga(mod,n0=n0,n=n)
  evsi.Jal[i] <- mean(pmax(0,llpred))-max(mean(llpred),0)
}
end<-Sys.time()
time.fit<-end-start
evsi.Jal

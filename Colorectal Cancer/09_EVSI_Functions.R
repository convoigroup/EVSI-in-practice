###Calibrating PSA Distribution from Published Data - Plus additional information
#Mean rate and confidence intervals come from published data on the log scale
mean.rate<- log(c(0.00836, 0.00990, 0.01156, 0.01333, 0.01521))
#Calculate the standard deviation from the bottom end of the confidence interval and assume a normal approximation
sd.rate<-(-log(c((0.00418), (0.00495), (0.00578), (0.00667), (0.00761)))+mean.rate)/1.96
age.fit<-c(50,55,60,65,70)

#Create a PSA distribution for the rate
#We assume that the rates for different ages are likely to be increasing and induce a positive correlation 
#between the different rates.
cor<-0.85
cov.rate<-matrix(c(sd.rate[1]*sd.rate[1],cor*sd.rate[1]*sd.rate[2],cor*sd.rate[1]*sd.rate[3],cor*sd.rate[1]*sd.rate[4],cor*sd.rate[1]*sd.rate[5],
                   cor*sd.rate[2]*sd.rate[1],sd.rate[2]*sd.rate[2],cor*sd.rate[2]*sd.rate[3],cor*sd.rate[2]*sd.rate[4],cor*sd.rate[2]*sd.rate[5],
                   cor*sd.rate[3]*sd.rate[1],cor*sd.rate[3]*sd.rate[2],sd.rate[3]*sd.rate[3],cor*sd.rate[3]*sd.rate[4],cor*sd.rate[3]*sd.rate[5],
                   cor*sd.rate[4]*sd.rate[1],cor*sd.rate[4]*sd.rate[2],cor*sd.rate[4]*sd.rate[3],sd.rate[4]*sd.rate[4],cor*sd.rate[4]*sd.rate[5],
                   cor*sd.rate[5]*sd.rate[1],cor*sd.rate[5]*sd.rate[2],cor*sd.rate[5]*sd.rate[3],cor*sd.rate[5]*sd.rate[4],sd.rate[5]*sd.rate[5]),
                 nrow=5,ncol=5)

library(mvtnorm) # library to generate from MVNormal
set.seed(1234) # Set seed to replicate results
n.sim.co<-5000 # How many PSA simulations
coefs<-array(NA, dim=c(n.sim.co,2))
for(j in 1:n.sim.co){
  #Sample from MVNormal to give a possible set of rates
  rates.samp<-as.numeric(exp((rmvnorm(1,mean.rate,cov.rate))))
  #Fit a Weibull hazard to the estimated rates
  fit <- nls(rates.samp ~ lambda_1*g*age.fit^(g-1), 
             start = c(lambda_1 = 0.001, g = 1),control=list(minFactor=2E-15,maxiter=7500))
  #Record the parameter estimates
  coefs[j,]<-coef(fit)}

#Mean Weibull survival curve for sanity check
curve(1-pweibull(x,shape=2.779083721693,scale=0.000002854961^(-1/2.779083721693)),xlim=c(0,100),ylim=c(0,1),col="red",lwd=2)
#Plot small number of the PSA survival curves to check that they are sensible
for(k in 1:150){
  curve(1-pweibull(x,shape=(coefs)[k,2],scale=(coefs)[k,1]^(-1/(coefs)[k,2])),xlim=c(0,100),add=TRUE)
}

#Find transformation that gives linear, approximately normal distribution between the estimated parameters.
#Log-transformation is not successful so not a log-normal distribution
plot((coefs^(1/15)))
hist((coefs[,1]^(1/15)))
hist((coefs[,2]^(1/15)))

#Simulate PSA distribution from a multivariate normal distribution 
#We use this rather than the parameters from the model fit so the Bayesian model can match the 
#prior exactly.
PSA.coefs.scale<-rmvnorm(n.sim,apply(coefs^(1/15),2,mean),cov(coefs^(1/15),use="complete.obs"))
#Rescale to natural scale
PSA.coefs<-PSA.coefs.scale^15
#Check the simulated values match our calibrated values.
plot(PSA.coefs)
points(coefs,col="red")
hist(PSA.coefs[,2],nclass=50,freq=FALSE)
hist(coefs[,2],nclass=50,col="red",add=TRUE,freq=FALSE)

hist(PSA.coefs[,1],nclass=100,freq=FALSE)
hist(coefs[,1],nclass=300,col="red",add=TRUE,freq=FALSE)


lamb_g_mean<-apply(coefs^(1/15),2,mean)
lamb_g_cov<-cov(coefs^(1/15),use="complete.obs")
lamb_g_rescale<-15

###Function to extract mortality tables
hmd.mx2 =  function (country, username, password, label = country) 
{
  path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
                "Mx_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  mx <- try(utils::read.table(con, skip = 2, header = TRUE, 
                              na.strings = "."), TRUE)
  close(con)
  if (class(mx) == "try-error") 
    stop("Connection error at www.mortality.org. Please check username, password and country label.")
  path <- paste("https://www.mortality.org/hmd/", country, "/STATS/", 
                "Exposures_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  pop <- try(utils::read.table(con, skip = 2, header = TRUE, 
                               na.strings = "."), TRUE)
  close(con)
  if (class(pop) == "try-error") 
    stop("Exposures file not found at www.mortality.org")
  obj <- list(type = "mortality", label = label, lambda = 0)
  obj$year <- sort(unique(mx[, 1]))
  n <- length(obj$year)
  m <- length(unique(mx[, 2]))
  obj$age <- mx[1:m, 2]
  mnames <- names(mx)[-c(1, 2)]
  n.mort <- length(mnames)
  obj$rate <- obj$pop <- list()
  for (i in 1:n.mort) {
    obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
    obj$rate[[i]][obj$rate[[i]] < 0] <- NA
    obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, 
                                                              obj$year)
  }
  names(obj$pop) = names(obj$rate) <- tolower(mnames)
  obj$age <- as.numeric(as.character(obj$age))
  if (is.na(obj$age[m])) 
    obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
  return(structure(obj, class = "demogdata"))
}







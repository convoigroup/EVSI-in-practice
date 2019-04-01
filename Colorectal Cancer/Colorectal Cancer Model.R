#### 08.1 Load packages ####
library(dampack)
library(psych)
library(readxl)
library(dplyr)
library(msm)
library(lhs)
library(BCEA)
library(R2OpenBUGS)
library(R2jags)
library(rjags)
library(mgcv)
library(dfoptim)
library(demography)

beta_params <- function(mean, sigma){
  alpha <- ((1-mean) / sigma^2 - 1 / mean) * mean^2 
  beta  <- alpha*(1 / mean - 1)
  params <- list(alpha = alpha, beta = beta)
  return(params)
}

## Generate Model data
#-------------------------------------------------------#
#### GENERATE Initial Set of Parameters from Wu 2006 ####
#-------------------------------------------------------#
gen_params_init <- function(){
  #### Estimate Weibull parameters of transition rate to small adenoma
  ## Rates from paper
  rates <- data.frame(age = c(50, 55, 60, 65, 70),
                      rate = c(0.00836, 0.00990, 0.01156, 0.01333, 0.01521),
                      se = c((0.01672-0.00418), (0.01980-0.00495), (0.02312-0.00578), (0.02666-0.00667), (0.03042-0.00761))/(2*1.96)
  )
  
  fit2<-lm(rate~age,rates,weights= 1/se^2)
  summary(fit2)
  
  ## Fit nonlinear model
  fit <- nls(rate ~ lambda_1*g*age^(g-1), rates, 
             start = c(lambda_1 = 0.001, g = 1), 
             weights = 1/se^2)

  ## Obtain parameters
  params <- coef(fit)
  lambda_1.fit <- params[1]
  g.fit        <- params[2]
  ## Obtain Variance-Covariance matrix
  m.vcov <- vcov(fit)

  #### Transition rates ####
  ## Parameter names
  v.names.params <- c(
    "lambda_1",
    "g",
    "lambda_2",
    "lambda_3",
    "lambda_4",
    "lambda_5",
    "lambda_6",
    # "lambda_7",
    # "lambda_8",
    "prev.adeno",
    "prop.adeno.sm"
  )
  n.params <- length(v.names.params)
  
  ## Deafult values from Wu 2006
  v.params0 <- c(lambda_1 = as.numeric(params[1]),
                 g        = as.numeric(params[2]),
                 lambda_2 = 0.0346,
                 lambda_3 = 0.0215,
                 lambda_4 = 0.3697,
                 lambda_5 = 0.2382,
                 lambda_6 = 0.4852,
                 # lambda_7 = 0.0302, # Non-calibrated parameter
                 # lambda_8 = 0.2099, # Non-calibrated parameter
                 prev.adeno    = 0.27, # Prevalence of adenoma at age 50 from Rutter et al. 2007
                 prop.adeno.sm = 0.71
  )
  return(v.params0)
}

#------------------------------------------------------#
#### GENERATE Base-case Set of Screening Parameters ####
#------------------------------------------------------#
gen_params_scr <- function(){
  v.params.nhm <- gen_params_init()
  v.params.scr <- c(
    ### NHM parameters
    v.params.nhm,
    ### Colonoscopy test characteristics
    ## Small Adenoma
    sens.COL = 0.815,
    spec.COL = 0.868,
    ## Large Adenoma and CRC
    sens.COL.CRC = 0.950,
    spec.COL.CRC = 0.868,
    ### Increase rate from Normal to Small Adenoma with history of lesion
    ## Hazard ratio (HR) for Low risk (LR)
    hr.Hist.LR = 2,
    ## Hazard ratio (HR) for High risk (HR)
    hr.Hist.HR = 3,
    ### Utilities
    ## Normal
    u.Normal = 1,
    ## Clinical CRC Early
    u.CRC_Early = 0.855, # From Reid Ness 1999 AJG (average of Local, 0.95, and regional, 0.76)
    ## Clinical CRC Late
    u.CRC_Late = 0.3,    # From Reid Ness 1999 AJG (Mets of 0.30)
    ### Cost parameters
    ## Colonoscopy
    c.COL = 10000,
    ## Cancer by stages
    c.crc.early = 21524,
    c.crc.late = 37000
  )
  return(v.params.scr)
}


#### 01.3 ASR mortality ####
mu.asr <- read_xlsx("ASR_Females_US_2008.xlsx")
lambda_asr <- mu.asr %>%
  dplyr::filter(Age >= 50 & Age < 100) %>%
  dplyr::select(hazard) %>%
  as.matrix()

#### 01.4 Generate initial Set of Unknown Parameters ####
v.params0 <- gen_params_init()
### Create name of parameters
v.names.params <- names(v.params0)

#### 01.5 External parameters ####
### Number of cycles
n.t   <- 50         # time horizon, 50 cycles
### Initial Age
age.init <- 50
## Vector of Ages
v.ages <- age.init:(age.init + n.t-1)

### Name of health states in microsimulation model
## Full names
v.name.states <- c("Normal",
                   "SmallAdeno",
                   "LargeAdeno",
                   "PreClinCRC_Early",
                   "PreClinCRC_Late",
                   "ClinCRC_Early",
                   "ClinCRC_Late",
                   "CRC_Death",
                   "OC_Death") 
## Short names
v.n <- c("N",
         "SA",
         "LA",
         "PrClE",
         "PrClL",
         "ClE",
         "ClL",
         "CRC_D",
         "OC_D")
n.s <- length(v.n) # Number of states

### Labels Target Ages
v.ages.prev <- seq(52, 97, by = 5)
v.ages.crc  <- c(seq(52, 82, by = 5), 93)

### External transition rates
lambda_7 <- 0.0302 # CRC mortality in early stage
lambda_8 <- 0.2099 # CRC mortality in late stage

### Prevalence of preclinical cancer by stage
prev.precl.early <- 0.0012
prev.precl.late  <- 0.0008

#### 01.6 Screening parameters ####
v.names.strategies <- c("No Screening", "Screening")
### State names and cycles
v.name.states.scr <- c("Normal",
                       "Normal_Hist_LR",
                       "Normal_Hist_HR",
                       "SmallAdeno",
                       "SmallAdeno_LR",
                       "SmallAdeno_HR",
                       "LargeAdeno",
                       "LargeAdeno_LR",
                       "LargeAdeno_HR",
                       "PreClinCRC_Early",
                       "PreClinCRC_Early_LR",
                       "PreClinCRC_Early_HR",
                       "PreClinCRC_Late",
                       "PreClinCRC_Late_LR",
                       "PreClinCRC_Late_HR",
                       "ClinCRC_Early",
                       "ClinCRC_Late",
                       "CRC_Death",
                       "OC_Death") 
n.s.scr <- length(v.name.states.scr) # Number of states
### Discounting factors
d.e <- 0.03
d.c <- 0.03


#### 08.2.1 Screening parameters ####
v.params.scr <- gen_params_scr()

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

#### Simulate PSA distributions ####
n.sim<-5000

gen_bounds <- function(){
  #--- Lower bounds ---#
  v.lb <- c(lambda_1 = 2e-06,
            g        = 2,
            lambda_2 = 0.01,
            lambda_3 = 0.01,
            lambda_4 = 0.2,
            lambda_5 = 0.2,
            lambda_6 = 0.3,
            # lambda_7 = 0.01, # Non-calibrated parameter
            # lambda_8 = 0.1,  # Non-calibrated parameter
            prev.adeno = 0.25,
            prop.adeno.sm = 0.0464/(0.0464 + 0.0764)# From Wu et al., 2006
  ) 
  #--- Upper bounds ---#
  v.ub <- c(lambda_1 = 2e-05,
            g        = 4,
            lambda_2 = 0.10,
            lambda_3 = 0.04,
            lambda_4 = 0.5,
            lambda_5 = 0.3,
            lambda_6 = 0.7,
            # lambda_7 = 0.06, # Non-calibrated parameter
            # lambda_8 = 0.4,  # Non-calibrated parameter
            prev.adeno = 0.35,
            prop.adeno.sm = 0.1856/(0.1856 + 0.0096) # From Wu et al., 2006
  ) 
  #--- Standard errors based on bounds ---#
  v.se <- (v.ub - v.lb)/(2*1.96)
  
  #--- Return output ---#
  return(out = list(v.lb = v.lb,
                    v.ub = v.ub,
                    v.se = v.se))
}

sample.prior.informed <- function(n) { # n <- 100
  ### Arguments:
  ##    - n: the number of samples desired 
  ### Depends on:
  ##    - `lhs` package
  ##    - `gen_params_init()` function that returns an initial parameter set
  ##    - `gen_bounds()` function that returns lower and upper bounds, and 
  ##      standard errors
  
  ### Obtain parameter initial values
  v.params.init <- gen_params_init()
  
  ### Obtain bounds
  l.bounds <- gen_bounds()
  
  ### Number of parameters
  n.params <- length(l.bounds$v.se)
  
  ### Draw LHS from Uniform[0,1] distributions
  lhs.grid <- randomLHS(n = n, k = n.params) 
  colnames(lhs.grid) <- v.names.params
  
  ### Transformed design 
  ## Transformed bounds
  lb.transf <- log(l.bounds$v.lb)
  ub.transf <- log(l.bounds$v.ub)
  
  # Find means and SDs of logit-normal and log-normal based on bounds 
  # assuming bounds are represent the 95% equal tailed interval for these 
  # distributions
  mu.transf <- (ub.transf + lb.transf)/2
  sd.transf <- (ub.transf - lb.transf)/(2*1.96)
  
  ### Log-Normal the rate parameters
  ## Get values in Normal scale
  normal.grid <- lhs.grid
  ## Apply Normal CDF
  for (i in 3:(n.params-2)){ 
    normal.grid[, i] <- qnorm(lhs.grid[,i], mu.transf[i], sd.transf[i])
  }
  ## Get values in Original scale
  transf.grid <- exp(normal.grid)
  
  ### Beta distributions to the proportion parameters
  ## Proportion of small adenomas
  prop.adeno.sm.params <- beta_params(mean = v.params.init["prop.adeno.sm"], 
                                      sigma = l.bounds$v.se["prop.adeno.sm"])
  transf.grid[, "prop.adeno.sm"] <- qbeta(lhs.grid[, "prop.adeno.sm"], 
                                          shape1 = prop.adeno.sm.params$alpha, 
                                          shape2 = prop.adeno.sm.params$beta) 
  
  library(mvtnorm)
  #Distribution for shape and scale of Weibull comes from uncertainty on the rates
  #These variables come from EVSI function script
  transf.grid[,c("lambda_1","g")]<-rmvnorm(n,lamb_g_mean,lamb_g_cov)^lamb_g_rescale
  
  ## Prevalance of adenoma
  #Prevalance is the proportion of patients who have had the event at age 50
  #This comes directly from the survival curves with are associated with uncertainty.
  transf.grid[, "prev.adeno"] <- pweibull(50,shape=transf.grid[,"g"],
                                          scale=transf.grid[,"lambda_1"]^(-1/transf.grid[,"g"]))
  
  ## Return transformed LHS grid
  return(transf.grid) 
}

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

#Simulate PSA distribution from a multivariate normal distribution 
#We use this rather than the parameters from the model fit so the Bayesian model can match the 
#prior exactly.
PSA.coefs.scale<-rmvnorm(n.sim,apply(coefs^(1/15),2,mean),cov(coefs^(1/15),use="complete.obs"))
#Rescale to natural scale
PSA.coefs<-PSA.coefs.scale^15

lamb_g_mean<-apply(coefs^(1/15),2,mean)
lamb_g_cov<-cov(coefs^(1/15),use="complete.obs")
lamb_g_rescale<-15

gen_psa <- function(seed = 7783){
  load("04_crc-nhm-det_posterior-IMIS.rData")
  #n.sim <- nrow(post_imis)
  set.seed = seed
  df.psa.params <- data.frame(
    ### NHM parameters
    post_imis,
    ### Colonoscopy test characteristics
    ## Samll Adenomas
    sens.COL = rbeta(n.sim, (269.4+104.1), (94.6+15.6)), # Knudsen 2016 and originally from  Van Rijn 2006
    spec.COL = rbeta(n.sim, 2475, 378), # Knudsen 2016 and originally from Schroy 2013
    ## Large Adenomas & CRC
    sens.COL.CRC = rbeta(n.sim, 59.3, 1.2), # Knudsen 2016 and originally from  Van Rijn 2006
    spec.COL.CRC = rbeta(n.sim, 2475, 378), # Knudsen 2016 and originally from Schroy 2013
    ### Increase rate from Normal to Small Adenoma whith history of lesion
    ## Hazard ratio (HR) for Low risk (LR)
    hr.Hist.LR = exp(rnorm(n.sim, mean = log(2), sd = (log(3) - log(1))/(2*1.96))), # Expert opinion
    ## Hazard ratio (HR) for High risk (HR)
    hr.Hist.HR = exp(rnorm(n.sim, mean = log(3), sd = (log(4) - log(2))/(2*1.96))), # Expert opinion
    ### Utilities
    ## Normal
    u.Normal = 1-exp(rnorm(n.sim, mean = log(0.01), sd = (log(0.02) - log(0.0001))/(2*1.96))),
    ## Clinical CRC Early
    u.CRC_Early = 1-exp(rnorm(n.sim, mean = log(1-0.855), sd = (log(0.3) - log(0.1))/(2*1.96))),
    ## Clinical CRC Late
    u.CRC_Late = 1-exp(rnorm(n.sim, mean = log(1-0.3), sd = (log(0.6) - log(0.4))/(2*1.96))),
    ### Cost parameters
    ## Colonoscopy
    c.COL = exp(rnorm(n.sim, mean = log(10000), sd = (log(11000) - log(9000))/(2*1.96))),
    ## Cancer by stages
    c.crc.early = exp(rnorm(n.sim, mean = log(21524), sd = (log(23000) - log(20000))/(2*1.96)))*1.2356, # Inlfation factor: 1.2356
    c.crc.late = exp(rnorm(n.sim, mean = log(37000), sd = (log(39000) - log(35000))/(2*1.96)))*1.2356 # Inlfation factor: 1.2356
  )
  return(df.psa.params)
}

set.seed(12)
### Sample from informed priors
df.prior.informed <- sample.prior.informed(n.sim)

### Substitute posterior of calibrated parameters with MAP estimates
df.psa.params <- gen_psa()
df.psa.params <- df.psa.params[1:n.sim, ]
df.psa.params.informed <- df.psa.params
df.psa.params.informed[, 1:length(v.params0)] <- df.prior.informed

### Functions to run health economic model
crc_screening <- function(v.params.scr, SCR = T){
  ### Source of Model: https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-6-136
  ### Arguments:  
  #     v.params: vector of deep-model parameters 
  #
  ### Uses external functions:
  # crc_nhm_tp: Computes a transition probability array from a continuous time model
  #
  with(as.list(v.params.scr), {
    
    #### Initial State vector ####
    v.init.scr <- c(1-(prev.adeno + prev.precl.early + prev.precl.late), # Normal
                    0,                              # Normal with HISTORY LR
                    0,                              # Normal with HISTORY HR
                    prev.adeno*prop.adeno.sm,       # SMALL Adenoma
                    0,                              # SMALL Adenoma with History LR
                    0,                              # SMALL Adenoma with History HR
                    prev.adeno*(1 - prop.adeno.sm), # LARGE Adenoma
                    0,                              # LARGE Adenoma with History LR
                    0,                              # LARGE Adenoma with History HR
                    prev.precl.early,               # Preclinical EARLY
                    0,                              # Preclinical EARLY with History LR
                    0,                              # Preclinical EARLY with History HR
                    prev.precl.late,                # Preclinical LATE
                    0,                              # Preclinical LATE with History LR
                    0,                              # Preclinical LATE with History HR
                    0, 0, 0, 0)
    names(v.init.scr) <- v.name.states.scr
    
    ### Obtain transition probability array of screening model
    a.P.scr <- crc_screening_tp(v.params.scr, SCR = SCR)
    
    ### Initialize Markov Trace
    m.M <- matrix(0, ncol = n.s.scr, nrow = (n.t+1))
    colnames(m.M) <- v.name.states.scr
    rownames(m.M) <- paste("Cycle", 0:(n.t), sep = "")
    ### Initialize Markov Incidence Array
    a.A <- array(0, dim = c(n.s.scr, n.s.scr, n.t + 1), 
                 dimnames = list(v.name.states.scr, v.name.states.scr, 
                                 paste("Cycle", 0:(n.t), sep = "")))
    ### Set starting distribution of population
    m.M[1, ] <- v.init.scr
    diag(a.A[, , 1]) <- v.init.scr
    
    ### Run Markov model
    for (t in 1:n.t) { # t <- 2
      m.M[t + 1, ] <- m.M[t, ] %*% a.P.scr[, , t]
      a.A[, , t + 1] <- m.M[t, ] * a.P.scr[, , t]
    }
    
    #--------------------------------------------#
    #### Generate Epidemiologic model outputs ####
    #--------------------------------------------#
    ### Population at risk
    v.n.atrisk <- rowSums(m.M[, c("Normal",
                                  "Normal_Hist_LR",
                                  "Normal_Hist_HR",
                                  "SmallAdeno",
                                  "SmallAdeno_LR",
                                  "SmallAdeno_HR",
                                  "LargeAdeno",
                                  "LargeAdeno_LR",
                                  "LargeAdeno_HR",
                                  "PreClinCRC_Early",
                                  "PreClinCRC_Early_LR",
                                  "PreClinCRC_Early_HR",
                                  "PreClinCRC_Late",
                                  "PreClinCRC_Late_LR",
                                  "PreClinCRC_Late_HR")]
    )
    ### Compute Prevalence of Adenoma across all ages
    prevAdeno <- rowSums(m.M[, c("SmallAdeno", 
                                 "SmallAdeno_LR",
                                 "SmallAdeno_HR",
                                 "LargeAdeno",
                                 "LargeAdeno_LR",
                                 "LargeAdeno_HR")])/v.n.atrisk
    names(prevAdeno) <- v.ages
    
    ### Compute Proportion of SMALL Adenomas across all ages
    propAdenoSmall <- rowSums(m.M[, c("SmallAdeno", 
                                      "SmallAdeno_LR",
                                      "SmallAdeno_HR")])/rowSums(m.M[, c("SmallAdeno", 
                                                                         "SmallAdeno_LR",
                                                                         "SmallAdeno_HR",
                                                                         "LargeAdeno",
                                                                         "LargeAdeno_LR",
                                                                         "LargeAdeno_HR")])
    names(propAdenoSmall) <- v.ages
    
    ### Transitions to Detected Cancer for ALL Stages by age
    transDetCRC <- colSums(colSums(a.A[c("Normal",
                                         "Normal_Hist_LR",
                                         "Normal_Hist_HR",
                                         "SmallAdeno",
                                         "SmallAdeno_LR",
                                         "SmallAdeno_HR",
                                         "LargeAdeno",
                                         "LargeAdeno_LR",
                                         "LargeAdeno_HR",
                                         "PreClinCRC_Early",
                                         "PreClinCRC_Early_LR",
                                         "PreClinCRC_Early_HR",
                                         "PreClinCRC_Late",
                                         "PreClinCRC_Late_LR",
                                         "PreClinCRC_Late_HR"), # From States
                                       c("ClinCRC_Early",
                                         "ClinCRC_Late"), # To States
                                       -1]))
    ### Transitions to Detected Cancer for EARLY Stages by age
    transDetCRC_early <- colSums(a.A[c("Normal",
                                       "Normal_Hist_LR",
                                       "Normal_Hist_HR",
                                       "SmallAdeno",
                                       "SmallAdeno_LR",
                                       "SmallAdeno_HR",
                                       "LargeAdeno",
                                       "LargeAdeno_LR",
                                       "LargeAdeno_HR",
                                       "PreClinCRC_Early",
                                       "PreClinCRC_Early_LR",
                                       "PreClinCRC_Early_HR",
                                       "PreClinCRC_Late",
                                       "PreClinCRC_Late_LR",
                                       "PreClinCRC_Late_HR"), # From States
                                     c("ClinCRC_Early"), # To States
                                     -1])
    ### Transitions to Detected Cancer for LATE Stages by age
    transDetCRC_late <- colSums(a.A[c("Normal",
                                      "Normal_Hist_LR",
                                      "Normal_Hist_HR",
                                      "SmallAdeno",
                                      "SmallAdeno_LR",
                                      "SmallAdeno_HR",
                                      "LargeAdeno",
                                      "LargeAdeno_LR",
                                      "LargeAdeno_HR",
                                      "PreClinCRC_Early",
                                      "PreClinCRC_Early_LR",
                                      "PreClinCRC_Early_HR",
                                      "PreClinCRC_Late",
                                      "PreClinCRC_Late_LR",
                                      "PreClinCRC_Late_HR"), # From States
                                    c("ClinCRC_Late"), # To States
                                    -1]) 
    
    ### Lifetime risk 
    ltRisk <- 1-exp(-sum(-log(1-transDetCRC))) # sum(transDetCa)
    
    ### Transitions to CRC Death for from ALL health states by age
    transDeathCRC <- colSums(a.A[c("Normal",
                                   "Normal_Hist_LR",
                                   "Normal_Hist_HR",
                                   "SmallAdeno",
                                   "SmallAdeno_LR",
                                   "SmallAdeno_HR",
                                   "LargeAdeno",
                                   "LargeAdeno_LR",
                                   "LargeAdeno_HR",
                                   "PreClinCRC_Early",
                                   "PreClinCRC_Early_LR",
                                   "PreClinCRC_Early_HR",
                                   "PreClinCRC_Late",
                                   "PreClinCRC_Late_LR",
                                   "PreClinCRC_Late_HR"), # From States
                                 c("CRC_Death"), # To States
                                 -1])
    
    ### Compute Incidence across all ages
    ## ALL Stages
    incCRC <- -log(1-transDetCRC/v.n.atrisk[-(n.t+1)])
    ## Early Stage
    incCRC_early <- -log(1-transDetCRC_early/v.n.atrisk[-(n.t+1)])
    ## Late Stage
    incCRC_late <- -log(1-transDetCRC_late/v.n.atrisk[-(n.t+1)])
    ## Name targets
    names(incCRC) <- names(incCRC_early) <- names(incCRC_late) <- v.ages
    
    ### Compute CRC mortality across all ages
    deathCRC <- -log(1-transDeathCRC/v.n.atrisk[-(n.t+1)])
    names(deathCRC) <- v.ages
    
    #-------------------------------------------------#
    #### Generate Cost-Effectiveness model outputs ####
    #-------------------------------------------------#
    ### Discounting vectors
    ## Effectiveness
    v.dwe <- 1 / ((1 + d.e) ^ (0:(n.t)))
    ## Costs
    v.dwc <- 1 / ((1 + d.c) ^ (0:(n.t)))
    
    ### Expected Life Years (LY)
    v.ly <- c(rep(1, (n.s.scr-2)), 0, 0)
    v.ev.ly <- m.M %*% v.ly
    ev.ly <- sum(v.ev.ly)
    ev.ly.disc <- sum(v.ev.ly * v.dwe)
    
    ### Expected Quality-Adjusted Life Years (QALY)
    ## Vector of utilities
    v.u <- c(rep(u.Normal, 15),
             u.CRC_Early,
             u.CRC_Late,
             0, 0)
    names(v.u) <- v.name.states.scr
    v.ev.qaly <- m.M %*% v.u
    ev.qaly <- sum(v.ev.qaly)
    ev.qaly.disc <- sum(v.ev.qaly * v.dwe)
    
    ### Expected Costs
    ## Matrix with state costs over time
    m.costs <- cbind(
      ## Normal
      (c(p.Screen, 0) * SCR) * c.COL,
      (c(p.Surv.LR, 0) * SCR) * c.COL,
      (c(p.Surv.HR, 0) * SCR) * c.COL,
      ## Small Adenoma
      (c(p.Screen, 0) * SCR) * c.COL,
      (c(p.Surv.LR, 0) * SCR) * c.COL,
      (c(p.Surv.HR, 0) * SCR) * c.COL,
      ## Large Adenoma
      (c(p.Screen, 0) * SCR) * c.COL,
      (c(p.Surv.LR, 0) * SCR) * c.COL,
      (c(p.Surv.HR, 0) * SCR) * c.COL,
      ## PreClinCRC_Early
      (c(p.Screen, 0) * SCR) * c.COL,
      (c(p.Surv.LR, 0) * SCR) * c.COL,
      (c(p.Surv.HR, 0) * SCR) * c.COL,
      ## PreClinCRC_Late
      (c(p.Screen, 0) * SCR) * c.COL,
      (c(p.Surv.LR, 0) * SCR) * c.COL,
      (c(p.Surv.HR, 0) * SCR) * c.COL,
      ## ClinCRC_Early
      c.crc.early, 
      ## ClinCRC_Late 
      c.crc.late,  
      0, 0)
    colnames(m.costs) <- v.name.states.scr
    v.ev.costs <- rowSums(m.M * m.costs)
    ev.costs <- sum(v.ev.costs)
    ev.costs.disc <- sum(v.ev.costs * v.dwc)
    
    #----------------------------#
    #### Return model outputs ####
    #----------------------------#
    
    return(list(m.M = m.M,
                a.A = a.A,
                ### Epidemiologic outcomes
                prevAdeno    = prevAdeno[as.character(v.ages.prev)], # Select only specific ages to macth target ages
                propAdenoSmall = propAdenoSmall[as.character(v.ages.prev)],
                incCRC       = incCRC[as.character(v.ages.crc)],
                incCRC_early = incCRC_early[as.character(v.ages.crc)],
                incCRC_late  = incCRC_late[as.character(v.ages.crc)],
                ltRisk       = ltRisk,
                deathCRC     = deathCRC[as.character(v.ages.crc)],
                ### CEA outputs
                v.ev.ly       = v.ev.ly,
                ev.ly         = ev.ly,
                ev.ly.disc    = ev.ly.disc,
                v.ev.qaly     = v.ev.qaly,
                ev.qaly       = ev.qaly,
                ev.qaly.disc  = ev.qaly.disc,
                v.ev.costs    = v.ev.costs,
                ev.costs      = ev.costs,
                ev.costs.disc = ev.costs.disc))
  }
  )
}

crc_screening_tp <- function(v.params.scr, SCR = T){
  ### Source of Model: https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-6-136
  ### Depends on:
  ###   - 
  with(as.list(v.params.scr), {
    ### Adjust screening and surveillance parameters based on screnning YES or NO
    ## Screening
    p.Screen  <- p.Screen * SCR
    ## Surveillance on Low Risk (LR)
    p.Surv.LR <- p.Surv.LR * SCR
    ## Surveillance on High Risk (LR)
    p.Surv.HR <- p.Surv.HR * SCR
    
    ### Adenoma testing results
    p.TruePos  <- sens.COL
    p.FalseNeg <- 1 - sens.COL
    p.TrueNeg  <- spec.COL
    p.FalsePos <- 1 - spec.COL 
    
    ### CRC testing results
    p.TruePosCRC  <- sens.COL.CRC
    p.FalseNegCRC <- 1 - sens.COL.CRC
    p.TrueNegCRC  <- spec.COL.CRC
    p.FalsePosCRC <- 1 - spec.COL.CRC 
    
    ### Obtain transition probability array of NHM
    a.P <- crc_nhm_tp(v.params.scr)
    
    #### Transition Probability Array of NHM with History ####
    a.P.hist.LR <- crc_nhm_hist_tp(v.params.scr, hist.LR = T)
    a.P.hist.HR <- crc_nhm_hist_tp(v.params.scr, hist.LR = F)
    
    #### Transition Probability Array with Screening ####
    a.P.scr <- array(0, dim = list(n.s.scr, n.s.scr, n.t),
                     dimnames = list(v.name.states.scr, v.name.states.scr, v.ages))
    ### Vectors with indexes of states for history paths
    ## NO History
    v.hist.NO <- c(1, 4, 7, 10, 13, 16:19)
    ## Low Risk
    v.hist.LR <- c(2, 5, 8, 11, 14, 16:19)
    ## Highs Risk
    v.hist.HR <- c(3, 6, 9, 12, 15, 16:19)
    #### Fill Transition Intensity Array
    ### From Normal (1)
    a.P.scr["Normal", v.hist.NO, ] <- t((1-p.Screen) * t(a.P["Normal", , ]) +             # If not screeened -> Follow NHM
                                          p.Screen * (p.TrueNeg * t(a.P["Normal", , ])) + # If screened and true negative -> Follow NHM
                                          p.Screen * (p.FalsePos * t(a.P["Normal", , ]))  # If screened and false positive -> Follow NHM 
    )
    # ## -> Normal
    # a.P.scr["Normal", "Normal", ] <- ((1-p.Screen) * a.P["Normal", "Normal", ] +            # If not screeened -> Follow NHM
    #                              p.Screen * (p.TrueNeg * a.P["Normal", "Normal", ]) + # If screened and true negative -> Follow NHM
    #                              p.Screen * (p.FalsePos * a.P["Normal", "Normal", ])) # If screened and false positive -> Follow NHM 
    # ## -> Normal_Hist
    # # a.P.scr["Normal", "Normal_Hist", ] <- ((1-p.Screen) * a.P["Normal", "Normal", ] +            # If not screeened -> Follow NHM
    # #                                     p.Screen * (p.TrueNeg * a.P["Normal", "Normal", ]) + # If screened and true negative -> Follow NHM
    # #                                     p.Screen * (p.FalsePos * a.P["Normal", "Normal", ])) # If screened and false positive -> Follow NHM 
    # ## -> SmallAdeno
    # a.P.scr["Normal", "SmallAdeno", ] <- ((1-p.Screen) * a.P["Normal", "SmallAdeno", ] +            # If not screeened -> Follow NHM
    #                                     p.Screen * (p.TrueNeg * a.P["Normal", "SmallAdeno", ]) + # If screened and true negative -> Follow NHM
    #                                     p.Screen * (p.FalsePos * a.P["Normal", "SmallAdeno", ])) # If screened and false positive -> Follow NHM 
    # ## -> LargeAdeno
    # a.P.scr["Normal", "LargeAdeno", ] <- ((1-p.Screen) * a.P["Normal", "LargeAdeno", ] +            # If not screeened -> Follow NHM
    #                                     p.Screen * (p.TrueNeg * a.P["Normal", "LargeAdeno", ]) + # If screened and true negative -> Follow NHM
    #                                     p.Screen * (p.FalsePos * a.P["Normal", "LargeAdeno", ])) # If screened and false positive -> Follow NHM 
    # ## -> PreClinCRC_Early
    # a.P.scr["Normal", "PreClinCRC_Early", ] <- ((1-p.Screen) * a.P["Normal", "PreClinCRC_Early", ] +            # If not screeened -> Follow NHM
    #                                         p.Screen * (p.TrueNeg * a.P["Normal", "PreClinCRC_Early", ]) + # If screened and true negative -> Follow NHM
    #                                         p.Screen * (p.FalsePos * a.P["Normal", "PreClinCRC_Early", ])) # If screened and false positive -> Follow NHM 
    # ## -> ClinCRC_Late
    # a.P.scr["Normal", "PreClinCRC_Late", ] <- ((1-p.Screen) * a.P["Normal", "PreClinCRC_Late", ] +            # If not screeened -> Follow NHM
    #                                               p.Screen * (p.TrueNeg * a.P["Normal", "PreClinCRC_Late", ]) + # If screened and true negative -> Follow NHM
    #                                               p.Screen * (p.FalsePos * a.P["Normal", "PreClinCRC_Late", ])) # If screened and false positive -> Follow NHM 
    # ## -> ClinCRC_Early
    # a.P.scr["Normal", "ClinCRC_Early", ] <- ((1-p.Screen) * a.P["Normal", "ClinCRC_Early", ] +            # If not screeened -> Follow NHM
    #                                               p.Screen * (p.TrueNeg * a.P["Normal", "ClinCRC_Early", ]) + # If screened and true negative -> Follow NHM
    #                                               p.Screen * (p.FalsePos * a.P["Normal", "ClinCRC_Early", ])) # If screened and false positive -> Follow NHM 
    # ## -> PreClinCRC_Late
    # a.P.scr["Normal", "ClinCRC_Late", ] <- ((1-p.Screen) * a.P["Normal", "ClinCRC_Late", ] +            # If not screeened -> Follow NHM
    #                                              p.Screen * (p.TrueNeg * a.P["Normal", "ClinCRC_Late", ]) + # If screened and true negative -> Follow NHM
    #                                              p.Screen * (p.FalsePos * a.P["Normal", "ClinCRC_Late", ])) # If screened and false positive -> Follow NHM
    # ## -> CRC_Death
    # a.P.scr["Normal", "CRC_Death", ] <- ((1-p.Screen) * a.P["Normal", "CRC_Death", ] +            # If not screeened -> Follow NHM
    #                                           p.Screen * (p.TrueNeg * a.P["Normal", "CRC_Death", ]) + # If screened and true negative -> Follow NHM
    #                                           p.Screen * (p.FalsePos * a.P["Normal", "CRC_Death", ])) # If screened and false positive -> Follow NHM
    # ## -> OC_Death
    # a.P.scr["Normal", "OC_Death", ] <- ((1-p.Screen) * a.P["Normal", "OC_Death", ] +            # If not screeened -> Follow NHM
    #                                           p.Screen * (p.TrueNeg * a.P["Normal", "OC_Death", ]) + # If screened and true negative -> Follow NHM
    #                                           p.Screen * (p.FalsePos * a.P["Normal", "OC_Death", ])) # If screened and false positive -> Follow NHM
    ### From Normal_Hist_LR (2)
    a.P.scr["Normal_Hist_LR", v.hist.LR, ] <- t((1-p.Surv.LR) * t(a.P.hist.LR["Normal_Hist", -1, ]) +             # If not screeened -> Follow NHM
                                                  p.Surv.LR * (p.TrueNeg * t(a.P.hist.LR["Normal_Hist", -1, ])) + # If screened and true negative -> Follow NHM
                                                  p.Surv.LR * (p.FalsePos * t(a.P.hist.LR["Normal_Hist", -1, ]))  # If screened and false positive -> Follow NHM 
    )
    ### From Normal_Hist_HR (3)
    a.P.scr["Normal_Hist_HR", v.hist.HR, ] <- t((1-p.Surv.HR) * t(a.P.hist.HR["Normal_Hist", -1, ]) +             # If not screeened -> Follow NHM
                                                  p.Surv.HR * (p.TrueNeg * t(a.P.hist.HR["Normal_Hist", -1, ])) + # If screened and true negative -> Follow NHM
                                                  p.Surv.HR * (p.FalsePos * t(a.P.hist.HR["Normal_Hist", -1, ]))  # If screened and false positive -> Follow NHM 
    )
    
    ### From SmallAdeno (4)
    a.P.scr["SmallAdeno", v.hist.NO, ] <- t((1-p.Screen) * t(a.P["SmallAdeno", , ]) +  # If not screeened -> Follow NHM
                                              p.Screen * (p.FalseNeg * t(a.P["SmallAdeno", , ])) # If screened and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["SmallAdeno", "Normal_Hist_LR", ] <- p.Screen * p.TruePos # If screened and true positive -> Normal_Hist
    
    ### From SmallAdeno_LR (5)
    a.P.scr["SmallAdeno_LR", v.hist.LR, ] <- t((1-p.Surv.LR) * t(a.P.hist.LR["SmallAdeno", -1, ]) +  # If not surveillance -> Follow NHM
                                                 p.Surv.LR * (p.FalseNeg * t(a.P.hist.LR["SmallAdeno", -1, ])) # If surveillance and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["SmallAdeno_LR", "Normal_Hist_LR", ] <- p.Surv.LR * p.TruePos # If screened and true positive -> Normal_Hist
    
    ### From SmallAdeno_HR (6)
    a.P.scr["SmallAdeno_HR", v.hist.LR, ] <- t((1-p.Surv.HR) * t(a.P.hist.HR["SmallAdeno", -1, ]) +  # If not surveillance -> Follow NHM
                                                 p.Surv.HR * (p.FalseNeg * t(a.P.hist.HR["SmallAdeno", -1, ])) # If surveillance and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["SmallAdeno_HR", "Normal_Hist_LR", ] <- p.Surv.HR * p.TruePos # If screened and true positive -> Normal_Hist
    
    ### From LargeAdeno (7)     
    a.P.scr["LargeAdeno", v.hist.NO, ] <- t((1-p.Screen) * t(a.P["LargeAdeno", , ]) +            # If not screeened -> Follow NHM
                                              p.Screen * (p.FalseNegCRC * t(a.P["LargeAdeno", , ])) # If screened and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["LargeAdeno", "Normal_Hist_HR", ] <- p.Screen * p.TruePosCRC # If screened and true positive -> Normal_Hist
    
    ### From LargeAdeno_LR (8)     
    a.P.scr["LargeAdeno_LR", v.hist.LR, ] <- t((1-p.Surv.LR) * t(a.P.hist.LR["LargeAdeno", -1, ]) +            # If not surveillance -> Follow NHM
                                                 p.Surv.LR * (p.FalseNegCRC * t(a.P.hist.LR["LargeAdeno", -1, ])) # If surveillane and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["LargeAdeno_LR", "Normal_Hist_HR", ] <- p.Surv.LR * p.TruePosCRC # If screened and true positive -> Normal_Hist
    
    ### From LargeAdeno_HR (9)     
    a.P.scr["LargeAdeno_HR", v.hist.HR, ] <- t((1-p.Surv.HR) * t(a.P.hist.HR["LargeAdeno", -1, ]) +            # If not surveillance -> Follow NHM
                                                 p.Surv.HR * (p.FalseNegCRC * t(a.P.hist.HR["LargeAdeno", -1, ])) # If surveillane and false negative -> Follow NHM
    )
    ## -> Normal_Hist
    a.P.scr["LargeAdeno_HR", "Normal_Hist_HR", ] <- p.Surv.HR * p.TruePosCRC # If screened and true positive -> Normal_Hist
    
    ### From PreClinCRC_Early (10) 
    a.P.scr["PreClinCRC_Early", v.hist.NO, ] <- t((1-p.Screen) * t(a.P["PreClinCRC_Early", , ]) + # If not screeened -> Follow NHM
                                                    p.Screen * (p.FalseNegCRC * t(a.P["PreClinCRC_Early", , ])) # If screened and false negative -> Follow NHM
                                                  
    )
    a.P.scr["PreClinCRC_Early", "ClinCRC_Early", ]  <- a.P.scr["PreClinCRC_Early", "ClinCRC_Early", ] + p.Screen * p.TruePosCRC
    
    ### From PreClinCRC_Early_LR (11) 
    a.P.scr["PreClinCRC_Early_LR", v.hist.LR, ] <- t((1-p.Surv.LR) * t(a.P.hist.LR["PreClinCRC_Early", -1, ]) + # If not screeened -> Follow NHM
                                                       p.Surv.LR * (p.FalseNegCRC * t(a.P.hist.LR["PreClinCRC_Early", -1, ])) # If screened and false negative -> Follow NHM
                                                     
    )
    a.P.scr["PreClinCRC_Early_LR", "ClinCRC_Early", ]  <- a.P.scr["PreClinCRC_Early_LR", "ClinCRC_Early", ] + p.Surv.LR * p.TruePosCRC
    
    ### From PreClinCRC_Early_HR (12) 
    a.P.scr["PreClinCRC_Early_HR", v.hist.HR, ] <- t((1-p.Surv.HR) * t(a.P.hist.HR["PreClinCRC_Early", -1, ]) + # If not screeened -> Follow NHM
                                                       p.Surv.HR * (p.FalseNegCRC * t(a.P.hist.HR["PreClinCRC_Early", -1, ])) # If screened and false negative -> Follow NHM
                                                     
    )
    a.P.scr["PreClinCRC_Early_HR", "ClinCRC_Early", ]  <- a.P.scr["PreClinCRC_Early_HR", "ClinCRC_Early", ] + p.Surv.HR * p.TruePosCRC
    
    ### From PreClinCRC_Late (13)
    a.P.scr["PreClinCRC_Late", v.hist.NO, ] <- t((1-p.Screen) * t(a.P["PreClinCRC_Late", , ]) + # If not screeened -> Follow NHM
                                                   p.Screen * (p.FalseNegCRC * t(a.P["PreClinCRC_Late", , ])) # If screened and false negative -> Follow NHM
                                                 
    )
    a.P.scr["PreClinCRC_Late", "ClinCRC_Late", ]  <- a.P.scr["PreClinCRC_Late", "ClinCRC_Late", ] + p.Screen * p.TruePosCRC
    
    ### From PreClinCRC_Late_LR (14)
    a.P.scr["PreClinCRC_Late_LR", v.hist.LR, ] <- t((1-p.Surv.LR) * t(a.P.hist.LR["PreClinCRC_Late", -1, ]) + # If not screeened -> Follow NHM
                                                      p.Surv.LR * (p.FalseNegCRC * t(a.P.hist.LR["PreClinCRC_Late", -1, ])) # If screened and false negative -> Follow NHM
                                                    
    )
    a.P.scr["PreClinCRC_Late_LR", "ClinCRC_Late", ]  <- a.P.scr["PreClinCRC_Late_LR", "ClinCRC_Late", ] + p.Surv.LR * p.TruePosCRC
    
    ### From PreClinCRC_Late_HR (15)
    a.P.scr["PreClinCRC_Late_HR", v.hist.LR, ] <- t((1-p.Surv.HR) * t(a.P.hist.HR["PreClinCRC_Late", -1, ]) + # If not screeened -> Follow NHM
                                                      p.Surv.HR * (p.FalseNegCRC * t(a.P.hist.HR["PreClinCRC_Late", -1, ])) # If screened and false negative -> Follow NHM
                                                    
    )
    a.P.scr["PreClinCRC_Late_HR", "ClinCRC_Late", ]  <- a.P.scr["PreClinCRC_Late_HR", "ClinCRC_Late", ] + p.Surv.HR * p.TruePosCRC
    
    ### From ClinCRC_Early (16)
    a.P.scr["ClinCRC_Early", v.hist.NO, ] <- a.P["ClinCRC_Early", , ]
    
    
    ### From ClinCRC_Late (17)
    a.P.scr["ClinCRC_Late", v.hist.NO, ]  <- a.P["ClinCRC_Late", , ]
    
    ### From CRC_Death
    a.P.scr["CRC_Death", "CRC_Death", ] <- 1
    
    ### From OC_Death
    a.P.scr["OC_Death", "OC_Death", ] <- 1
    
    
    ### Check if Transition Probability array is valid
    valid <- apply(a.P.scr, 3, function(x) all.equal(sum(rowSums(x)), n.s.scr))
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n.t)))) {
      stop("This is not a valid transition Matrix")
    }
    
    return(a.P.scr) # Return transition probability array
  }
  )
}

crc_nhm_tp <- function(v.params){
  ### Source of Model: https://bmccancer.biomedcentral.com/articles/10.1186/1471-2407-6-136
  ### Depends on:
  ###   - Package `msm` to compute the matrix exponential
  with(as.list(v.params), {
    
    #### State names and cycles ####
    v.name.states <- c("Normal",
                       "SmallAdeno",
                       "LargeAdeno",
                       "PreClinCRC_Early",
                       "PreClinCRC_Late",
                       "ClinCRC_Early",
                       "ClinCRC_Late",
                       "CRC_Death",
                       "OC_Death") 
    n.s <- length(v.name.states) # Number of states
    
    #### Continuous time model ####
    ### Transition Intensity Array
    a.G <- array(0, dim = list(n.s, n.s, n.t), 
                 dimnames = list(v.name.states, v.name.states, v.ages))
    
    ## Age-dependent rate
    lambda_1_t <- lambda_1*g*v.ages^(g-1)
    
    ## Fill Transition Intensity Array
    a.G["Normal",     "Normal", ]                 <- -(lambda_1_t + lambda_asr)
    a.G["Normal",     "SmallAdeno", ]             <- lambda_1_t
    a.G["Normal",     "OC_Death", ]               <- lambda_asr
    a.G["SmallAdeno", "SmallAdeno", ]             <- -(lambda_2 + lambda_asr)
    a.G["SmallAdeno", "LargeAdeno", ]             <- lambda_2  
    a.G["SmallAdeno",     "OC_Death", ]           <- lambda_asr    
    a.G["LargeAdeno", "LargeAdeno", ]             <- -(lambda_3 + lambda_asr)
    a.G["LargeAdeno", "PreClinCRC_Early", ]       <- lambda_3 
    a.G["LargeAdeno",     "OC_Death", ]           <- lambda_asr    
    a.G["PreClinCRC_Early", "PreClinCRC_Early", ] <- -(lambda_4 + lambda_5 + lambda_asr)
    a.G["PreClinCRC_Early", "PreClinCRC_Late", ]  <- lambda_4
    a.G["PreClinCRC_Early", "ClinCRC_Early", ]    <- lambda_5
    a.G["PreClinCRC_Early",     "OC_Death", ]     <- lambda_asr    
    a.G["PreClinCRC_Late", "PreClinCRC_Late", ]   <- -(lambda_6 + lambda_asr)
    a.G["PreClinCRC_Late", "ClinCRC_Late", ]      <- lambda_6
    a.G["PreClinCRC_Late",     "OC_Death", ]      <- lambda_asr
    a.G["ClinCRC_Early", "ClinCRC_Early", ]       <- -(lambda_7 + lambda_asr)
    a.G["ClinCRC_Early", "CRC_Death", ]           <- lambda_7
    a.G["ClinCRC_Early",     "OC_Death", ]        <- lambda_asr
    a.G["ClinCRC_Late", "ClinCRC_Late", ]         <- -(lambda_8 + lambda_asr)
    a.G["ClinCRC_Late", "CRC_Death", ]            <- lambda_8
    a.G["ClinCRC_Late",     "OC_Death", ]         <- lambda_asr
    
    ### Compute Transition proability 
    a.P <- sapply(seq_along(v.ages), 
                  function(x) MatrixExp(a.G[,,x]), 
                  simplify="array")
    
    ### Check if Transition Probability array is valid
    valid <- apply(a.P, 3, function(x) all.equal(sum(rowSums(x)), n.s))
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n.t)))) {
      stop("This is not a valid transition Matrix")
    }
    
    return(a.P) # Return transition probability array
  }
  )
}

crc_nhm_hist_tp <- function(v.params.scr, hist.LR = T){
  ### Computes a transition probability array accounting for history
  ### Inputs:
  ###   - v.params.scr: Vector of parameters of screening model
  ###   - hist.LR: Indicator variable selecting history for low risk (=T) or high risk (=F)
  ### Depends on:
  ###   - Package `msm` to compute the matrix exponential
  with(as.list(v.params.scr), {
    
    hr.Hist <- hr.Hist.LR*(hist.LR) + hr.Hist.HR*(1-hist.LR)
    
    #### Continuous time model ####
    v.name.states.hist <- c("Normal",
                            "Normal_Hist",
                            "SmallAdeno",
                            "LargeAdeno",
                            "PreClinCRC_Early",
                            "PreClinCRC_Late",
                            "ClinCRC_Early",
                            "ClinCRC_Late",
                            "CRC_Death",
                            "OC_Death") 
    n.s.hist <- length(v.name.states.hist)
    
    ### Transition Intensity Array
    a.G.hist <- array(0, dim = list(n.s.hist, n.s.hist, n.t), 
                      dimnames = list(v.name.states.hist, v.name.states.hist, v.ages))
    
    ## Age-dependent rate
    lambda_1_t <- lambda_1*g*v.ages^(g-1)
    
    ## Fill Transition Intensity Array
    a.G.hist["Normal",     "Normal", ]                 <- -(lambda_1_t + lambda_asr)
    a.G.hist["Normal",     "SmallAdeno", ]             <- lambda_1_t
    a.G.hist["Normal",     "OC_Death", ]               <- lambda_asr
    
    a.G.hist["Normal_Hist", "Normal_Hist", ]           <- -(lambda_1_t*hr.Hist + lambda_asr)
    a.G.hist["Normal_Hist", "SmallAdeno", ]            <- lambda_1_t*hr.Hist
    a.G.hist["Normal_Hist", "OC_Death", ]              <- lambda_asr
    
    a.G.hist["SmallAdeno", "SmallAdeno", ]             <- -(lambda_2 + lambda_asr)
    a.G.hist["SmallAdeno", "LargeAdeno", ]             <- lambda_2  
    a.G.hist["SmallAdeno",     "OC_Death", ]           <- lambda_asr    
    a.G.hist["LargeAdeno", "LargeAdeno", ]             <- -(lambda_3 + lambda_asr)
    a.G.hist["LargeAdeno", "PreClinCRC_Early", ]       <- lambda_3 
    a.G.hist["LargeAdeno",     "OC_Death", ]           <- lambda_asr    
    a.G.hist["PreClinCRC_Early", "PreClinCRC_Early", ] <- -(lambda_4 + lambda_5 + lambda_asr)
    a.G.hist["PreClinCRC_Early", "PreClinCRC_Late", ]  <- lambda_4
    a.G.hist["PreClinCRC_Early", "ClinCRC_Early", ]    <- lambda_5
    a.G.hist["PreClinCRC_Early",     "OC_Death", ]     <- lambda_asr    
    a.G.hist["PreClinCRC_Late", "PreClinCRC_Late", ]   <- -(lambda_6 + lambda_asr)
    a.G.hist["PreClinCRC_Late", "ClinCRC_Late", ]      <- lambda_6
    a.G.hist["PreClinCRC_Late",     "OC_Death", ]      <- lambda_asr
    a.G.hist["ClinCRC_Early", "ClinCRC_Early", ]       <- -(lambda_7 + lambda_asr)
    a.G.hist["ClinCRC_Early", "CRC_Death", ]           <- lambda_7
    a.G.hist["ClinCRC_Early",     "OC_Death", ]        <- lambda_asr
    a.G.hist["ClinCRC_Late", "ClinCRC_Late", ]         <- -(lambda_8 + lambda_asr)
    a.G.hist["ClinCRC_Late", "CRC_Death", ]            <- lambda_8
    a.G.hist["ClinCRC_Late",     "OC_Death", ]         <- lambda_asr
    
    ### Compute Transition proability 
    a.P.hist <- sapply(seq_along(v.ages), 
                       function(x) MatrixExp(a.G.hist[,,x]), 
                       simplify="array")
    
    ### Check if Transition Probability array is valid
    valid <- apply(a.P.hist, 3, function(x) all.equal(sum(rowSums(x)), n.s.hist))
    if (!isTRUE(all.equal(as.numeric(sum(valid)), as.numeric(n.t)))) {
      stop("This is not a valid transition Matrix")
    }
    
    return(a.P.hist) # Return transition probability array
  }
  )
}

cea_crc_screening <- function(v.params.scr){
  out.scr <- crc_screening(v.params.scr = v.params.scr, 
                           SCR = T)
  out.noscr <- crc_screening(v.params.scr = v.params.scr, 
                             SCR = F)
  v.costs <- c(NoScr = out.noscr$ev.costs.disc,
               Scr = out.scr$ev.costs.disc)
  v.qalys <- c(NoScr = out.noscr$ev.qaly.disc,
               Scr = out.scr$ev.qaly.disc)
  return(list(Costs = v.costs,
              QALYs = v.qalys))
}

### Initialize matrices for Costs and QALYs
m.c.informed <- m.e.informed <-matrix(NaN, nrow = n.sim, ncol = 2)
colnames(m.c.informed) <- colnames(m.e.informed) <- v.names.strategies

## Run health economic model
for(i in 1:n.sim){
  out.cea <- cea_crc_screening(v.params.scr = df.psa.params.informed[i, ])
  ### Fill Costs and QALYs matrices
  m.c.informed[i, ] <- out.cea$Costs
  m.e.informed[i, ] <- out.cea$QALYs
}

wtp.vec<-seq(0,120000,by=500)
m<-bcea(m.e.informed,m.c.informed,ref=2,interventions = v.names.strategies,wtp=wtp.vec)
plot(m)

wtp<-75000

#INB for specific willingness-to-pay
INB<-m$ib[which(m$k==wtp),]

#EVPPI
evi.key<-evppi(c("lambda_1","g","prev.adeno"),df.psa.params.informed[1:5000,],m)
plot(evi.key)
evi.key$evppi[which(evi.key$k==wtp)]

#Fitted Values for specific WTP
save<-wtp*evi.key$fitted.effects[,1]-evi.key$fitted.costs[,1]

## Data Generation
#Sample sizes
n.<-c(5,40,100,200,500,750,1000,1500)

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
# Extract age distribution from Canadian life tables.
tables<-hmd.mx2("CAN",username="anna.heath@sickkids.ca",password="1546527368")
probs.age<-tables$pop$total[,"2011"]/sum(tables$pop$total[,"2011"])
ages.vec<-0:110
ages<-rmultinom(15000,1,probs.age)
ages.long<-array(NA,dim=15000)
for(i in 1:15000){
  ages.long[i]<-ages.vec[which(ages[,i]==1)]}
ages.long<-ages.long[which((ages.long>=25)&(ages.long<=90))]

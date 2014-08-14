###############################################################################
###### Double  switch constrained model from Morales et al. 2004 ##############
###### The script fits the model to simulated data ############################
###############################################################################

library(CircStats)
library(rube)


source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Data_preparation.R")
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Functions.R")
(wd<-setwd("I:/Jan_BUGS"))

##################################################
###### Double switch constained model model ######
##################################################

### Define winbugs model
doubleConstrainedModel="model{
  
  ## priors
  
  b[1] ~ dgamma(0.01,0.01)   ## shape parameter for slow movement
  b[2] ~ dgamma(0.01,0.01)I(1.1,)   ## shape parameter for fast movement
  
  a[2] ~ dgamma(0.01, 0.01)  ## scale parameter for fast movement
  eps ~ dnorm(0.0, 0.01)I(0.0,)  ## a nonnegative variate    
  a[1]  <- a[2] + eps		## scale parameter for slow movement
  
  rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
  rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement
  
  mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
  mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement
  
  q[1] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 1 at time t-1
  q[2] ~ dunif(0,1)		## probability of being in state 1 at t given that individual was in state 2 at time t-1
  
  phi[1] ~ dunif(0,1)
  phi[2] <- 1-phi[1]
  idx[1] ~ dcat(phi[])		## asign state for first observation 
  
  Pi <- 3.14159265359		## define pi
  
  
  for (t in 2:NObs) {
    
    nu[t,1] <- q[idx[t-1]]
    nu[t,2] <- 1-q[idx[t-1]]
    idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index
    
    ## likelihood for steps
    Length[t] ~ dgamma(b[idx[t]], a[idx[t]])	# Weibull distriution for step length
    
    ## likelihood for turns.  
    ## use the ones trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution
    
    ones[t] <- 1
    ones[t] ~ dbern(wC[t])
    ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
    wC[t] <- ( 1/(2*Pi)*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t]-mu[idx[t]])) )/500
    
  }
}"

####################################################################
###### Simulate data and test double switch constrained model ######
####################################################################

### Simulation of step length  and turning angle for two behavioural movement modes without switching
# Simulate step length
simstep1<-rgamma(500, shape=0.543, scale = 1.028)
simstep2<-rgamma(500, shape=4, scale = 2)
simstep<-c(simstep1, simstep2)
simstep<-sample(simstep)

# Simulate turning angle
simangle1<-rwrpcauchy(500, 3.171, 0.164)-pi
simangle2<-rwrpcauchy(500, 0.638, 0.4)-pi
simangle<-c(simangle1, simangle2)
simangle<-sample(simangle)

###simulation of switching step length and angle via two phased distributions
# Simulate step length
set.seed(2)
lala<-runif(200)
simstep_tp<-sapply(lala < 0.5, twoPhasedGamma, shape1=0.543, shape2=4, rate1=1.028, rate2=2) #probability set to 0.5

# Simulate angle
set.seed(2)
lala<-runif(200)
simangle_tp<-sapply(lala < 0.5, twoPhasedWrpC, mu1=3.171, mu2=0.638, rho1=0.164, rho2=0.4)-pi #probability set to 0.5

### Bundle data. Change according to which simulation to use
simBUGS <- list(NObs = length(simstep),
                Length = simstep + 0.000001, 
                Angle = simangle)

### Function to generate broad starting values
double.constrained.inits <- function(){
  list( b = runif(n=2, min=1, max=4),
        eps = abs(rnorm(n=1, mean=0, sd=5)),
        #a = runif(n=1, min=1, max=2),
        rho = runif(n=2, min=0.0001, max=0.9999),
        mu = runif(n=2, min=-2, max=2),
        q = runif(n=2, min=0, max=1)
  )   
}

### Function to generate rather fixed starting values
double.constrained.inits <- function(){
  list( b = c(runif(n=1, min=0.5, max=0.6),runif(n=1, min=4.0, max=4.1)),
        #a = c(runif(n=1, min=0.0001, max=5),runif(n=1, min=0.0001, max=5)),
        eps = abs(rnorm(n=1, mean=0, sd=1)),
        rho = c(runif(n=1, min=0.1, max=0.2),runif(n=1, min=0.4, max=0.5)),
        mu = c(runif(n=1, min=2.9, max=3),runif(n=1, min=0.4, max=0.5)) 
        #nu = runif(1, min=0.0001, max=0.9999)
  )
}


### Parameters to be monitored (= to estimate)
double.constrained.params <- c("a", "b", #2-item lists of the two values for the Weibull
                               "rho", "mu", #2-item lists of the two values for the Wrapped Cauchy
                               "q", #probability of staying in state 1 (q[1]) or switching (q[2])
                               "phi") #probability of being in state 1 (phi[1]) or state 2 (phi[2])

### MCMC settings
nc <- 4    			# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
#ns <- 1000

### test with rube before run
summary(rube(model = doubleConstrainedModel, data = simBUGS, inits = double.constrained.inits))


### Run WinBUGS
double.constrained.sim <- rube(
  model = doubleConstrainedModel,
  data = simBUGS,
  inits = double.constrained.inits,
  parameters.to.save = double.constrained.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  #n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,			
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")


### Plot probability density function with curve of estimated parameters
# For step length
hist(simstep, breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.constrained.sim$sims.array[,,3]), 
             rate=mean(double.constrained.sim$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
curve(dgamma(x, shape=mean(double.constrained.sim$sims.array[,,4]), 
             rate=mean(double.constrained.sim$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)

# For angle
hist(simangle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(double.constrained.sim$sims.array[,,5]), 
                 mu=mean(double.constrained.sim$sims.array[,,7])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(double.constrained.sim$sims.array[,,6]), 
                 mu=mean(double.constrained.sim$sims.array[,,8])), type = "l", col="blue", add=T, lwd=2)

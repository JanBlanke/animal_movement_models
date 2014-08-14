###############################################################################################
###### This is the modified double switch model with covariates from Morales et al. 2004 ######
###### The script fits the model to both simulated data and gazelle data ######################
###############################################################################################

library(CircStats)
library(rube) # Wrapper for R2WinBUGS. Download zip from http://www.stat.cmu.edu/~hseltman/rube/ and install from local zip file 

source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Data_preparation_equidist.R") # Run script that prepares gazella data
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Functions.R") # Run script that contains functions for two phased data simulation
(wd<-setwd("C:/WinBUGS_Analyse")) # Use short path to make WinBUGS happy

#################################################
###### Double switch model with covariates ######
#################################################

### Define WinBUGS model

double.covar.model="model{
  
  ### define variables

  Pi <- 3.14159265359    ## define pi

  ### priors
  
  #b[1] ~ dunif(0.25, 20)   ## shape parameter for slow movement
  #b[2] ~ dunif(0.25, 20)   ## shape parameter for fast movement
  b[1] ~ dgamma(0.01,0.01)
  b[2] ~ dgamma(0.01,0.01)

  #a[1] ~ dgamma(0.01, 0.01)    ## Scale parameter for slow movement
  #a[2] ~ dgamma(0.01, 0.01)    ## Scale parameter for fast movement
  #a[1] ~ dunif(0.1, 10) 
  #a[2] ~ dunif(0.1, 10)  
  aTwo ~ dgamma(0.01, 0.01)
  #aTwo ~ dunif(0.01, 10)
  a[2] <- aTwo
  eps ~ dnorm(0.0, 0.01)I(0.0,)  ## A nonnegative variate    
  a[1]  <- aTwo + eps

  rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
  rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

  #mu[1] ~ dunif(0, 6.283185)    ## mean direction of turns for slow movement for simulation 
  #mu[2] ~ dunif(0, 6.283185)	  ## mean direction of turns for fast movement for simulation
  mu[1] ~ dunif(-3.14159265359, 3.14159265359)    ## mean direction of turns for slow movement 
  mu[2] ~ dunif(-3.14159265359, 3.14159265359)

  beta[1] ~ dnorm(0,0.1)    ## intercept of linear predictor
  beta[2] ~ dnorm(0,0.1)    ## intercept of linear predictor

  m[1] ~ dnorm(0,0.1)   ## slope of linear predictor. use for only one covariate to keep m simple
  m[2] ~ dnorm(0,0.1)   ## slope of linear predictor. use for only one covariate to keep m simple

  ## Use for more than one covariate to make m a 2-D array
  #for(i in 1:2){
  #  for(j in 1:2){
  #    m[j,i] ~ dnorm(0,0.1)	# coefficients to relate distance to habitat i to switching rate
  #  }
  #} 

  ## First step
  #for (i in 1:2){
  m[1] <- 0
  #} 
  

  ### first step

  phi[1] ~ dunif(0,1)
  phi[2] <- 1-phi[1]
  idx[1] ~ dcat(phi[])		## asign state for first observation 
  
  ### all other steps
  
  for (t in 2:NObs) {    
    logit.q[t] <- exp(beta[idx[t-1]] + m[idx[t-1]] * Covar1[t] )     
    q[t] <- logit.q[t]/(1 + logit.q[t])    

    nu[t,1]<-q[t]
    nu[t,2]<-1-q[t]    
    idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index
    
    ### likelihood for steps

    Length[t] ~ dgamma(b[idx[t]], a[idx[t]])	## Gamma distribution for step length
    
    ### likelihood for turns using ones trick to sample from the Wrapped Cauchy distribution
        
    ones[t] <- 1
    ones[t] ~ dbern(wC[t])
    ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
    wC[t] <- ((1/(2*Pi))*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t] - mu[idx[t]])))/500
    Angle[t] ~ dunif(-3.14159265359, 3.14159265359) # For imputing NAs
  }

}"


#######################################################################
###### Simulate data and test double switch model with covariates######
#######################################################################

### Functions

twoStateGamma <- function(state, shape1, shape2, rate1, rate2){
  if(state == 1) rgamma(1, shape1, rate1) else if(state == 2) rgamma(1, shape2, rate2)  
}

twoStateWrpC <- function(state, rho1, rho2, mu1, mu2){
  if(state == 1) rwrpcauchy(1, mu1, rho1) else if(state == 2) rwrpcauchy(1, mu2, rho2)
}

### First step of simulation from x1y1 to x2y2

# Define vectors
x<-c(); y<-c(); dx<-c(); dy<-c(); covar<-c(); state<-c(); simstep<-c(); simangle<-c(); abs.angle<-c(); prob<-c()

# Set the distribution parameters

rate1 <- 0.75
rate2 <- 0.25

shape1 <- 1
shape2 <- 3

rho1 <- 0.05
rho2 <- 0.75

mu1 <- 2
mu2 <- 4

# intercept and slope. they need to be vectors with 2 values. the first corresponds to the 1->2 transition and the second to the 2->1 transitions
m <- c(2.3, -1.4)
beta <- c(0.6, -0.3)


prob[1] <- 0
state[1] <- ifelse(runif(1) <= 0.5, 1, 2)
simstep[1]<- twoStateGamma(state[1], shape1, shape2, rate1, rate2)
simangle[1]<-runif(1, 0, 2*pi)
abs.angle[1]<-simangle[1]

x[1]<-0 # x coordinate of starting position
y[1]<-0 # y coordinate of starting position

dx[1]<-simstep[1]*sin(simangle[1]) 
dy[1]<-simstep[1]*cos(simangle[1])

covar[1]<-0

x[2]<-x[1]+dx[1]
y[2]<-y[1]+dy[1]

### Further steps of simulation from x2y2 to xnyn

step.count<-250

for (i in 2:step.count){
  # get the previous state. the parameters of the linear predictor will depend on this, just as the WinBUGS model.
  prev.state <- state[i - 1]
  
  # calculate covariate
  covar[i]<-sin(.5*x[i])*sin(.5*y[i])
  
  # linear predictor. parameters must be indexed by state
  linpred <- beta[prev.state] + m[prev.state] * covar[i]
  
  # logit link
  prob[i]<-exp(linpred)/(1+exp(linpred))
  
  # decision from which distribution to sample
  state[i]<-ifelse(runif(1)<=prob[i],1,2)
  
  # sample from distributions based on decision
  simstep[i]<-twoStateGamma(state[i],shape1, shape2, rate1, rate2)
  simangle[i]<-twoStateWrpC(state[i], rho1, rho2, mu1, mu2)
  abs.angle[i]<-abs.angle[i-1]+simangle[i]
  
  dx[i] <- simstep[i] * sin(abs.angle[i]) 
  dy[i] <- simstep[i] * cos(abs.angle[i])
  
  x[i+1]<-x[i]+dx[i]
  y[i+1]<-y[i]+dy[i]
  
}

### Save and explore simulation result
covar.sim<-as.data.frame(cbind(x,y,c(dx, NA),c(dy,NA),c(abs.angle,NA),c(simangle,NA),c(simstep, NA),c(covar,NA),c(prob,NA),c(state,NA)))
colnames(covar.sim)<-c("x","y","dx","dy","abs.angle","simangle","simstep","covar","prob","state")
covar.sim<-na.omit(covar.sim)

plot(x,y, type="l")

### Bundle simulated data
simBUGS <- list(NObs = length(covar.sim$simstep),
                Length = covar.sim$simstep,
                Angle = covar.sim$simangle,
                Covar = covar.sim$covar)

### Function to generate broad starting values
#double.covar.inits <- function(){
#  list(  b = runif(n=2, min=1, max=2),
#         #eps = runif(n=1, min=0.01, max=0.5),
#         #aTwo = runif(n=1, min=0.15, max=0.35),
#         rho = runif(n=2, min=0.0001, max=0.9999),
#         q = runif(n=2, min=0, max=1),
#         a = runif(n=2, min=0.1, max=1),
#         mu = runif(n=2, min=0, max=6))
#            
#  }

### Function to generate rather fixed starting values
double.covar.inits.fix <- function(){
  list(  b = c(runif(n=1, min=0.75, max=1.25),runif(n=1, min=3.75, max=4.25)),
         a = c(runif(n=1, min=0.25, max=0.75),runif(n=1, min=0.1, max=0.5)),
         #eps = runif(n=1, min=2, max=4),
         #aTwo = runif(n=1, min=0.5, max=1.5),
         rho = c(runif(n=1, min=0.01, max=0.1), runif(n=1, min=0.5, max=0.9)),
         mu = c(runif(n=1, min=1.5, max=2.5), runif(n=1, min=3.5, max=4.5)),
         m = c(rnorm(1,2.3,0.25), rnorm(1,-1.4,0.25)),
         beta = c(rnorm(1, 0.6, 0.25), rnorm(1, -0.3, 0.25)) 
         
  )   
}


### Parameters to be monitored (= to estimate)
double.covar.params <- c("a", "b","rho","mu","m", "beta")

### MCMC settings. Represent settings used by Morales et al. 2004
nc <- 3      		# Number of chains
ni <- 10000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
#ns <- 1000

### Test with rube before run
summary(rube(model = double.covar.model, data = simBUGS, inits = double.covar.inits.fix))

### Run WinBUGS
double.covar.sim <- rube(
  model = double.covar.model,
  data = simBUGS,
  inits = double.covar.inits.fix,
  parameters.to.save = double.covar.params,
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

p3(double.covar.sim)
summary(double.covar.sim)

###Summary & inspection
double.switch.sim$summary
summary(double.switch.sim)
p3(double.switch.sim) #Interactive plot for visualizing WinBUGS output


### Plot probability density functions with curve of estimated simulation parameters
# For step length
hist(simstep, breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.switch.sim$sims.array[,,3]), # ...$sims.array[rows,chains,parameters]
             rate=mean(double.switch.sim$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
curve(dgamma(x, shape=mean(double.switch.sim$sims.array[,,4]), 
             rate=mean(double.switch.sim$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)

# For step length, for each state 
par(mfrow=c(2,1))
hist(rgamma(250, shape=shape1, rate=rate1), breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.switch.sim$sims.array[,,4]), 
             rate=mean(double.switch.sim$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)
hist(rgamma(250, shape=shape2, rate=rate2), breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.switch.sim$sims.array[,,3]), # ...$sims.array[rows,chains,parameters]
             rate=mean(double.switch.sim$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
par(mfrow=c(1,1))

# For turning angle
hist(simangle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(double.switch.sim$sims.array[,,5]), 
                 mu=mean(double.switch.sim$sims.array[,,7])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(double.switch.sim$sims.array[,,6]), 
                 mu=mean(double.switch.sim$sims.array[,,8])), type = "l", col="blue", add=T, lwd=2)

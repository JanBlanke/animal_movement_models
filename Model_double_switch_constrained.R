###############################################################################
###### Double  switch constrained model from Morales et al. 2004 ##############
###### The script fits the model to the gazelle data ##########################
###############################################################################

library(CircStats)
library(rube)


source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Data_preparation_equdist.R") # Run script that prepares gazella data
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Functions.R") # Run script that contains functions for two phased data simulation
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Misc_scripts/ggplot2_function_graphstyle.R") # Source graphical settings for ggplot2

(wd<-setwd("C:/WinBUGS_Analyse")) # Use short path to make WinBUGS happy

##################################################
###### Double switch constained model model ######
##################################################

### Define winbugs model
doubleConstrainedModel="model{
  
  ## priors
  
  b[1] ~ dgamma(0.01,0.01)   ## shape parameter for slow movement
  b[2] ~ dgamma(0.01,0.01) #I(1.1,)   ## shape parameter for fast movement
  
  aTwo ~ dgamma(0.01, 0.01)	## scale parameter for fast movement
  a[2] <- aTwo
  eps ~ dnorm(0.0, 0.01)I(0.0,)  ## a nonnegative variate    
  a[1]  <- aTwo + eps		## scale parameter for slow movement
  
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
    
    Angle[t]~ dunif(-3.14159265359,3.14159265359) # For imputing NAs
  }
}"



### Bundle data
BUGS602191 <- list(NObs = nrow(ltraj.602191A),
                   Length = ltraj.602191A$dist/1000, 
                   Angle = ltraj.602191A$rel.angle)

### Generate starting values
double.constrained.inits <- function(){
  list( b = runif(n=2, min=1, max=5),
        eps = abs(rnorm(n=1, mean=0, sd=10)),
        #a[2] = runif(n=1, min=1, max=5),
        rho = runif(n=2, min=0.0001, max=0.9999),
        mu = runif(n=2, min=-2, max=2),
        q = runif(n=2, min=0, max=1)
        )   
}

double.constrained.inits.602191 <- function(){
  list(  b = c(runif(n=1, min=0, max=6),runif(n=1, min=1, max=2)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.4, max=0.8),
         eps = runif(n=1, min=4, max=6),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.1, max=0.3)),
         mu = c(runif(n=1, min=0.4, max=0.8), runif(n=1, min=-0.1, max=0.1)),
         q = runif(n=2, min=0, max=1)          
  )   
}


### Parameters to be monitored (= to estimate)
double.constrained.params <- c("a", "b","rho","mu", "q", "idx")

### MCMC settings
nc <- 3     		# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
#ns <- 1000

### test with rube before run
summary(rube(model = doubleConstrainedModel,inits = double.constrained.inits.602191, data = BUGS602191 ))

### Run WinBUGS
double.constrained.G602191 <- rube(
  model = doubleConstrainedModel,
  data = BUGS602191,
  inits = double.constrained.inits.602191,
  parameters.to.save = double.constrained.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  #n.sims = ns,
  debug = FALSE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,  		
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")


#### Derive states. Insert specific individuals and change columns accordingly
dim(double.constrained.G602189$sims.matrix)
idx.sample<-double.constrained.G602189$sims.matrix[,11:697] ## Create matrix containing the samples in which each column is an idx parameter
idx.states<-idx.estimate(idx.sample=idx.sample, ncols=687,nrows=6000) ## Posterior summary for idx to gain specific states

### Create dataframe of data and idx for further investigation. Insert specific individuals!
data.602191.idxed<-as.data.frame(cbind(ltraj.602191A$dist/1000, ltraj.602191A$rel.angle,as.Date(ltraj.602191A$date),idx.states))
colnames(data.602191.idxed)<-c("steps", "angles","date", "idx")










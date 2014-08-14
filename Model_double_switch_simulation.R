###############################################################################
###### Double  switch model from Morales et al. 2004 ##########################
###### The script fits the model to both simulated data #######################
###############################################################################

### Load libraries
library(CircStats)
library(rube) # Wrapper for R2WinBUGS. Download zip from http://www.stat.cmu.edu/~hseltman/rube/ and install from local zip file 

### Load required scripts
#source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Data_preparation_no_winter.R") # Run script that prepares gazella data
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Functions.R") # Run script that contains functions for two phased data simulation
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Misc_scripts/ggplot2_function_graphstyle.R") # Source graphical settings for ggplot2

### Set working directory
(wd<-setwd("C:/WinBUGS_Analyse")) # Use short path to make WinBUGS happy


#################################
###### Double switch model ######
#################################

### Define WinBUGS model
double.switch.model=" model{

   ### Priors 

   b[1] ~ dgamma(0.01,0.01)   ## Shape parameter for slow movement
   b[2] ~ dgamma(0.01,0.01)   ## Shape parameter for fast movement
   #b[1] ~ dunif(0.25, 20)    
   #b[2] ~ dunif(0.25, 20)    
   
   #a[1] ~ dgamma(0.01, 0.01)    ## Scale parameter for slow movement
   #a[2] ~ dgamma(0.01, 0.01)    ## Scale parameter for fast movement
   #a[1] ~ dunif(0.1, 10) 
   #a[2] ~ dunif(0.1, 10)  
   aTwo ~ dgamma(0.01, 0.01)
   #aTwo ~ dunif(0.01, 10)
   a[2] <- aTwo
   eps ~ dnorm(0.0, 0.01)I(0.0,)  ## A nonnegative variate    
   a[1]  <- aTwo + eps  	

   rho[1] ~ dunif(0,1) 		## Mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1) 		## Mean cosine of turns for fast movement

   #mu[1] ~ dunif(-3.14159265359, 3.14159265359)  ## Mean direction of turns for slow movement 
   #mu[2] ~ dunif(-3.14159265359, 3.14159265359)  ## Mean direction of turns for fast movement
   mu[1] ~ dunif(0, 6.283185)  
   mu[2] ~ dunif(0, 6.283185)	  
   
   q[1] ~ dunif(0,1) 		## Probability of being in state 1 at t given that individual was in state 1 at time t-1    1->1
   q[2] ~ dunif(0,1) 		## Probability of being in state 1 at t given that individual was in state 2 at time t-1    2->1    

   phi[1] ~ dunif(0,1) 
   phi[2] <- 1-phi[1]
   idx[1] ~ dcat(phi[])   ## Sample state for first observation 
  
   Pi <- 3.14159265359		## Define pi


   for (t in 2:NObs) {

      nu[t,1] <- q[idx[t-1]]
      nu[t,2] <- 1-q[idx[t-1]]
      idx[t] ~ dcat(nu[t,])   ##  Idx is the latent variable and the parameter index

      ### Likelihood for steps
      Length[t] ~ dgamma(b[idx[t]], a[idx[t]])    ## Gamma distribution for step length

      ### Likelihood for turns. Use the 'ones' trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution
      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## Below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ((1/(2*Pi))*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t] - mu[idx[t]])))/500
      
      #Angle[t]~ dunif(-3.14159265359,3.14159265359) # For imputing NAs
  }
}"


########################################################
###### Simulate data and test double switch model ######
########################################################

### Simulate data to feed the double switch model

## Set the switching probabilities as defined above (and in Morales)
q21 <- 0.5
q12 <- 0.25

## Set the distribution parameters
shape1 <- 1
shape2 <- 1.5

rate1 <- 0.5
rate2 <- 0.25

rho1 <- 0.05
rho2 <- 0.75

mu1 <- 1.7
mu2 <- 2.7

hist(rgamma(250,shape = shape2, rate = rate2),breaks=50, col="black")
hist(rgamma(250,shape = shape1, rate = rate1),breaks=30, col="blue", add=T)

## Set the number of samples you want to draw.
samp.size <- 250

## Initialize variables
states <- NULL
simstep <- NULL
simangle <- NULL

## Randomly draw the first state and corresponding step length and turn angle before entering the loop 
states[1] <- sample(1:2, size=1)
simstep[1] <- twoStateGamma(states[1], shape1, shape2, rate1, rate2)
simangle[1] <- twoStateWrpC(states[1], rho1, rho2, mu1, mu2)

## Loop to update the state when state transitions are explicitly modeled (as defined in the function "state.samp") and draw step lengths and turn angles accordingly

for(i in 2:samp.size){
  states[i] <- state.samp(states[i - 1], q21, q12)
  simstep[i] <- twoStateGamma(states[i], shape1, shape2, rate1, rate2)
  simangle[i] <- twoStateWrpC(states[i], rho1, rho2, mu1, mu2)}


### Bundle data. 
{simBUGS <- list(NObs = length(simstep),
                 Length = simstep, 
                 Angle = simangle)}

### Function to generate rather fixed starting values
double.switch.inits.fix <- function(){
  list(  b = c(runif(n=1, min=0.5, max=1.5),runif(n=1, min=3.5, max=4.5)),
         #a = c(runif(n=1, min=0.25, max=0.75),runif(n=1, min=0.1, max=0.5)),
         eps = runif(n=1, min=0.1, max=0.5),
         #eps2 = runif(n=1, min=1, max=3),
         aTwo = runif(n=1, min=0.25, max=0.75),
         #bTwo = runif(n=1, min=0.75, max=1.25),
         rho = c(runif(n=1, min=0.01, max=0.1), runif(n=1, min=0.6, max=0.9)),
         mu = c(runif(n=1, min=1.5, max=2.5), runif(n=1, min=3.5, max=4.5)),
         q = c(runif(n=1, min=0.25, max=0.75), runif(n=1, min=0.1, max=0.5))
  )   
}


### Parameters to be monitored (= to estimate)
double.switch.params <- c("a", "b","rho","mu","q", "idx")

### MCMC settings. Represent settings used by Morales et al. 2004
nc <- 3      		# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate


### Test with rube before run
summary(rube(model = double.switch.model, data = simBUGS, inits = double.switch.inits.fix))

### Run WinBUGS
double.switch.simulation <- rube(
  model = double.switch.model,
  data = simBUGS,
  inits = double.switch.inits.fix,
  parameters.to.save = double.switch.params,
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

###Summary & inspection
summary(double.switch.simulation)
p3(double.switch.sim) #Interactive plot for visualizing WinBUGS output

### Derive states
idx.sample<-double.switch.simulation$sims.matrix[,11:260]
idx.states<-idx.estimate(idx.sample=idx.sample, nrows=4500, ncols=250)

### Create data frame of simulation output for further investigation
sim.output<-as.data.frame(cbind(simstep, simangle,idx.states))

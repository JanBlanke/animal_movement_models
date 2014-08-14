###############################################################################################
###### This is the modified double switch model with covariates from Morales et al. 2004 ###### 
###### The script fits the model to both simulated data and gazelle data ######################
###############################################################################################

library(CircStats)
library(rube) # Wrapper for R2WinBUGS. Download zip from http://www.stat.cmu.edu/~hseltman/rube/ and install from local zip file 

#source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Data_preparation.R") # Run script that prepares gazella data
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Functions.R") # Run script that contains functions for two phased data simulation
(wd<-setwd("I:/Jan_BUGS")) # Use short path to make WinBUGS happy

#################################################
###### Double switch model with covariates ######
#################################################

### Define WinBUGS model
{double.covar.model="model{

  ## priors

  b[1] ~ dunif(0.25, 20)    ## shape parameter for slow movement
  b[2] ~ dunif(0.25, 20)

  a[1] ~ dunif(0.01, 10) 
  a[2] ~ dunif(0.01, 10)

  rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
  rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement

  mu[1] ~ dunif(0, 6.283185)  ## mean direction of turns for slow movement 
  mu[2] ~ dunif(0, 6.283185)	## mean direction of turns for fast movement

  ## priors for habitat types
  for (i in 1:10) {
  mu.phi[i] ~ dnorm(0.0, 0.01)
  }
  
  Pi <- 3.14159265359		## define pi


for (t in 1:NObs) {


    ## movement state is related to current habitat type
    mu.type[t] <- mu.phi[typ[t]] ## typ is a variable that indicates habitat type at current location
    logit.nu[t] ~ dnorm(mu.type[t], 0.01)
    nu_h[t] <- exp(logit.nu[t])/(1 + exp(logit.nu[t]))  ## probability of being in movement type 1
    nu[t,1] <- nu_h[t]
    nu[t,2] <- 1 - nu_h[t]
    idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index

   ## likelihood for steps
   Length[t] ~ dweib(b[idx[t]], a[idx[t]])	# Weibull distriution for step length

   ## likelihood for turns.  
   ## use the ones trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

   ones[t] <- 1
   ones[t] ~ dbern(wC[t])
   ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
   wC[t] <- ((1/(2*Pi))*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t] - mu[idx[t]])))/500
 
  }
}"
}

#######################################################################
###### Simulate data and test double switch model with covariates######
#######################################################################
#Set the switching probabilities as defined above (and in Morales)
{q21 <- 0.5
 q12 <- 0.25
 
 #Set the distribution parameters. Notice that "mu1" and "mu2" are set to Pi in the function calls below.
 shape1 <- 1
 shape2 <- 4
 
 rate1 <- 0.5
 rate2 <- 0.25
 
 rho1 <- 0.05
 rho2 <- 0.75
 
 mu1 <- 1.4
 mu2 <- 4.6}

#Set the number of samples you want to draw.
samp.size <- 250

#Initialize variables
states <- NULL
simstep <- NULL
simangle <- NULL

#Randomly draw the first state and corresponding step length and turn angle before entering the loop 
states[1] <- sample(1:2, size=1)
simstep[1] <- twoStateGamma(states[1], shape1, shape2, rate1, rate2)
simangle[1] <- twoStateWrpC(states[1], rho1, rho2, mu1, mu2)

#Loop to update the state when state transitions are explicitly modeled (as defined in the function "state.samp"),
#and draw step lengths and turn angles accordingly
for(i in 2:samp.size){
  states[i] <- state.samp(states[i - 1], q21, q12)
  simstep[i] <- twoStateGamma(states[i], shape1, shape2, rate1, rate2)
  simangle[i] <- twoStateWrpC(states[i], rho1, rho2, pi, pi)
  
}

gaga<-sample(1:10, 250, replace=TRUE)

### Bundle simulated data
simBUGS <- list(NObs = length(simstep),
                Length = simstep,
                Angle = simangle,
                typ = gaga)

### Function to generate broad starting values
double.covar.inits <- function(){
  list(  b = runif(n=2, min=1, max=2),
         #eps = runif(n=1, min=0.01, max=0.5),
         #aTwo = runif(n=1, min=0.15, max=0.35),
         rho = runif(n=2, min=0.0001, max=0.9999),
         q = runif(n=2, min=0, max=1),
         a = runif(n=2, min=0.1, max=1),
         mu = runif(n=2, min=0, max=6))
            
  }

### Function to generate rather fixed starting values
double.covar.inits.fix <- function(){
  list(  b = c(runif(n=1, min=0.75, max=1.25),runif(n=1, min=3.75, max=4.25)),
         a = c(runif(n=1, min=0.5, max=1),runif(n=1, min=0.1, max=0.5)),
         #eps = runif(n=1, min=2, max=4),
         #aTwo = runif(n=1, min=0.5, max=1.5),
         rho = c(runif(n=1, min=0.0, max=0.5), runif(n=1, min=0.5, max=0.9)),
         mu = c(runif(n=1, min=0, max=2), runif(n=1, min=2, max=4))
         #q = c(runif(n=1, min=0.1, max=0.4), runif(n=1, min=0.4, max=0.99))
         )   
}


### Parameters to be monitored (= to estimate)
double.covar.params <- c("a", "b","rho","mu")

### MCMC settings. Represent settings used by Morales et al. 2004
nc <- 3      		# Number of chains
ni <- 5000			# Number of draws from posterior (for each chain)
nb <- 2500			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
#ns <- 1000

### Test with rube before run
summary(rube(model = double.covar.model, data = simBUGS, inits = double.covar.inits.fix))

### Run WinBUGS
{double.covar.sim <- rube(
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
  bugs.directory ="C:/Program Files/WinBUGS14")}

###Summary & inspection
summary(double.switch.sim) ###!!!
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



#####################################################
###### Fit double switch model to gazelle data ######
#####################################################
{
### Bundle data
G602189BUGS <- list(NObs = nrow(ltraj602189),
                    Length = ltraj602189$dist/1000 + 0.000001, 
                    Angle = ltraj602189$rel.angle)

### Generate starting values
double.switch.inits <- function(){
  list(  b = runif(n=2, min=1, max=4),
         eps = runif(n=1, min=0.01, max=0.5),
         aTwo = runif(n=1, min=0.15, max=0.35),
         rho = runif(n=2, min=0.0001, max=0.9999),
         q = runif(n=2, min=0, max=1)
         )   
}

### Parameters to be monitored (= to estimate)
double.switch.params <- c("a", "b","rho","q")


### MCMC settings
nc <- 4    			# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
ns <- 1000      # Number of simulations

### Run WinBUGS
summary(rube(model = double.switch.model,data = G602189BUGS,inits = double.switch.inits))

double.switch.G602189<- rube(
  model = double.switch.model,
  data = G602189BUGS,
  inits = double.switch.inits,
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

### Plot the posterior distribution of the estimated parameters
samples<-double.switch.G602189$sims.matrix # Samples is a matrix with all the MCMC samples; each variable is a column
hist(samples[,1],main="Posterior of shape", breaks=100) # Change the column of "samples" for plots of further variables

### Plot probability density function with curve of estimated parameters
# For step length
hist(ltraj602189$dist/1000, breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.switch.G602189$sims.array[,,3]), 
             scale=mean(double.switch.G602189$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
curve(dgamma(x, shape=mean(double.switch.G602189$sims.array[,,4]), 
             scale=mean(double.switch.G602189$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)

hist(ltraj602189$rel.angle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(double.switch.G602189$sims.array[,,5]), 
                 mu=pi), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(double.switch.G602189$sims.array[,,6]), 
                 mu=pi), type = "l", col="blue", add=T, lwd=2)


write.csv(double.switch.G602189$summary, "double_switch_G602189.csv")
}

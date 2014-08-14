##############################################################################
###### This is the modified double model from Morales et al. 2004.############
###### The script fits the model to both simulated data and gazelle data######
##############################################################################



library(CircStats)
library(rube)# Wrapper for R2WinBUGS. Download zip from http://www.stat.cmu.edu/~hseltman/rube/ and install from local zip file 


source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Data_preparation.R")
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Functions.R")
(wd<-setwd("I:/Jan_BUGS")) # Use short path to make WinBUGS happy

##########################
###### Double model ######
##########################

### Define double model from the publication of Morales et al. 2004

double.model="model{
  
  ## priors
  
  b[1] ~ dunif(0.25, 20) #dgamma(0.01,0.01)   ## shape parameter for slow movement
  b[2] ~ dunif(0.25, 20) #dgamma(0.01,0.01)   ## shape parameter for fast movement
  
  a[1] ~ dunif(0.1, 10)
  a[2] ~ dunif(0.1, 10)    #dgamma(0.01, 0.01)  ## scale parameter for fast movement
  #eps ~ dnorm(0.0, 0.01)I(0.0,)  ## a nonnegative variate
  #a[1]  <- a[2] + eps		## scale parameter for slow movement
  
  rho[1] ~ dunif(0,1)		## mean cosine of turns for slow movement
  rho[2] ~ dunif(0,1)		## mean cosine of turns for slow movement
  
  #mu[1] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for slow movement 
  #mu[2] ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns for fast movement
  
  Pi <- 3.14159265359		## define pi
  
  
  for (t in 1:NObs) {
    
    nu[t,1] ~ dunif(0,1)    ## probability of being in movement state 1 at time t
    nu[t,2] <- 1 - nu[t,1]
    idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index
    
    ## likelihood for steps
    Length[t] ~ dgamma(b[idx[t]], a[idx[t]])	# Gamma distribution for step length
    
    ## likelihood for turns.  
    ## use the ones trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution
    
    ones[t] <- 1
    ones[t] ~ dbern(wC[t])
    ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
    wC[t] <- ((1/(2*Pi))*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t] - Pi)))/500
  }
}"



#################################################
###### Simulate data and test double model ######
#################################################

### Simulation of step length  and turning angle for two behavioural movement modes without switching
# Simulate step length
#simstep1<-rgamma(100, shape=1.1, scale = 1.028)
#simstep2<-rgamma(100, shape=5, scale = 1.5)
#simstep<-c(simstep1, simstep2)
#simstep<-sample(simstep)

# Simulate turning angle
#simangle1<-rwrpcauchy(100, pi, 0.164)-pi
#simangle2<-rwrpcauchy(100, pi, 0.5)-pi
#simangle<-c(simangle1, simangle2)
#simangle<-sample(simangle)

### Simulation of switching step length and angle via two phased distributions
# Simulate step length
lala<-runif(250)
simstep<-sapply(lala < 0.5, twoPhasedGamma, shape1=1, shape2=4, rate1=0.33, rate2=0.25) #probability set to 0.5

# Simulate angle
simangle<-sapply(lala < 0.5, twoPhasedWrpC, mu1=pi, mu2=pi, rho1=0.05, rho2=0.5) #probability set to 0.5

### Bundle data. Substitute "simstep" to apply simulation with switching
simBUGS <- list(NObs = length(simstep),
                Length = simstep, 
                Angle = simangle
                )

### Function to generate broad starting values
double.inits <- function(){
  list( b = c(runif(n=1, min=0.25, max=3),runif(n=1, min=2, max=6)),
        a = c(runif(n=1, min=0.1, max=1), runif(n=1, min=0.1, max=1)),
        #eps = abs(rnorm(n=1, mean=0, sd=1)),
        rho = c(runif(n=1, min=0.1, max=0.9),runif(n=1, min=0.1, max=0.9))
        #mu = c(runif(n=1, min=2.9, max=3),runif(n=1, min=0.4, max=0.5)) 
        #nu = runif(1, min=0.4, max=0.5)
        )
}


### Function to generate rather fixed starting values
double.inits.fix <- function(){
  list( b = c(runif(n=1, min=0.75, max=1.5),runif(n=1, min=3.5, max=4.5)),
        a = c(runif(n=1, min=0.25, max=0.5), runif(n=1, min=0.15, max=0.35)),
        #eps = abs(rnorm(n=1, mean=0, sd=1)),
        rho = c(runif(n=1, min=0.0, max=0.1),runif(n=1, min=0.45, max=0.55))
        #mu = c(runif(n=1, min=2.9, max=3),runif(n=1, min=0.4, max=0.5)) 
        #nu = runif(1, min=0.4, max=0.5)
        )
}


### Parameters to be monitored (= to estimate)
double.params <- c("a","b", "rho","idx")

### MCMC settings. Represent settings used by Morales et al. 2004
nc <- 3      		# Number of chains
ni <- 5000			# Number of draws from posterior (for each chain)
nb <- 2500			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
ns <- 1000     # Number of simulations


### Test with rube before run WinBUGS
summary(rube(model = double.model, data = simBUGS, inits = double.inits))

### Run WinBUGS
double.sim<- rube(
  model = double.model, 
  data = simBUGS,
  inits=double.inits,
  parameters.to.save = double.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  #n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,
  over.relax=TRUE,
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")  

###Summary & inspection
summary(double.sim)
p3(double.sim)

### Plot probability density functions with estimated parameters
# For step length
hist(simstep, breaks=100, freq=F) # Substitute simstep1/simstep2 to show one component
curve(dgamma(x, shape=mean(double.sim$sims.array[,,3]), # ...$sims.array[rows,chains,parameters]
               rate=mean(double.sim$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
curve(dgamma(x, shape=mean(double.sim$sims.array[,,4]), 
               rate=mean(double.sim$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)

# For turning angle
hist(simangle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(double.sim$sims.array[,,5]), 
                 mu=mean(double.sim$sims.array[,,7])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(double.sim$sims.array[,,6]), 
                 mu=mean(double.sim$sims.array[,,8])), type = "l", col="blue", add=T, lwd=2)

### save summary to file
#write.csv(double.sim$summary,file = "double_sim.csv")



##############################################
###### Fit double model to gazelle data ######
##############################################

### Bundle data
G602189BUGS <- list(NObs = nrow(ltraj.602189),
                    Length = ltraj.602189$dist/1000 + 0.000001, 
                    Angle = ltraj.602189$rel.angle)

### Function to generate broad starting values
double.inits <- function(){
  list( b = runif(n=2, min=0.5, max=10),
        #a = runif(n=1, min=1, max=3),
        #eps = abs(rnorm(n=1, mean=0, sd=10)),
        rho = runif(n=2, min=0.0001, max=0.9999)
        #mu = runif(n=2, min=-3.1, max=3.1) 
        #nu = runif(1, min=0.0001, max=0.9999)
        )
}

### Parameters to be monitored (= to estimate)
double.params <- c("a","b", # 2-item lists of the two parameters for the gamma distribution
                   "rho") #2-item lists of the two parameters for the Wrapped Cauchy


### MCMC settings
nc <- 3  				# Number of chains
ni <- 50000			# Number of draws from posterior (for each chain)
nb <- 25000			# Number of draws to discard as burn-in
nt <- 10					# Thinning rate
ns <- 1000      # Number of simulations

### Run WinBUGS

summary(rube(model = double.model,data = G602189BUGS,inits = double.inits))

double.G602189<- rube(
  model = double.model,
  data = G602189BUGS,
  inits = double.inits,
  parameters.to.save = double.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  #n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,  
  over.relax=TRUE,
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")

###Summary & inspection
summary(double.G602189)
p3(double.sim
   
### Plot probability density function of estimated parameters
hist(ltraj602189$dist/1000, breaks=100, freq=F)
curve(dgamma(x, shape=mean(double.G602189$sims.array[,,3]), 
             scale=mean(double.G602189$sims.array[,,1])), type = "l", col="blue", add=T,lwd=2)
curve(dgamma(x, shape=mean(double.G602189$sims.array[,,4]), 
             scale=mean(double.G602189$sims.array[,,2])), type = "l", col="red", add=T,lwd=2)

hist(ltraj602189$rel.angle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(double.G602189$sims.array[,,5]), 
                 mu=mean(double.G602189$sims.array[,,7])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(double.G602189$sims.array[,,6]), 
                 mu=mean(double.G602189$sims.array[,,8])), type = "l", col="blue", add=T, lwd=2)






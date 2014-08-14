library(CircStats)
library(rube)

(wd<-setwd("I:/Jan_BUGS"))
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Data_preparation.R")
source("//Nzc-ap03/crcgis/Projects/Jan_Blanke/6_Masters_thesis/Rscripts/Functions.R")

#################################
###### Triple switch model ######
#################################

### Define winbugs model
tripleSwitchModel="model{

   ## priors
   ## shape parameters for step length
   shape[1] ~ dgamma(0.01,0.01)
   shape[2] ~ dgamma(0.01,0.01)
   shape[3] ~ dgamma(0.01,0.01) 

   eps1 ~  dnorm(0.0, 0.01)I(0.0,)
   eps2 ~ dnorm(0.0, 0.01)I(0.0,)
   ## rate parameters for step length
   lowerrate ~ dgamma(0.01, 0.01)
   rate[3] <- lowerrate 
   rate[2] <- rate[3] + eps1 
   rate[1] <- rate[2] + eps2

   ## mean cosine for turns
   rho[1] ~ dunif(0,1)
   rho[2] ~ dunif(0,1)
   rho[3] ~ dunif(0,1)

   ## mean direction for turns
   mu[1] ~ dunif(-3.14159265359, 3.14159265359) 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)
   mu[3] ~ dunif(-3.14159265359, 3.14159265359)


   ##  priors for the probability of switching from anything to 1
   qq[1] ~ dunif(0,1)
   qq[2] ~ dunif(0,1)
   qq[3] ~ dunif(0,1)
   q[1] ~ dunif(0,1)
   q[2] ~ dunif(0,1)
   q[3] ~ dunif(0,1)
   ## asign state for first observation 
   
  phi[1]~dunif(0,1)
  zwerg<-1-phi[1]
  phi[2]~dunif(0,zwerg)
  phi[3]<-zwerg - phi[2]

  idx[1] ~ dcat(phi[])    
  
   Pi <- 3.14159265359		## define pi


   for (t in 2:NObs) {

      nu[t,1] <- q[idx[t-1]]
      nu[t,2] <- (1 -q [idx[t-1]] ) * qq[idx[t-1]] 
      nu[t,3] <- (1 -q [idx[t-1]] ) * (1-qq[idx[t-1]] )

      idx[t] ~ dcat(nu[t,])   ##  idx is the latent variable and the parameter index

      ## likelihood for steps
      Length[t] ~ dgamma(shape[idx[t]], rate[idx[t]])	# Weibull distriution for step length

      ## likelihood for turns.  
      ## use the ones trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ( 1/(2*Pi)*(1-rho[idx[t]]*rho[idx[t]])/(1+rho[idx[t]]*rho[idx[t]]-2*rho[idx[t]]*cos(Angle[t]-mu[idx[t]])) )/500

  }
}"


##########################################
###### Simulate data and test model ######
##########################################

### Simulation of step length  and turning angle for three behavioural movement modes without switching
# Simulate step length
simstep1<-rgamma(100, shape=0.543, rate = 1.028)
simstep2<-rgamma(100, shape=15, rate = 6)
simstep3<-rgamma(100, shape=0.7, rate=4)
simstep<-c(simstep1, simstep2, simstep3)
simstep<-sample(simstep)

hist(simstep1, breaks=100, col=2)
hist(simstep2, breaks=100, col=3, add=T)
hist(simstep3, breaks=100, col=4, add=T)
hist(simstep, breaks=100)

# Simulate turning angle
simangle1<-rwrpcauchy(100, 3.171, 0.164)-pi
simangle2<-rwrpcauchy(100, 0.638, 0.4)-pi
simangle3<-rwrpcauchy(100, 1.5, 0.8)-pi
simangle<-c(simangle1, simangle2, simangle3)
simangle<-sample(simangle)

hist(simangle1, breaks=100, col=2)
hist(simangle2, breaks=100, col=3, add=T)
hist(simangle3, breaks=100, col=4, add=T)
hist(simangle, breaks=100)


### Simulation of switching step length and angle
# Simulate step length
set.seed(2)
lala<-runif(700)
simstep<-sapply(lala < 0.5, twoPhasedGamma, shape1=0.543, shape2=15, rate1=1.028, rate2=6) #probability set to 0.5

hist(simstep, breaks=100)

# Simulate turning angle
set.seed(2)
lala<-runif(700)
simangle<-sapply(lala < 0.5, twoPhasedWrpC, mu1=3.171, mu2=0.638, rho1=0.164, rho2=0.4)-pi #probability set to 0.5

hist(simangle, breaks=100)


### Bundle data
simBUGS <- list(NObs = length(simstep),
                Length = simstep + 0.000001, 
                Angle = simangle)

### Generate starting values
triple.switch.inits <- function(){
  list( shape = runif(n=3, min=1, max=5),
        eps1 = abs(rnorm(n=1, mean=0, sd=10)),
        eps2 = abs(rnorm(n=1, mean=0, sd=10)),
        lowerrate = runif(n=1, min=1, max=5),
        rho = runif(n=3, min=0.0001, max=0.9999),
        mu = runif(n=3, min=-2, max=2),
        q = runif(n=3, min=0, max=1),
        qq = runif(n=3, min=0, max=1)
        #phi = runif(n=1, min=0, max=1)
        #			nu = structure(
        #					.Data=c(rep(0, times=(2*192))),
        #					.Dim=c(192,2))
        )  
  
}


### Parameters to be monitored (= to estimate)
triple.switch.params <- c("shape", "rate", #2-item lists of the two values for the Weibull
                          "rho", "mu", #2-item lists of the two values for the Wrapped Cauchy
                          "q","qq", #probability of staying in state 1 (q[1]) or switching (q[2])
                          "phi") #probability of being in state 1 (phi[1]) or state 2 (phi[2])

### MCMC settings
nc <- 4    			# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
#ns <- 1000

### Run WinBUGS
triple.switch.sim <- rube(
  model = tripleSwitchModel,
  data = simBUGS,
  inits = triple.switch.inits,
  parameters.to.save = triple.switch.params,
  #n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  #n.burnin = nb,
  #n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,			
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")


### Plot the posterior distribution for the estimated parameters
samples<-triple.switch.sim$sims.matrix # Samples is a matrix with all the MCMC samples; each variable is a column
hist(samples[,3],main="Posterior of shape1", breaks=100) # Change the column of samples for further posterior distributions

### Plot probability density functions with curve of estimated parameters
# For step length
hist(simstep, breaks=100, freq=F)
curve(dgamma(x, shape=mean(triple.switch.sim$sims.list$shape[1]), 
             rate=mean(triple.switch.sim$sims.list$rate[1])), type = "l", col="red", add=T, lwd=2)
curve(dgamma(x, shape=mean(triple.switch.sim$sims.list$shape[2]), 
             rate=mean(triple.switch.sim$sims.list$rate[2])), type = "l", col="blue", add=T, lwd=2)
curve(dgamma(x, shape=mean(triple.switch.sim$sims.list$shape[3]), 
             rate=mean(triple.switch.sim$sims.list$rate[3])), type = "l", col="green", add=T, lwd=2)


# For turning angle
hist(simangle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(triple.switch.sim$sims.list$rho[1]), 
                 mu=mean(triple.switch.sim$sims.list$mu[1])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(triple.switch.sim$sims.list$rho[2]), 
                 mu=mean(triple.switch.sim$sims.list$mu[2])), type = "l", col="blue", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(triple.switch.sim$sims.list$rho[3]), 
                 mu=mean(triple.switch.sim$sims.list$mu[3])), type = "l", col="green", add=T, lwd=2)



#####################################################
###### Fit double switch model to gazelle data ######
#####################################################

### Bundle data
G602189BUGS <- list(NObs = nrow(ltraj602189),
                    Length = ltraj602189$dist/1000 + 0.000001, 
                    Angle = ltraj602189$rel.angle)

### Generate starting values
triple.switch.inits <- function(){
  list( shape = runif(n=3, min=1, max=5),
        eps1 = abs(rnorm(n=1, mean=0, sd=10)),
        eps2 = abs(rnorm(n=1, mean=0, sd=10)),
        lowerrate = runif(n=1, min=1, max=5),
        rho = runif(n=3, min=0.0001, max=0.9999),
        mu = runif(n=3, min=-2, max=2),
        q = runif(n=3, min=0, max=1),
        qq = runif(n=3, min=0, max=1)
        #phi = runif(n=1, min=0, max=1)
        #  		nu = structure(
        #					.Data=c(rep(0, times=(2*192))),
        #					.Dim=c(192,2))
        )  
  
}

### Parameters to be monitored (= to estimate)
triple.switch.params <- c("shape", "rate", #2-item lists of the two values for the Weibull
                          "rho", "mu", #2-item lists of the two values for the Wrapped Cauchy
                          "q","qq", #probability of staying in state 1 (q[1]) or switching (q[2])
                          "phi") #probability of being in state 1 (phi[1]) or state 2 (phi[2])


### MCMC settings
nc <- 4    			# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
ns <- 1000      # Number of simulations

triple.switch.602189 <- rube(
  model = tripleSwitchModel,
  data = G602189BUGS,
  inits = triple.switch.inits,
  parameters.to.save = triple.switch.params,
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
  bugs.directory ="C:/Programme/WinBUGS14")

### Plot the posterior distribution for the estimated parameters
samples<-triple.switch.602189$sims.matrix # Samples is a matrix with all the MCMC samples; each variable is a column
hist(samples[,3],main="Posterior of shape1", breaks=100) # Change the column of samples for further posterior distributions

### Plot probability density functions with curve of estimated parameters
# For step length
hist(ltraj602189$dist/1000, breaks=100, freq=F)
curve(dgamma(x, shape=mean(triple.switch.602189$sims.list$shape[1]), 
             rate=mean(triple.switch.602189$sims.list$rate[1])), type = "l", col="red", add=T, lwd=2)
curve(dgamma(x, shape=mean(triple.switch.602189$sims.list$shape[2]), 
             rate=mean(triple.switch.602189$sims.list$rate[2])), type = "l", col="blue", add=T, lwd=2)
curve(dgamma(x, shape=mean(triple.switch.602189$sims.list$shape[3]), 
             rate=mean(triple.switch.602189$sims.list$rate[3])), type = "l", col="green", add=T, lwd=2)


# For turning angle
hist(ltraj602189$rel.angle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(triple.switch.602189$sims.list$rho[1]), 
                 mu=mean(triple.switch.602189$sims.list$mu[1])), type = "l", col="red", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(triple.switch.602189$sims.list$rho[2]), 
                 mu=mean(triple.switch.602189$sims.list$mu[2])), type = "l", col="blue", add=T, lwd=2)
curve(dwrpcauchy(x, rho=mean(triple.switch.602189$sims.list$rho[3]), 
                 mu=mean(triple.switch.602189$sims.list$mu[3])), type = "l", col="green", add=T, lwd=2)






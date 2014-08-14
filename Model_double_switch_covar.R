###############################################################################################
###### This is the modified double switch model with covariates from Morales et al. 2004 ######
###### The script fits the model to the gazelle data ##########################################
###############################################################################################

library(rube)

#source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Data_preparation_equdist.R") # Run script that prepares gazella data
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Functions.R") # Run script that contains functions for two phased data simulation
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

  rho[1] ~ dunif(0,1)  	## mean cosine of turns for slow movement
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
  #m[1,i] <- 0
  #} 
    
  m[1]<-0

  ### first step

  phi[1] ~ dunif(0,1)
  phi[2] <- 1-phi[1]
  idx[1] ~ dcat(phi[])		## asign state for first observation 
  
  ### all other steps
  
  for (t in 2:NObs) {    
    logit.q[t] <- exp(beta[idx[t-1]] + m[idx[t-1]] * cov1[t])     
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


### Treatment
#ltraj.602189A<-na.omit(ltraj.602189A)
#ltraj.602190A<-na.omit(ltraj.602190A)
#ltraj.602191A<-na.omit(ltraj.602191A)

### Center and standardize
# For 602189A
#sdize.602189A.temp<-sdize(ltraj.602189A.growth$temp)
sdize.602189A.NDVI<-scale(ltraj.602189A.growth$NDVI, center=T, scale=T)
#sdize.602189A.NDVI[c(191,193,557)]<-0
sdize.602189A.NDVI.buffer<-scale(ltraj.602189A.winter$NDVI.jan.buffer, center=T, scale=T)
#sdize.602189A.NDVI.buffer[557]<-0

# For 602190A
# For 602189A
#sdize.602189A.temp<-sdize(ltraj.602189A.growth$temp)
sdize.602190A.NDVI<-scale(ltraj.602190A.winter$NDVI.jan, center=T, scale=T)
#sdize.602190A.NDVI[546]<-0
sdize.602190A.NDVI.buffer<-scale(ltraj.602190A.winter$NDVI.jan.buffer, center=T, scale=T)
#sdize.602190A.NDVI.buffer[546]<-0

# For 602191A
#sdize.602191A.temp<-sdize(ltraj.602191A$temp)
#sdize.602191A.NDVI<-sdize(ltraj.602191A$NDVI)

### Bundle data
BUGS602189 <- list(NObs = nrow(ltraj.602189A.growth),
                    Length = ltraj.602189A.growth$dist/1000, 
                    Angle = ltraj.602189A.growth$rel.angle,
                    cov1 = as.vector(sdize.602189A.NDVI)#,
                    #cov2 = as.vector(sdize.602189A.NDVI.buffer)
                   )

BUGS602190 <- list(NObs = nrow(ltraj.602190A.winter),
                    Length = ltraj.602190A.winter$dist/1000, 
                    Angle = ltraj.602190A.winter$rel.angle,
                   cov1 = as.vector(sdize.602190A.NDVI),
                   cov2 = as.vector(sdize.602190A.NDVI.buffer))

BUGS602191 <- list(NObs = nrow(ltraj602191A),
                   Length = ltraj602191A$dist/1000, 
                   Angle = ltraj602191A$rel.angle,
                   #Covar1 = sdize.602191A.temp,
                   NDVI = sdize.602191A.NDVI)

### Generate starting values
# For gazelle 602189
double.covar.inits <- function(){
  list(  b = c(runif(n=1, min=1.1, max=1.3),runif(n=1, min=1.5, max=1.8)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.2, max=0.4)),
         aTwo = runif(n=1, min=0.2, max=0.5),
         eps = runif(n=1, min=0.1, max=0.6),
         rho = c(runif(n=1, min=0, max=0.1), runif(n=1, min=0.3, max=0.7)),
         mu = c(runif(n=1, min=0.1, max=0.5), runif(n=1, min=0.01, max=0.05)),
         #m = c(rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1)),
         beta = c(rnorm(1,0,1), rnorm(1,0,1))          
  )   
}

# For gazelle 602190
double.covar.inits <- function(){
  list(  b = c(runif(n=1, min=1.2, max=1.4),runif(n=1, min=1.2, max=1.4)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.2, max=0.4)),
         aTwo = runif(n=1, min=0.3, max=0.6),
         eps = runif(n=1, min=0.1, max=0.7),
         rho = c(runif(n=1, min=0.0, max=0.2), runif(n=1, min=0.3, max=0.6)),
         mu = c(runif(n=1, min=0.1, max=0.4), runif(n=1, min=-0.4, max=-0.2)),
         #m = c(rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1)),
         beta = c(rnorm(1,0,1), rnorm(1,0,1))          
  )   
}

# For gazelle 602191
double.covar.inits <- function(){
  list(  b = c(runif(n=1, min=0.5, max=1),runif(n=1, min=0.8, max=1.2)),
         #a = c(runif(n=1, min=5.4, max=5.6),runif(n=1, min=0.5, max=0.7)),
         aTwo = runif(n=1, min=0.2, max=0.6),
         eps = runif(n=1, min=0.1, max=0.8),
         rho = c(runif(n=1, min=0, max=0.2), runif(n=1, min=0.2, max=0.4)),
         mu = c(runif(n=1, min=1.4, max=1.8), runif(n=1, min=0, max=0.1)),
         #m = c(rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1),rnorm(1,0,1)),
         beta = c(rnorm(1,0,1), rnorm(1,0,1))          
  )   
}

### Parameters to be monitored (= to estimate)
double.covar.params <- c("a", "b","rho","mu","beta","m") 


### MCMC settings
nc <- 3    			# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
ns <- 1000      # Number of simulations

### Run WinBUGS
# Run gazelle 602189
# Growth!!!
summary(rube(model = double.covar.model,data = BUGS602189,inits = double.covar.inits))
double.covar.G602189.winter<- rube(
  model = double.covar.model,
  data = BUGS602189,
  inits = double.covar.inits,
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
summary(double.covar.G602189)


# Run gazelle 602190
summary(rube(model = double.covar.model,data = BUGS602190,inits = double.covar.inits))
### GROWTH!!!
double.covar.G602190.winter <- rube(
  model = double.covar.model,
  data = BUGS602190,
  inits = double.covar.inits,
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
summary(double.covar.G602190)

# Run gazelle 602191
summary(rube(model = double.covar.model,data = BUGS602191,inits = double.covar.inits))
double.covar.G602191<- rube(
  model = double.covar.model,
  data = BUGS602191,
  inits = double.covar.inits,
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
summary(double.covar.G602191)



### Save objects and summary
setwd("C:/Users/Janni/Documents/1_Masterarbeit/Ergebnisse/Final_results")

save(double.covar.G602189.growth, file="602189_growth_thomas_NDVI.RData")

write.csv(double.covar.G602190.winter$summary, "G602190_winter_NDVI_3_versions.csv")

setwd("C:/Users/Janni/Desktop")
save.image("extract.RData")



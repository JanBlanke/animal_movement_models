########################################################### 
###### Double  switch model from Morales et al. 2004 ######
###### The script fits the model to gazelle data ##########
###########################################################

### Load libraries
library(CircStats)
library(rube) # Wrapper for R2WinBUGS. Download zip from http://www.stat.cmu.edu/~hseltman/rube/ and install from local zip file 

### Load required scripts
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Data_preparation_equdist.R") # Run script that prepares gazella data
source("C:/Users/Janni/Documents/1_Masterarbeit/R-skripte/Functions.R") # Run script that contains functions for two phased data simulation

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
   
   #a[1] ~ dgamma(0.01, 0.01)    ## Scale parameter for slow movement
   #a[2] ~ dgamma(0.01, 0.01)    ## Scale parameter for fast movement
   
   aTwo ~ dgamma(0.01, 0.01)
   #aTwo ~ dunif(0.01, 10)
   a[2] <- aTwo
   eps ~ dnorm(0.0, 0.01)I(0.0,)  ## A nonnegative variate    
   a[1]  <- aTwo + eps		

   rho[1] ~ dunif(0,1) 		## Mean cosine of turns for slow movement
   rho[2] ~ dunif(0,1) 		## Mean cosine of turns for fast movement

   mu[1] ~ dunif(-3.14159265359, 3.14159265359)  ## Mean direction of turns for slow movement 
   mu[2] ~ dunif(-3.14159265359, 3.14159265359)  ## Mean direction of turns for fast movement
      
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
      
      Angle[t]~ dunif(-3.14159265359,3.14159265359) # For imputing NAs
  }
}"


#####################################################
###### Fit double switch model to gazelle data ######
#####################################################

### Exclude false last step
ltraj.602189A<-ltraj.602189A[1:703,]
ltraj.602190A<-ltraj.602190A[1:687,]
ltraj.602191A<-ltraj.602191A[1:685,]


### Bundle data
BUGS602189 <- list(NObs = nrow(ltraj.602189A),
                    Length = ltraj.602189A$dist/1000, 
                    Angle = ltraj.602189A$rel.angle)

BUGS602190 <- list(NObs = nrow(ltraj.602190A),
                    Length = ltraj.602190A$dist/1000, 
                    Angle = ltraj.602190A$rel.angle)

BUGS602191 <- list(NObs = nrow(ltraj.602191A),
                    Length = ltraj.602191A$dist/1000, 
                    Angle = ltraj.602191A$rel.angle)

BUGS618669B <- list(NObs = nrow(ltraj.618669B),
                    Length = ltraj.618669B$dist/1000, 
                    Angle = ltraj.618669B$rel.angle)

BUGS618663B <- list(NObs = nrow(ltraj.618663B),
                    Length = ltraj.618663B$dist/1000, 
                    Angle = ltraj.618663B$rel.angle)


### Generate starting values 
## For  gazelle 602189
double.switch.inits.602189 <- function(){
  list(  b = c(runif(n=1, min=0.8, max=1),runif(n=1, min=0.9, max=1.1)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.2, max=0.4),
         eps = runif(n=1, min=0.1, max=0.5),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.2, max=0.4)),
         mu = c(runif(n=1, min=0.1, max=2), runif(n=1, min=0.01, max=0.03)),
         q = runif(n=2, min=0, max=1)          
  )   
}

double.switch.inits.602189 <- function(){
  list(  b = c(runif(n=1, min=0, max=1),runif(n=1, min=0, max=1.1)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.1, max=0.4),
         eps = runif(n=1, min=0.1, max=1),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.2, max=0.4)),
         mu = c(runif(n=1, min=0.1, max=2), runif(n=1, min=0.01, max=0.03)),
         q = runif(n=2, min=0, max=1)          
  )   
}



## For  gazelle 602190
double.switch.inits.602190 <- function(){
  list(  b = c(runif(n=1, min=1, max=3),runif(n=1, min=.8, max=1.2)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=.8, max=1),
         eps = runif(n=1, min=.5, max=1.5),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.4, max=0.6)),
         mu = c(runif(n=1, min=0.1, max=0.5), runif(n=1, min=-0.5, max=-0.1)),
         q = runif(n=2, min=0, max=1)          
  )   
}

## For  gazelle 602191
double.switch.inits.602191 <- function(){
  list(  #b = c(runif(n=1, min=4, max=6),runif(n=1, min=0, max=2)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.4, max=0.8),
         eps = runif(n=1, min=4, max=6),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.1, max=0.3)),
         mu = c(runif(n=1, min=0.4, max=0.8), runif(n=1, min=-0.1, max=0.1)),
         q = runif(n=2, min=0, max=1)          
  )   
}

## For  gazelle 618669B
double.switch.inits.618669 <- function(){
  list(  b = c(runif(n=1, min=2, max=5),runif(n=1, min=2, max=3)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.2, max=0.5),
         eps = runif(n=1, min=1, max=2),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.09, max=0.15)),
         mu = c(runif(n=1, min=0.6, max=1), runif(n=1, min=0.8, max=1.2)),
         q = runif(n=2, min=0, max=1)          
  )   
}

## For  gazelle 618663B
double.switch.inits.618663 <- function(){
  list(  b = c(runif(n=1, min=1, max=3),runif(n=1, min=5, max=7)),
         #a = c(runif(n=1, min=0.7, max=0.9),runif(n=1, min=0.1, max=0.5)),
         aTwo = runif(n=1, min=0.2, max=0.5),
         eps = runif(n=1, min=2, max=3),
         rho = c(runif(n=1, min=0.1, max=0.3), runif(n=1, min=0.09, max=0.15)),
         mu = c(runif(n=1, min=0.6, max=1), runif(n=1, min=0.8, max=1.2)),
         q = runif(n=2, min=0, max=1)          
  )   
}


### Parameters to be monitored (= to estimate)
double.switch.params <- c("a", "b","rho","mu", "q", "idx")

### MCMC settings
nc <- 3    			## Number of chains
ni <- 20000			## Number of draws from posterior (for each chain)
nb <- 5000			## Number of draws to discard as burn-in
nt <- 10				## Thinning rate
ns <- 1000      ## Number of simulations

### Check with rube
summary(rube(model = double.switch.model, data = BUGS602189, inits = double.switch.inits.602189))
summary(rube(model = double.switch.model, data = BUGS602190, inits = double.switch.inits.602190))
summary(rube(model = double.switch.model, data = BUGS602191, inits = double.switch.inits.602191))
summary(rube(model = double.switch.model, data = BUGS618669B, inits = double.switch.inits.618669))
summary(rube(model = double.switch.model, data = BUGS618663B, inits = double.switch.inits.618663))

### Run WinBUGS

## Gazelle 602189A
double.switch.G602189<- rube(
  model = double.switch.model,
  data = BUGS602189,
  inits = double.switch.inits.602189,
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
summary(double.switch.G602189)

## Gazelle 602190A
double.switch.G602190<- rube(
  model = double.switch.model,
  data = BUGS602190,
  inits = double.switch.inits.602190,
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
summary(double.switch.G602190)

## Gazelle 602191A
double.switch.G602191<- rube(
  model = double.switch.model,
  data = BUGS602191,
  inits = double.switch.inits.602191,
  parameters.to.save = double.switch.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  #n.sims = ns,
  debug = F,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,      
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")
summary(double.switch.G602191)

## Gazelle 618669B
double.switch.G618669B <- rube(
  model = double.switch.model,
  data = BUGS618669B,
  inits = double.switch.inits.618669,
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
summary(double.switch.G618669B)

## Gazelle 618663B
double.switch.G618663B <- rube(
  model = double.switch.model,
  data = BUGS618663B,
  inits = double.switch.inits.618663,
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
summary(double.switch.G618663B)

### Derive states. Insert specific individuals and change columns accordingly
dim(double.switch.G602191$sims.matrix)
idx.sample<-double.switch.G602191$sims.matrix[,11:695] ## Create matrix containing the samples in which each column is an idx parameter
idx.states<-idx.estimate(idx.sample=idx.sample, ncols=685,nrows=4500) ## Posterior summary for idx to gain specific states

### Create dataframe of data and idx for further investigation. Insert specific individuals!
data.602191.idxed<-as.data.frame(cbind(ltraj.602191A$dist/1000, ltraj.602191A$rel.angle,as.Date(ltraj.602191A$date),idx.states))
colnames(data.602191.idxed)<-c("steps", "angles","date", "idx")


### Save results and objects. Insert specific individuals!
setwd("C:/Users/Janni/Documents/1_Masterarbeit/Ergebnisse/Final_results")
## Write summary as csv
write.csv(double.switch.G602189$summary, "602189A_ds_non_zero_modal.csv")
## Save winbugs/rube object in working directory
save(double.switch.G602190, file="602190A_ds_non_zero_modal.RData")


#double.switch.G602191->double.switch.G602191.winter
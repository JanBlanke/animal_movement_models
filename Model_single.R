library(rube)

(wd<-setwd("C:/WinBUGS_Analyse")) # Use short path to make WinBUGS happy

##########################
###### Single model ######
##########################

### Define single model according to Ecology journal
single.model="model{

  ##PRIORS
   b ~ dgamma(0.01,0.01)
   #b ~ dunif(0.25, 20) 	## shape parameter for step length distribution
   #a ~ dunif(0.01, 10)	## rate parameter for step length distribution
   a ~ dgamma(0.01,0.01)

   rho ~ dunif(0,1)		## mean cosine of turns   
   mu ~ dunif(-3.14159265359, 3.14159265359)	## mean direction of turns
   #mu ~ dunif(0, 6.283185)
    
  
   Pi <- 3.14159265359		## define pi


	##MODEL
   for (t in 1:NObs) {

      ## likelihood for steps
      Length[t] ~ dgamma(b, a)	# Modification: Gamma instead  of Weibull distriution for step length

      ## likelihood for turns.  
      ## use the ones trick (see WinBUGS manual) to sample from the Wrapped Cauchy distribution

      ones[t] <- 1
      ones[t] ~ dbern(wC[t])
      ## below is the pdf for Wrapped Cauchy distribution, divided by 500 (arbitrary) to ensure that wC[t] will be less than one
      wC[t] <- ((1/(2*Pi))*(1-rho*rho)/(1+rho*rho-2*rho*cos(Angle[t]-mu)))/500
      
      Angle[t]~ dunif(-3.14159265359,3.14159265359) # For imputing NAs

  }
}"


##########################################
###### Simulate data and test model ######
##########################################

###simulation of single step length and turning angle
simstep<-rgamma(160, shape=0.543, rate = 1.028)
simangle<-rwrpcauchy(160, pi,0.164)


### Bundle data
simBUGS <- list(NObs = length(simstep),
                Length = simstep, 
                Angle = simangle)

### Function to generate starting values
single.inits <- function(){
  list( b = runif(n=1, min=0.25, max=10),
        a = runif(n=1, min=0.1, max=10),
        rho = runif(n=1, min=0.0001, max=0.9999)
        #mu = runif(n=1, min=-3.14159265359, max=3.14159265359)
        ) 
}

### Function to generate strict and fixed starting values
single.inits <- function(){
  list( shape = 0.5,
        rate = 1,
        rho = 0.164,
        mu=3)  
}

### Parameters to be monitored (= to estimate)
single.params <- c("b","a", #the two parameters for the Weibull distribution
                   "rho", "mu") #the two values for the Wrapped Cauchy

### MCMC settings
nc <- 3      		# Number of chains
ni <- 20000				# Number of draws from posterior (for each chain)
nb <- 5000					# Number of draws to discard as burn-in
nt <- 10					# Thinning rate
ns <- 1000


### test with rube before run
summary(rube(model = single.model, data = simBUGS, inits = single.inits))

### Run WinBUGS
single.sim <- rube(
  model = single.model,
  data = simBUGS,
  inits = single.inits,
  parameters.to.save = single.params,
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

summary(single.sim)

##############################
###### Run gazelle data ######
##############################

# Objects to use
ltraj.602189A<-ltraj.602189A[1:703,]
#ltraj.602189A<-na.omit(ltraj.602189A)

ltraj.602190A
ltraj.602190A<-ltraj.602190A[1:687,]
ltraj.602190A<-na.omit(ltraj.602190A)

ltraj.602191A<-ltraj.602191A[1:685,]

ltraj.618669B<-ltraj.618669B[1:249,]

ltraj.618663B<-ltraj.618663B[1:197,]

### Bundle data
G602189BUGS <- list(NObs = nrow(ltraj.602189A),
                    Length = ltraj.602189A$dist/1000, 
                    Angle = ltraj.602189A$rel.angle)

G602190BUGS <- list(NObs = nrow(ltraj.602190A),
                   Length = ltraj.602190A$dist/1000, 
                   Angle = ltraj.602190A$rel.angle)

G602191BUGS <- list(NObs = nrow(ltraj.602191A),
                    Length = ltraj.602191A$dist/1000, 
                    Angle = ltraj.602191A$rel.angle)

G618669BBUGS <- list(NObs = nrow(ltraj.618669B),
                    Length = ltraj.618669B$dist/1000, 
                    Angle = ltraj.618669B$rel.angle)

G618663BBUGS <- list(NObs = nrow(ltraj.618663B),
                     Length = ltraj.618663B$dist/1000, 
                     Angle = ltraj.618663B$rel.angle)



### Function to generate starting values
single.inits <- function(){
  list( b = runif(n=1, min=0.25, max=10),
        a = runif(n=1, min=0.1, max=10),
        rho = runif(n=1, min=0.01, max=0.9),
        mu = runif(n=1, min=-2, max=2)
         ) 
} 

### Parameters to be monitored (= to estimate)
single.params <- c("b","a", #the two parameters for the gamma distribution
                   "rho", "mu") #the two values for the Wrapped Cauchy

### MCMC settings
nc <- 3      		# Number of chains
ni <- 20000			# Number of draws from posterior (for each chain)
nb <- 5000			# Number of draws to discard as burn-in
nt <- 10				# Thinning rate
ns <- 1000      # Number of simulations

### Run WinBUGS
summary(rube(model = single.model,data = G602189BUGS, inits = single.inits))
summary(rube(model = single.model,data = G602190BUGS, inits = single.inits))
summary(rube(model = single.model,data = G602191BUGS, inits = single.inits))
summary(rube(model = single.model,data = G618669BBUGS, inits = single.inits))
summary(rube(model = single.model,data = G618663BBUGS, inits = single.inits))



single.602189<- rube(
    model = single.model,
    data = G602189BUGS,
    inits = single.inits,
    parameters.to.save = single.params,
    n.thin = nt,
    n.chains = nc,
    n.iter = ni,
    n.burnin = nb,
    n.sims = ns,
    debug = TRUE,
    digits = 5,
    wd = getwd(),
    DIC = TRUE,  		
    program = "winbugs",
    bugs.directory ="C:/Program Files/WinBUGS14")
summary(single.602189)

single.602190<- rube(
  model = single.model,
  data = G602190BUGS,
  inits = single.inits,
  parameters.to.save = single.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,    	
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")
summary(single.602190)

single.602191<- rube(
  model = single.model,
  data = G602191BUGS,
  inits = single.inits,
  parameters.to.save = single.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,      
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")
summary(single.602191)

single.618669B<- rube(
  model = single.model,
  data = G618669BBUGS,
  inits = single.inits,
  parameters.to.save = single.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,      
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")
summary(single.618669B)

single.618663B<- rube(
  model = single.model,
  data = G618663BBUGS,
  inits = single.inits,
  parameters.to.save = single.params,
  n.thin = nt,
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.sims = ns,
  debug = TRUE,
  digits = 5,
  wd = getwd(),
  DIC = TRUE,      
  program = "winbugs",
  bugs.directory ="C:/Program Files/WinBUGS14")
summary(single.618663B)



### Plot probability density function of estimated parameters on histogram of gazelle data
hist(ltraj618669B$dist/1000, breaks=100, freq=F)
curve(dgamma(x, shape=mean(single.618669B.nona$sims.list$b), 
             rate=mean(single.618669B.nona$sims.list$a)), type = "l", col="red", add=T, lwd=2)

hist(ltraj618669B$rel.angle, breaks=100, freq=F)
curve(dwrpcauchy(x, rho=mean(single.618669B.nona$sims.list$rho), 
                 mu=mean(single.618669B.nona$sims.list$mu)), type = "l", col="blue", add=T, lwd=2)


### Save objects and summary
setwd("C:/Users/Janni/Desktop")

save(single.618663B, file="single_model_618663.RData")

write.csv(single.618669B$summary, "single_model_618669B_summary.csv")






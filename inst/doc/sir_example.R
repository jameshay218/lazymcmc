## ----setup, include=FALSE------------------------------------------------
library(lazymcmc)
knitr::opts_chunk$set(echo = TRUE, fig.width=7,fig.height=6,
                      message=FALSE, results=FALSE,warning=FALSE)

## ----SIR-----------------------------------------------------------------
## Define the set of ODEs for the model. Return such that we can solve with deSolve
SIR_odes <- function(t, x, params) {
    S <- x[1]
    I <- x[2]
    R <- x[3]
    beta <- params[1]
    gamma <- params[2]
    dS <- -beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I
    list(c(dS,dI,dR))
}
## Make an R0 of something flu like
R0 <- 1.5

## People recover after 5 days, so recovery rate is 1/5
gamma <- 1/5

## So contact rate for this would be beta = gamma*R0
beta <- R0*gamma

## Times to solve model over
t <- seq(0,365,by=0.1)

## Note starting conditions for population size - we're working per capita here
results <- as.data.frame(deSolve::ode(y=c(S=1,I=0.0001,R=0),
                             times=t, func=SIR_odes,
                             parms=c(beta,gamma)))

## Plot the results
plot(results$I~results$time,type='l',col="red",
     ylim=c(0,1),xlab="Time (days)", ylab="Per capita",
     main="SIR dynamics")
lines(results$R~results$time,col="blue")
lines(results$S~results$time,col="green")
legend("topright", title="Key",
       legend=c("Susceptible","Infected","Recovered"),
       col=c("green","red","blue"),
       lty=c(1,1,1))

## ----data_generation-----------------------------------------------------
## Starting population size
S0 <- 999
I0 <- 1
R0 <- 0

## Changing variable name because of clash with recovereds
basic_repro_number <- 1.5
## So contact rate for this would be beta = gamma*basic_repro_number/N
beta <- basic_repro_number*gamma/sum(S0, I0, R0)

## Times to solve model over
t <- seq(1,52*7,by=1)

## Note starting conditions for population size - we're working per capita here
results <- as.data.frame(deSolve::ode(y=c(S=S0,I=I0,R=R0),
                             times=t, func=SIR_odes,
                             parms=c(beta,gamma)))

## Real data is discrete and has errors. 
## Let's generate this from the infected population:
incidence <- c(0,-diff(results$S))
## Get weekly incidence
weekly_incidence <- colSums(matrix(incidence,nrow=7))
## 
real_data <- data.frame(times=seq(1,52*7,by=7),inc=rpois(length(weekly_incidence),weekly_incidence))

## Plot these data with the observations (dots) overlaid
plot(real_data,type='l',col="red",
     xlab="Time (days)", ylab="Number of new infected people",
     main="Observed weekly incidence")


## ----likelihood----------------------------------------------------------
dat <- real_data
## Note that S0, I0, R0, t and SIR_odes are declared outside of
## this function, but are still accessible in the wider R
## environment
likelihood_func <- function(pars){
  beta <- pars[1]
  gamma <- pars[2]
  results <- as.data.frame(deSolve::ode(y=c(S=S0,I=I0,R=R0),
                             times=t, func=SIR_odes,
                             parms=c(beta,gamma)))
  incidence <- c(0,-diff(results$S))
  
  ## Get weekly incidence
  weekly_incidence <- colSums(matrix(incidence,nrow=7))
  
  lik <- sum(dpois(dat$inc,weekly_incidence,log=TRUE))
  lik
}
real_pars <- c(beta,gamma)

## Likelihood of true pars
print(likelihood_func(real_pars))

## ----lazymcmc, results='hide'--------------------------------------------
parTab <- data.frame(names=c("R0","gamma"),
                     values=c(1.5,0.2),
                     fixed=c(0,0),
                     steps=c(0.1,0.1),
                     lower_bound=c(1,0),
                     upper_bound=c(10,1))

mcmcPars <- c("iterations"=2000,"popt"=0.44,"opt_freq"=100,
              "thin"=1,"adaptive_period"=1000,"save_block"=100)


## Putting model solving code in a function for later use
solve_model <- function(pars, t, S0, I0, R0, SIR_odes){
  N <- sum(S0, I0, R0)
  basic_repro_number <- pars["R0"]
  gamma <- pars["gamma"]
  beta <- basic_repro_number*gamma/N
  results <- as.data.frame(deSolve::ode(y=c(S=S0,I=I0,R=R0),
                               times=t, func=SIR_odes,
                               parms=c(beta,gamma)))
  incidence <- c(0,-diff(results$S))
  return(incidence)
}

## We need to put our likleihood function in a closure environment
## It's important to write the function in this form!
create_lik <- function(parTab, data, PRIOR_FUNC,...){
  par_names <- parTab$names
  
  ## Extract observed incidence
  inc <- data$inc
  
  ## Get times to solve model over
  tstep <- 1
  max_weeks <- 52
  t <- seq(tstep,max_weeks*7,by=tstep)

  ## We will pass S0, I0, R0 and SIR_odes 
  N <- S0 + I0 + R0
  ## using the `...` bit.
  likelihood_func <- function(pars){
    names(pars) <- par_names
    incidence <- solve_model(pars, t, S0, I0, R0, SIR_odes)
    ## Get weekly incidence
    weekly_incidence <- colSums(matrix(incidence,nrow=7/tstep))
    
    lik <- sum(dpois(x=inc,lambda=weekly_incidence,log=TRUE))
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
}

## Starting points, chosen pseudo at random
## seeding chains for SIR models is hard with deSolve,
## so I've chosen points near the true values.
startTab <- parTab
startTab$values <- c(2.3, 0.17)

output <- run_MCMC(parTab=startTab, data=dat, mcmcPars=mcmcPars, filename="SIR_fitting",
                   CREATE_POSTERIOR_FUNC = create_lik, mvrPars = NULL, PRIOR_FUNC=NULL,
                   S0=999,I0=1,R0=0,SIR_odes=SIR_odes)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[,c("R0","gamma")]))

## ----correlation---------------------------------------------------------
plot(chain$R0~chain$gamma)

## ----multivariate, results='hide'----------------------------------------
## Prior function
my_prior <- function(pars){
  ## Diffuse gaussian priors on both parameters
  r0_prior <- dnorm(pars[1],1.5,100,1)
  gamma_prior <- dnorm(pars[2],0.2,100,1)
  return(r0_prior + gamma_prior)
}

## Use the previous chain to get a good starting
## covariance matrix
startTab$values <- get_best_pars(chain)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

## Remove restrictions on parameter bounds -
## not important with the multivariate sampler
startTab$lower_bound <- -Inf
startTab$upper_bound <- Inf

## Run with mvrPars
mcmcPars <- c("iterations"=20000,"popt"=0.234,"opt_freq"=500,
              "thin"=1,"adaptive_period"=7500,"save_block"=100)

output2 <- run_MCMC(parTab=startTab, data=dat, mcmcPars=mcmcPars, filename="SIR_fitting",
                   CREATE_POSTERIOR_FUNC = create_lik, mvrPars = mvrPars, PRIOR_FUNC=my_prior,
                   S0=999,I0=1,R0=0,SIR_odes=SIR_odes)
chain2 <- read.csv(output2$file)
chain2 <- chain2[chain2$sampno >= mcmcPars["adaptive_period"],]
plot(coda::as.mcmc(chain2[,c("R0","gamma")]))

## ----final---------------------------------------------------------------
## Get the maximum likelihood parameters and estimate
## the model trajectory for each week
best_pars <- get_best_pars(chain2)
times <- seq(1,52*7,by=1)
best_trajectory <- solve_model(best_pars,times,S0, I0, R0, SIR_odes)
best_incidence <- colSums(matrix(best_trajectory,nrow=7))
best_dat <- data.frame(t=seq(1,52,by=1), y=best_incidence)

## Bit of code to generate prediction intervals
n_samps <- 100
trajectories <- matrix(nrow=n_samps,ncol=52)
samps <- sample(nrow(chain2),n_samps)
for(i in 1:n_samps){
  pars <- as.numeric(chain2[samps[i],c("R0","gamma")])
  names(pars) <- c("R0","gamma")
  tmp <- solve_model(pars, times, S0, I0, R0, SIR_odes)
  trajectories[i,] <- colSums(matrix(tmp, nrow=7))
}
bounds <- apply(trajectories,2,function(x) quantile(x, c(0.025,0.975)))

plot(dat$inc~seq(1,52,by=1),col="red",pch=16, 
     xlab="Time (weeks)",ylab="Incidence",ylim=c(0,100))
lines(best_dat,col="blue")
lines(bounds[1,],col="blue",lty=2,lwd=0.5)
lines(bounds[2,],col="blue",lty=2,lwd=0.5)
legend("topright",c("Data","Model","95% credible intervals"),
       col=c("red","blue","blue"),lty=c(0,1,2),pch=c(16,NA,NA))


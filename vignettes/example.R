mean <- 5
sd <- 2
data <- rnorm(10000, mean=mean, sd=sd)

parTab <- data.frame(values=c(5,2),
                     names=c("mu","sd"),
                     fixed=0,
                     lower_bound=c(0,0),
                     upper_bound=c(10,10),
                     steps=c(0.1,0.1))

my_creation_function <- function(parTab, data, PRIOR_FUNC, ...){
  ##############################
  ## This is where you would manipulate all
  ## of your model arguments. For example,
  ## unpackaging the stuff from `...`,
  ## or transforming model parameters
  ##############################
  
  ## Somewhat redundant example
  parameter_names <- parTab$names
  
  ##############################
  ## This is where you would put your own model code
  ##############################
  f <- function(pars){
    names(pars) <- parameter_names
    mu <- pars["mu"]
    sd <- pars["sd"]
    
    ## Note the use of closures to capture the `data` argument
    lik <- sum(dnorm(data, mu, sd, TRUE))
    if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(pars)
    lik
  }
  return(f)
}

my_prior <- function(pars){
  a <- dnorm(pars["mu"],5,10,1)
  b <- dnorm(pars["sd"],2,10,1)
  return(a + b)
}

mcmcPars <- list("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=5000,"save_block"=100,"temperature" = c(1,5),
              "parallel_tempering_iter" = 10)

n_row_covMat <- sum(parTab$fixed == 0)
covMat <- diag(nrow(parTab))
mvrPars <- list(covMat,2.38/sqrt(n_row_covMat),w=0.8)
mvrPars <- rep(list(mvrPars), length(mcmcPars[["temperature"]]))

startTab <- parTab
startTab$values <- c(3.5,1)
startTab <- rep(list(startTab),length(mcmcPars[["temperature"]]))
startTab[[length(mcmcPars[["temperature"]])]]$values <- c(3,4)

# output <- run_MCMC(parTab=startTab, data=data, mcmcPars=mcmcPars, filename="test", 
#                    CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=mvrPars, 
#                    PRIOR_FUNC = my_prior, OPT_TUNING=0.2)

n_replicates <- 3
startTab <- rep(list(startTab),n_replicates)
mvrPars <- rep(list(mvrPars), n_replicates)
output <- run_MCMC_loop(startTab =startTab, data=data, mcmcPars=mcmcPars, filename=paste0("test",1:3), 
                   CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=mvrPars, 
                   PRIOR_FUNC = my_prior)

# lazymcmc

[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

> MCMC for the lazy (forked from [jameshay218/lazymcmc](https://github.com/jameshay218/lazymcmc))

A package for very generic MCMC in R. Using either univariate uniform or multivariate gaussian proposals, explores the posterior distribution for an arbritary set of model parameters specified in R. The focus is to allow an intuitive, flexible interface rather than fancy on fancy optimisations and syntax.

**The documentation and examples are a work in progress and the code is not fully tested**

## Motivation
There are a lot of great MCMC packages out there (see eg. [rstan](http://mc-stan.org) or the various BUGS iterations). However, sometimes you've written a basic model in R and just want something quick and dirty to fit to some data. The aim of this package is to let you do just that. This package allows you to run an MCMC chain with a lot of flexibility, but without much focus on efficiency (though it is pretty efficient!). The code itself is fast and the adaptive multivariate sample should work for most jobs, but it's not going to be as efficient (in terms of convergence speed and proposals) as something more complex like Hamiltonian MCMC. 

If you have a some R code that, given a set of parameters, generates a model and calculates a likelihood of observing some data, then this package will estimate your posterior distribution.

Here are the package features that you might find useful:
- Simple parameter input allows you to easily specify which parameters you want to estimate and which are fixed, as well as lower and upper bounds (ie. uniform priors)
- Specifying the likelihood function using closures means that your model function can be as complicated as you like
- Flexibility with MCMC parameters, such as the adaptive period, amount of thinning and desired acceptance rate
- MCMC iterations are saved to a .csv file for later analysis. You can specify the frequency of writing to disk
- This is a very minimal R package, so it's great for deploying on things like the [DIDE cluster](https://github.com/dide-tools/didehpc).

## To do
Some features that might be useful to add. Happy to accept requests for more:
- Allow block updates of model parameters (at the moment, it's either Gibbs sampling or a multivariate proposal for the entire parameter vector)
- Add a slice sampler for those pesky multimodal distributions
- More extensive input checks
- Support with `foreach` to run multiple chains at once

## Installation
Installation is straightforward:
```r
devtools::install_github("ada-w-yan/lazymcmc")
```
The outputs of the MCMC work well with the `coda` package. The only real dependency is the `MASS` package for the multivariate normal proposals. If you only want to use the univariate sampler, then you don't need to worry about this.

## Basic example
Perhaps a more relevant example is fitting an SIR model, which can be found in the vignettes [here](https://jameshay218.github.io/lazymcmc/inst/doc/sir_example.html). Below is a very basic example fitting a gaussian distribution to some toy data. It looks like a lot of text, but if you read through it and follow the steps exactly the same logic can be applied to any model.

Firstly, make sure you have the package installed etc.
```r
if(!require(lazymcmc)) devtools::install_github("ada-w-yan/lazymcmc")
library(lazymcmc)
```

The first step is to specify your model. Let's start with something simple - our model is just a gaussian distribution that we draw random numbers from (ie. `rnorm`). We'll generate some random data from it, and then try to re-estimate the mean and standard deviation
```r
mean <- 5
sd <- 2
data <- rnorm(10000, mean=mean, sd=sd)
```

Here's the first key part of the package syntax - the `parTab` structure. This is just a data frame that specifies limitations of the model parameters as follows:
```r
## Generate a table of parameters (either fixed or to be estimated)

## in the parTab data frame for n parameters:
## parTab$values: numeric vector of length n: starting values for parameters
## parTab$names: character vector of length n: parameter names
## parTab$fixed: numeric vector of length n: 1 = fix this parameter, 0 = fit this parameter
## parTab$lower_bound: numeric vector of length n: do not allow proposed parameters to go below this value
## parTab$upper_bound: numeric vector of length n: do not allow proposed parameters to go above this value
## parTab$steps: numeric vector of length n: initial step size for each parameter
## (to be adapted during adaptive period)
parTab <- data.frame(values=c(5,2),
                     names=c("mu","sd"),
                     fixed=0,
                     lower_bound=c(0,0),
                     upper_bound=c(10,10),
                     steps=c(0.1,0.1))
                    
## note that if PRIOR_FUNC == NULL, the prior is assumed to be a uniform
## distribution with these lower and upper bounds
## can also have model-specific columns to be passed to CREATE_POSTERIOR_FUNC
```
All we've done is create a data frame specifying that we have two parameters (`mu` and `sd`); that we want to estimate both of these (`fixed=0`); that both are restricted to the range -Inf to Inf; and that the initial step size of the univariate sampler will be 0.1.

The next step is the most important bit for the user. You'll need to specify a function in `lazymcmc` syntax. Don't worry - there's nothing special about this. All you need to do is wrap whatever your model generating code is (here, it's just `rnorm`) into another function, such that the only thing your model code needs to calculate a likelihood is the clean vector of model parameters (ie. `parTab$values`). The idea is to use [closures](https://www.r-bloggers.com/closures-in-r-a-useful-abstraction/) to return this simplified function, which is what the sampler will use to calculate likelihoods. Thanks to [Rich](https://github.com/richfitz) for showing me this idea. Here's an example:

```r
## You MUST specify the arguments to this function as parTab, data then PRIOR_FUNC. 
## Use the `...` to pass additional arguments through.
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

## Note that we've given it a prior function too. lazymcmc will 
## combine this function pointer with the likelihood function.
my_prior <- function(pars){
  a <- dnorm(pars["mu"],5,10,1)
  b <- dnorm(pars["sd"],2,10,1)
  return(a + b)
}

## To test that it's working
posterior <- my_creation_function(parTab, data, my_prior)
print(posterior(parTab$values))
```
Note that you could just as easily have set `PRIOR_FUNC <- NULL`, and changed the line in `my_creation_function` to just return the likelihood. If no priors are specified, then the uniform prior of the sampling restrictions will apply from `parTab`. We now have everything we need to run the MCMC sampler:

```r
## Update the proposal step size every 1000 iterations (opt_freq) for the first 5000 iterations 
## (adaptive_period) to aim for an acceptance rate of 0.44 (popt). After the adaptive period, 
## run for a further 10000 (iterations) steps. Save every nth rows, where n is "thin" (ie. 1 here).
## Write to disk every 100 iterations (post thinning). Note that the adaptive period is also saved
## to disk
mcmcPars <- c("iterations"=10000,"popt"=0.44,"opt_freq"=1000,
              "thin"=1,"adaptive_period"=5000,"save_block"=100)

## The MCMC code uses the parameter table. Here, we should specify some random starting
## points in the "values" column.
startTab <- parTab
startTab$values <- c(3.5,1)

## You could have done something like this:
## startTab$values <- runif(nrow(startTab), startTab$lower_bound, startTab$upper_bound)

output <- run_MCMC(parTab=startTab, data=data, mcmcPars=mcmcPars, filename="test", 
                   CREATE_POSTERIOR_FUNC=my_creation_function, mvrPars=NULL, 
                   PRIOR_FUNC = my_prior, OPT_TUNING=0.2)

# plot results (exclude adaptive period)
chain <- read.csv(output$file)
plot(coda::as.mcmc(chain[chain$sampno > mcmcPars["adaptive_period"],]))
```
Because `mvrPars=NULL`, we used the univariate sampler. To use the multivariate sampler, we need to specify an additional list of arguments for the multivariate gaussian proposal function:

```r
## A built in function to find the iteration from the MCMC chain
## with best support (maximum likelihood)
bestPars <- get_best_pars(chain)

## Calculate the covariance matrix of the model parameters from the 
## previous MCMC chain. Need to exclude the first and last column
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],2:(ncol(chain)-1)]
covMat <- cov(chain)

## List of multivariate proposal parameters.
## Specifying 3 things; 1) the nxn covariance matrix of all model parameters (including fixed ones);
## 2) the starting step scaler for the proposals. The equation below should work well as an
## initial guess, but it will get updated as the sampler progresses; 
## 3) the sampler updates the covariance matrix (by opt_freq) but gives some weight to
## the old covariance matrix to preserver markov properties. A weighting of 
## 0.8 means that the new covariance matrix after an adaptive step is 
## `covMat <- (1-w)*oldCovMat + w*newCovMat`
mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)

startTab <- parTab
startTab$values <- bestPars
output2 <- run_MCMC(startTab, data, mcmcPars, filename, my_creation_function, mvrPars, PRIOR_FUNC = my_prior  ,0.2)
```

And that's it. If you've worked through this example and understood each step, then you should be able to write your own functions to slot in. Good luck!

## License

GPL-3 © 

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{run_MCMC}
\alias{run_MCMC}
\title{Adaptive Metropolis-within-Gibbs Random Walk Algorithm.}
\usage{
run_MCMC(
  parTab,
  data = NULL,
  mcmcPars,
  filename,
  CREATE_POSTERIOR_FUNC = NULL,
  mvrPars = NULL,
  PRIOR_FUNC = NULL,
  OPT_TUNING = 0.2,
  ...
)
}
\arguments{
\item{parTab}{the parameter table controlling information such as bounds, initial values etc}

\item{data}{the data frame of data to be fitted}

\item{mcmcPars}{named vector named vector with parameters for the MCMC procedure. iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.}

\item{filename}{the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions}

\item{CREATE_POSTERIOR_FUNC}{pointer to posterior function creator used to calculate a likelihood. See the main example - this should return your likelihood function (that only takes a single vector of parameters as an argument).}

\item{mvrPars}{a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)}

\item{PRIOR_FUNC}{user function of prior for model parameters. Should take values, names and local from param_table}

\item{OPT_TUNING}{constant used to indicate what proportion of the adaptive period should be used to build the covariance matrix, if needed}
}
\value{
a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) the last used covariance matrix; 3) the last used scale size
}
\description{
The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
}

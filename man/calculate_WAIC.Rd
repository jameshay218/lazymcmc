% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyses.R
\name{calculate_WAIC}
\alias{calculate_WAIC}
\title{WAIC calculation}
\usage{
calculate_WAIC(chain, parTab, dat, f)
}
\arguments{
\item{chain}{the MCMC chain to be tested}

\item{parTab}{the parameter table used for this chain}

\item{dat}{the data}

\item{f}{the posterior calculating function}
}
\value{
a single WAIC value
}
\description{
Given an MCMC chain, calculates the WAIC
}

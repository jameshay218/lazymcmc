#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
mvr_proposal <- function(values, fixed, covMat){
    proposed <- values
    proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
    return(proposed)
}

#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size.
#' @param values a vector of the parameters to be explored
#' @param lower_bounds a vector of the low allowable bounds for the proposal
#' @param upper_bounds a vector of the upper allowable bounds for the proposal
#' @param steps a vector of step sizes for the proposal
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
univ_proposal <- function(values, lower_bounds, upper_bounds,steps, index){
    mn <- lower_bounds[index]
    mx <- upper_bounds[index]

    rtn <- values

    x <- toUnitScale(values[index],mn,mx)

    ## 5th index is step size
    stp <- steps[index]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    ## Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    ## Cyclical boundary conditions
    ##if (x < 0) x <- 1 + x
    ##if (x > 1) x <- x - 1

    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}


#' MCMC proposal function - univariate normal
#'
#' Proposal function for MCMC random walk, drawing from normal distribution centered on the current value
#' @param values a vector of the parameters to be explored
#' @param sds a vector of standard deviations of the proposal function
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
univ_proposal_normal <- function(values, sds, index){
    rtn <- values
    rtn[index] <- rnorm(1, values[index], sds[index])
    rtn
}

#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @export
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    #if(step > 1) step <- 1
    return(step)
}


#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -100000.
#' @param f the function to be protected
#' @return the protected function
#' @export
protect <- function(f){
    function(...){
        tryCatch(f(...),error=function(e){
            message("caught error: ", e$message)
            list.out <- list(lik = -1000000000, misc = NA)
            list.out
        })
    }
}

#' Get maximum likelihood parameters
#'
#' From an MCMC chain produced by \code{\link{run_MCMC}}, find the row with the highest log likelihood. Return the parameters (excluding the first and last columns, which are sampno and lnlike respectively) that give the highest likelihood.
#' @param chain the MCMC chain as saved by \code{\link{run_MCMC}}
#' @return a named vector of model parameters
#' @export
get_best_pars <- function(chain){
    tmpNames <- colnames(chain)[2:(ncol(chain) - 1)]
    bestPars <- as.numeric(chain[which.max(chain[, "lnlike"]),
                                 2:(ncol(chain) - 1)])
    names(bestPars) <- tmpNames
    return(bestPars)
}
#' Get parameters based on index
#'
#' From an MCMC chain produced by \code{\link{run_MCMC}}, return the parameters from the specified sample number (excluding the first and last columns, which are sampno and lnlike respectively).
#' @param chain the MCMC chain as saved by \code{\link{run_MCMC}}
#' @param index value for sampno, matching the sampno column from the MCMC chain, to be returned
#' @return a named vector of model parameters
#' @export
get_index_pars <- function(chain, index) {
    tmp_names <- colnames(chain)[2:(ncol(chain) - 1)]
    pars <- as.numeric(chain[chain$sampno == index, 2:(ncol(chain) - 1)])
    names(pars) <- tmp_names
    return(pars)
}

toUnitScale <- function(x, min, max){
    return((x-min)/(max-min))
}
fromUnitScale <- function(x,min,max){
    return(min + (max-min)*x)
}

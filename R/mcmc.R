#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{univ_proposal}} or \code{\link{mvr_proposal}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @param data the data frame of data to be fitted
#' @param mcmcPars named vector named vector with parameters for the MCMC procedure. iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param CREATE_POSTERIOR_FUNC pointer to posterior function creator used to calculate a likelihood. See the main example - this should return your likelihood function (that only takes a single vector of parameters as an argument).
#' @param mvrPars a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param PRIOR_FUNC user function of prior for model parameters. Should take values, names and local from param_table
#' @param OPT_TUNING constant used to indicate what proportion of the adaptive period should be used to build the covariance matrix, if needed
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) the last used covariance matrix; 3) the last used scale size
#' @export
run_MCMC <- function(parTab,
                     data=NULL,
                     mcmcPars,
                     filename,
                     CREATE_POSTERIOR_FUNC=NULL,
                     mvrPars=NULL,
                     PRIOR_FUNC=NULL,
                     OPT_TUNING=0.2,
                     ...){
    ## check that input parameters are correctly formatted
    parTab_check <- lazymcmc::param_table_check(parTab)
    if(parTab_check[[1]] == TRUE) return(parTab_check[[2]])
    mcmcPar_check <- lazymcmc::mcmc_param_check(mcmcPars, mvrPars)
    if(mcmcPar_check[[1]] == TRUE) return(mcmcPar_check[[2]])

    ## Allowable error in scale tuning
    TUNING_ERROR <- 0.1
    
    ## Extract MCMC parameters
    iterations <- mcmcPars["iterations"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    adaptive_period_init <- adaptive_period
    save_block <- mcmcPars["save_block"]

    ## added functionality by ada-w-yan: adjusting adaptive period depending on 
    ## the acceptance ratio.  If acceptance ratio within the last opt_freq 
    ## iterations of adaptive period not close enough to optimal,
    ## extend adaptive period for another adaptive_period iterations.
    ## To avoid the algorithm running infinitely, an absolute upper bound of
    ## max_adaptive_period adaptive iterations is passed in.
    ## if mcmcPars does not contain an element named adaptiveLeeway,
    ## never extend the adaptive period.
    ## 'close enough' to optimal means within poptRange as specified below.
    ## currently only works for univariate proposals.
    if(match("adaptiveLeeway",names(mcmcPars)) {
      adaptiveLeeway <- mcmcPars["adaptiveLeeway"]
      max_adaptive_period <- ["max_adaptive_period"]
    } else {
      adaptiveLeeway  <- 0
      max_adaptive_period <- adaptive_period
    }
    poptRange <- popt * c((1 - adaptiveLeeway),(1 + adaptiveLeeway))
    poptRange <- pmax(0,poptRange)
    poptRange <- pmin(1,poptRange)

    param_length <- nrow(parTab)
    unfixed_pars <- which(parTab$fixed == 0)
    unfixed_par_length <- nrow(parTab[parTab$fixed== 0,])
    current_pars <- parTab$values
    par_names <- as.character(parTab$names)

    ## Parameter constraints
    lower_bounds <- parTab$lower_bound
    upper_bounds <- parTab$upper_bound
    steps <- parTab$steps
    fixed <- parTab$fixed
    
    ## Arrays to store acceptance rates
    ## If univariate proposals
    if(is.null(mvrPars)){
        tempaccepted <- tempiter <- integer(param_length)
        reset <- integer(param_length)
        reset[] <- 0
    } else { # If multivariate proposals
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]][unfixed_pars,unfixed_pars]
        scale <- mvrPars[[2]]
        w <- mvrPars[[3]]
    }

    posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,data, 
                                                    PRIOR_FUNC,...))
    
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")

    ## Create empty chain to store every iteration for the adaptive period
    opt_chain <- matrix(nrow=adaptive_period,ncol=unfixed_par_length)
    chain_index <- 1
    
    ## Initial conditions ------------------------------------------------------
    ## Initial likelihood
    posterior.out <- posterior_simp(current_pars)
    ## added feature by ada-w-yan: for each recorded iteration,
    ## we can now write a vector with miscellaneous output to file in addition
    ## to the parameter values and likelihood
    ## (for example, predicted model output)
    ## usage: posterior_simp(proposal) should either return 
    ## the likelihood as a numeric vector of length 1, 
    ## or a list with elements
    ## list$lik: the likelihood as a numeric vector of length 1
    ## list$misc: any additional output as a vector
    ## the length of list$misc has to be the same for all proposal values
    if(is.atomic(posterior.out)){
      probab <- posterior.out
      misc <- numeric()
    } else {
      probab <- posterior.out$lik
      misc <- unname(posterior.out$misc)
    }
    misc_length <- length(misc)

    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+2+misc_length)

    ## Set up initial csv file
    if(is.atomic(posterior.out) || is.null(names(posterior.out$misc))){
      misc_colnames <- rep("misc",misc_length)
      if(misc_length > 0){
        misc_colnames <- paste0(misc_colnames,1:misc_length)
      }
    } else {
      misc_colnames <- names(posterior.out$misc)
    }
    chain_colnames <- c("sampno",par_names,misc_colnames,"lnlike")
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,current_pars,misc,probab)

    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2
    par_i <- 1
    i <- 1
    pcurUnfixed <- -1
    
    while (i <= (iterations+adaptive_period)){
        ## If using univariate proposals
        if(is.null(mvrPars)) {
            ## For each parameter (Gibbs)
            j <- unfixed_pars[par_i]
            par_i <- par_i + 1
            if(par_i > unfixed_par_length) par_i <- 1
            proposal <- univ_proposal(current_pars, lower_bounds, upper_bounds, steps,j)
            tempiter[j] <- tempiter[j] + 1
            ## If using multivariate proposals
        } else {
            proposal <- mvr_proposal(current_pars, unfixed_pars, scale*covMat)
            tempiter <- tempiter + 1
        }
        ## Propose new parameters and calculate posterior
        ## Check that all proposed parameters are in allowable range
        if(!any(
                proposal[unfixed_pars] < lower_bounds[unfixed_pars] |
                proposal[unfixed_pars] > upper_bounds[unfixed_pars]
            )
           ){
            ## Calculate new likelihood and find difference to old likelihood
            posterior.out <- posterior_simp(proposal)
            if(is.atomic(posterior.out)){
              new_probab <- posterior.out
              new_misc <- numeric()
            } else {
              new_probab <- posterior.out$lik
              new_misc <- posterior.out$misc
            }
            
            log_prob <- min(new_probab-probab,0)
         
            ## Accept with probability 1 if better, or proportional to
            ## difference if not
            if(is.finite(log_prob) && log(runif(1)) < log_prob){
                current_pars <- proposal
                probab <- new_probab
                misc <- new_misc
                
                ## Store acceptances
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else {
                    tempaccepted <- tempaccepted + 1
                }
            }
        }
        
        
        ## If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        if(i %% thin ==0){
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-1-misc_length)] <- current_pars
            if(misc_length > 0){
              save_chain[no_recorded,(ncol(save_chain)-1-misc_length+1):(ncol(save_chain)-1)] <- unname(misc)
            }
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }

       
        
        ## If within adaptive period, need to do some adapting!
        if(i <= adaptive_period){
            ## Current acceptance rate
            pcur <- tempaccepted/tempiter
            ## Save each step
            opt_chain[chain_index,] <- current_pars[unfixed_pars]
           
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                ## If using univariate proposals
                if(is.null(mvrPars)){
                    ## For each non fixed parameter, scale the step size
                    for(x in unfixed_pars) steps[x] <- scaletuning(steps[x],popt,pcur[x])
                    message(cat("Pcur: ", pcur[unfixed_pars],sep="\t"))
                    message(cat("Step sizes: ", steps[unfixed_pars],sep="\t"))
                    tempaccepted <- tempiter <- reset

                } else {       ## If using multivariate proposals
                    if(chain_index > OPT_TUNING*adaptive_period & chain_index < (0.8*adaptive_period)){
                        oldCovMat <- covMat
                        ## Creates a new covariance matrix, but weights it with the old one
                        covMat <- cov(opt_chain[1:chain_index,])
                        covMat <- w*covMat + (1-w)*oldCovMat
                    }
                    ## Scale tuning for last 20% of the adpative period
                    if(chain_index > (0.8)*adaptive_period){
                        scale <- scaletuning(scale, popt,pcur)
                    }
                    tempiter <- tempaccepted <- 0

                    message(cat("Pcur: ", pcur,sep="\t"))
                    message(cat("Scale: ", scale,sep="\t"))
                }
            }
            chain_index <- chain_index + 1
        }
      
      ## added functionality by ada-w-yan
      ## if at end of adaptive period,
      ## decide whether to extend adaptive period
        if(i == adaptive_period){
            ## update current acceptance probability
            pcurUnfixed <- pcur[unfixed_pars]
            ## if current acceptance probability not close enough to optimal, 
            ## extend adaptive period
            if(((max(pcurUnfixed) > poptRange[2]) || (min(pcurUnfixed) < poptRange[1])) &&
               adaptive_period < max_adaptive_period){
              ## update total number of adaptive iterations run
              adaptive_period <- adaptive_period + adaptive_period_init
              ## expand matrix in which adaptive iterations are stored
              opt_chain <- rbind(opt_chain, matrix(nrow=adaptive_period_init,ncol=unfixed_par_length))
            } else {
              p_accept_adaptive <- pcurUnfixed
            }
        }
        
        if(i %% save_block == 0){
            message(cat("Current iteration: ", i, sep="\t"))
            ## Print out optimisation frequencies
        }
        
        if(no_recorded == save_block){
            write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        sampno <- sampno + 1
        i <- i + 1
    }
    
    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }

    if(is.null(mvrPars)){
        covMat <- NULL
        scale <- NULL
    } else {
        steps <- NULL
    }
    
    p_accept <- tempaccepted/tempiter
    p_accept <- p_accept[unfixed_pars]
    
    ## by ada-w-yan: we now output the actual adaptive period used,
    ## the acceptance probability during the last opt_freq iterations
    ## of the adaptive period, and the acceptance probability during the
    ## non-adaptive iterations
    return(list("file"=mcmc_chain_file,"covMat"=covMat,"scale"=scale, 
                "steps"=steps, "adaptive_period" = adaptive_period,
                "p_accept_adaptive" = p_accept_adaptive,
                "p_accept" = p_accept))
}
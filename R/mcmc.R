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
                     seed,
                     ...){

  ## check that input parameters are correctly formatted
  parallel_tempering_flag <- is.list(parTab[[1]])
  if(parallel_tempering_flag){ # if parallel tempering
    
    parTab_check <- lapply(parTab,lazymcmc::param_table_check)
    parTab_check_flag <- unlist(lapply(parTab_check, function(x) x[[1]]))
    if(any(parTab_check_flag)) return(parTab_check[[which(parTab_check_flag)[1]]][[2]])
    
    mcmcPar_check <- lapply(mvrPars, function(x) lazymcmc::mcmc_param_check(mcmcPars, x))
    mcmcPar_check_flag <- unlist(lapply(mcmcPar_check, function(x) x[[1]]))
    if(any(mcmcPar_check_flag)) return(mcmcPar_check[[which(mcmcPar_check_flag)[1]]][[2]])
    if(!("parallel_tempering_iter" %in% names(mcmcPars))) {
      stop("parallel_tempering_iter not given")
    }
    
  } else {
    parTab_check <- lazymcmc::param_table_check(parTab)
    if(parTab_check[[1]] == TRUE) return(parTab_check[[2]])
    
    mcmcPar_check <- lazymcmc::mcmc_param_check(mcmcPars, mvrPars)
    if(mcmcPar_check[[1]] == TRUE) return(mcmcPar_check[[2]])
  }
  
  ## Allowable error in scale tuning
  TUNING_ERROR <- 0.1
  
  ## Extract MCMC parameters
  iterations <- mcmcPars[["iterations"]]
  popt <- mcmcPars[["popt"]]
  opt_freq<- mcmcPars[["opt_freq"]]
  thin <- mcmcPars[["thin"]]
  adaptive_period<- mcmcPars[["adaptive_period"]]
  adaptive_period_init <- adaptive_period
  save_block <- mcmcPars[["save_block"]]
  if("temperature" %in% names(mcmcPars)){
    temperatures <- mcmcPars[["temperature"]]
  } else {
    temperatures <- 1
  }
  
  if("parallel_tempering_iter" %in% names(mcmcPars)){
    parallel_tempering_iter <- mcmcPars[["parallel_tempering_iter"]]
  } else {
    parallel_tempering_iter <- iterations + adaptive_period + 1
  }
  
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
  if("adaptiveLeeway" %in% names(mcmcPars)) {
    adaptiveLeeway <- mcmcPars[["adaptiveLeeway"]]
    max_adaptive_period <- mcmcPars[["max_adaptive_period"]]
  } else {
    adaptiveLeeway  <- 0
    max_adaptive_period <- adaptive_period
  }
  poptRange <- popt * c((1 - adaptiveLeeway),(1 + adaptiveLeeway))
  poptRange <- pmax(0,poptRange)
  poptRange <- pmin(1,poptRange)
  
  if(!parallel_tempering_flag){ # only one set of starting values
    current_pars <- parTab$values
    steps <- list(parTab$steps)
    start_pars <- parTab$values
  } else {
    current_pars <- parTab[[1]]$values
    # store starting values for parallel chains
    start_pars <- lapply(parTab, function(x) x$values)
    steps <- lapply(parTab, function(x) x$steps)
    parTab <- parTab[[1]]
  }
  
  param_length <- nrow(parTab)
  
  unfixed_pars <- which(parTab$fixed == 0)
  
  unfixed_par_length <- nrow(parTab[parTab$fixed== 0,])
  
  par_names <- as.character(parTab$names)
  
  ## Parameter constraints
  lower_bounds <- parTab$lower_bound
  upper_bounds <- parTab$upper_bound
  fixed <- parTab$fixed
  
  ## Arrays to store acceptance rates
  ## If univariate proposals
  if(is.null(mvrPars)){
    tempaccepted <- tempiter <- integer(param_length)
    reset <- integer(param_length)
    reset[] <- 0
  } else { # If multivariate proposals
    tempaccepted <- tempiter <- 0
    if(!parallel_tempering_flag){
      mvrPars <- list(mvrPars)
    }
    covMat <- lapply(mvrPars, function(x) x[[1]][unfixed_pars,unfixed_pars])
    scale <- vapply(mvrPars, function(x) x[[2]], double(1))
    w <- mvrPars[[1]][[3]]
    reset <- 0
  }
  
  posterior_simp <- protect(CREATE_POSTERIOR_FUNC(parTab,data, 
                                                  PRIOR_FUNC,...))
  
  ## Setup MCMC chain file with correct column names
  mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
  
  ## Create empty chain to store every iteration for the adaptive period
  opt_chain <- matrix(nrow=adaptive_period,ncol=unfixed_par_length)
  opt_chain <- rep(list(opt_chain),length(temperatures))
  chain_index <- 1
  
  ## Initial conditions ------------------------------------------------------
  ## Initial likelihood
  posterior_out <- posterior_simp(current_pars)
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
  if(is.atomic(posterior_out)){
    probab <- posterior_out
    misc <- numeric()
  } else {
    probab <- posterior_out$lik
    misc <- unname(posterior_out$misc)
  }
  misc_length <- length(misc)
  
  ## Create empty chain to store "save_block" iterations at a time
  save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=param_length+2+misc_length)
  
  ## Set up initial csv file
  if(is.atomic(posterior_out) || is.null(names(posterior_out$misc))){
    misc_colnames <- rep("misc",misc_length)
    if(misc_length > 0){
      misc_colnames <- paste0(misc_colnames,1:misc_length)
    }
  } else {
    misc_colnames <- names(posterior_out$misc)
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
  p_accept_adaptive <- NA
  offset <- 0
  potential_swaps <- swaps <- 0
  
  ## set seed
  if(!missing(seed)){
    # using an integer
    if(length(seed) == 1){
      set.seed(seed)
      # or a previous state
    } else {
      .Random.seed <- seed
    }
    
  }

  create_run_MCMC_single_iter_fn <- function(unfixed_pars,unfixed_par_length,
                                             lower_bounds,upper_bounds,steps,scale,
                                             covMat,mvrPars,temperature){
    f <- function(par_i,current_pars,misc,probab,tempaccepted,tempiter){
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
        posterior_out <- posterior_simp(proposal)
        if(is.atomic(posterior_out)){
          new_probab <- posterior_out
          new_misc <- numeric()
        } else {
          new_probab <- posterior_out$lik
          new_misc <- posterior_out$misc
        }
        
        log_prob <- min(new_probab-probab,0)
        
        ## Accept with probability 1 if better, or proportional to
        ## difference if not
        if(is.finite(log_prob) && log(runif(1)) < log_prob/temperature){
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
      list("par_i" = par_i, "current_pars" = current_pars,
           "misc" = misc,
           "probab" = probab, "tempaccepted" = tempaccepted,
           "tempiter" = tempiter)
    }
    f
  }

  run_MCMC_single_iter <- lapply(seq_along(temperatures),
                                 function(x) create_run_MCMC_single_iter_fn
                                 (unfixed_pars,unfixed_par_length,
                                 lower_bounds,upper_bounds,
                                   steps[[x]],scale[x],
                                   covMat[[x]],mvrPars,temperatures[x]))
  
  # initialise MCMC
  mcmc_list <- list("par_i" = par_i, "current_pars" = current_pars,
                    "misc" = misc, "probab" = probab, "tempaccepted" = tempaccepted,
                    "tempiter" = tempiter)
  
  # replicate list for parallel tempering
  mcmc_list <- rep(list(mcmc_list),length(temperatures))
  # start values for parallel tempering
  if(parallel_tempering_flag){
    mcmc_list <- Map(function(x,y) modifyList(x,list(current_pars = y)), mcmc_list, start_pars)
  }

  # main body of running MCMC

  while (i <= (iterations+adaptive_period)){

    mcmc_list <- Map(do.call, run_MCMC_single_iter, mcmc_list)
    
    # perform parallel tempering
    
    if(i %% parallel_tempering_iter == 0){
      parallel_tempering_list <- parallel_tempering(mcmc_list, temperatures, offset)
      mcmc_list <- parallel_tempering_list$mcmc_list
      swaps <- swaps + parallel_tempering_list$swaps
      potential_swaps <- potential_swaps + 1
      offset <- 1 - offset
    }
    
    ## If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
    ## save this block of chain to file and reset chain
    if(i %% thin ==0){
      current_pars <- mcmc_list[[1]][["current_pars"]]
      probab <- mcmc_list[[1]][["probab"]]
      misc <- mcmc_list[[1]][["misc"]]
      
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

      ## Save each step
      
      save_opt_chain <- function(opt_chain,mcmc_list,unfixed_pars,chain_index){
        current_pars <- mcmc_list[["current_pars"]]
        opt_chain[chain_index,] <- current_pars[unfixed_pars]
        opt_chain
      }
      
      opt_chain <- Map(function(x,y) save_opt_chain(x,y,unfixed_pars,chain_index),opt_chain,mcmc_list)
      ## If in an adaptive step
      if(chain_index %% opt_freq == 0){
        
        reset_acceptance <- function(mcmc_list, reset){
          mcmc_list[["tempaccepted"]] <- mcmc_list[["tempiter"]] <- reset
          mcmc_list
        }
        
        ## If using univariate proposals
        if(is.null(mvrPars)){
          
          ## Current acceptance rate
          pcur <- lapply(mcmc_list, function(x) x[["tempaccepted"]] / x[["tempiter"]])
          
          ## For each non fixed parameter, scale the step size

          scale_univariate <- function(steps, popt, pcur, unfixed_pars){
            steps[unfixed_pars] <- vapply(unfixed_pars,function(x) scaletuning(steps[x],popt,pcur[x]),
                                          double(1))
            steps
          }

          steps <- Map(function(x,y) scale_univariate(x, popt, y, unfixed_pars), steps, pcur)
          
          message(cat("Pcur: ", pcur[[1]][unfixed_pars],sep="\t"))
          message(cat("Step sizes: ", steps[[1]][unfixed_pars],sep="\t"))
          
          mcmc_list <- lapply(mcmc_list, function(x) reset_acceptance(x, reset))
          
        } else {       ## If using multivariate proposals
          pcur <- vapply(mcmc_list, function(x) x[["tempaccepted"]] / x[["tempiter"]], double(1))
          if(chain_index > OPT_TUNING*adaptive_period & chain_index < (0.8*adaptive_period)){
            
            scale_mvr <- function(opt_chain,covMat,w,chain_index){
              oldCovMat <- covMat
              ## Creates a new covariance matrix, but weights it with the old one
              covMat <- cov(opt_chain[1:chain_index,])
              covMat <- w*covMat + (1-w)*oldCovMat
              covMat
            }
            
            covMat <- Map(function(x,y) scale_mvr(x,y,w,chain_index), opt_chain, covMat)
          }
          ## Scale tuning for last 20% of the adpative period
          if(chain_index > (0.8)*adaptive_period){
            scale <- Map(function(x,y) scaletuning(x, popt,y), scale, pcur)
            scale <- unlist(scale)
          }
          
          mcmc_list <- lapply(mcmc_list, function(x) reset_acceptance(x, reset))
          
          message(cat("Pcur: ", pcur[1],sep="\t"))
          message(cat("Scale: ", scale[1],sep="\t"))
        }
        
        # calibrate temperatures
        swap_ratio <- swaps / potential_swaps
        temperatures <- calibrate_temperatures(temperatures, swap_ratio)
        swaps <- potential_swaps <- 0
      }
      
      chain_index <- chain_index + 1
      
      # remake run_MCMC_single_iter for new step sizes
      run_MCMC_single_iter <- lapply(seq_along(temperatures),
                                     function(x) create_run_MCMC_single_iter_fn
                                     (unfixed_pars,unfixed_par_length,
                                       lower_bounds,upper_bounds,
                                       steps[[x]],scale[x],
                                       covMat[[x]],mvrPars,temperatures[x]))
    }
    
    ## added functionality by ada-w-yan
    ## if at end of adaptive period,
    ## decide whether to extend adaptive period
    if(is.null(mvrPars) && i == adaptive_period){
      ## update current acceptance probability
      pcurUnfixed <- pcur[[1]][unfixed_pars]
      ## if current acceptance probability not close enough to optimal, 
      ## extend adaptive period
      if(((max(pcurUnfixed) > poptRange[2]) || (min(pcurUnfixed) < poptRange[1])) &&
         adaptive_period < max_adaptive_period){
        ## update total number of adaptive iterations run
        adaptive_period <- adaptive_period + adaptive_period_init
        ## expand matrix in which adaptive iterations are stored
        expand_matrix <- matrix(nrow=adaptive_period_init,ncol=unfixed_par_length)
        opt_chain <- lapply(opt_chain, function(x) rbind(x, expand_matrix))
      } else {
        p_accept_adaptive <- pcurUnfixed[[1]]
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
  
  p_accept <- mcmc_list[[1]][["tempaccepted"]] / mcmc_list[[1]][["tempiter"]]
  if(is.null(mvrPars)){
    p_accept <- p_accept[unfixed_pars]
  }
  
  current_pars <- lapply(mcmc_list, function(x)x$current_pars)
  ## by ada-w-yan: we now output the actual adaptive period used,
  ## the acceptance probability during the last opt_freq iterations
  ## of the adaptive period, and the acceptance probability during the
  ## non-adaptive iterations
  return(list("file" = mcmc_chain_file,
              "current_pars" = current_pars,
              "covMat" = covMat,
              "scale" = scale, 
              "steps" = steps, 
              "adaptive_period" = adaptive_period,
              "p_accept_adaptive" = p_accept_adaptive,
              "p_accept" = p_accept,
              "temperatures" = temperatures,
              "seed" = .Random.seed))
}

#' Wrapper to run multiple parallel MCMC chains until convergence
#'
#' Runs run_MCMC, calculates convergence diagnostics, 
#' repeats until convergence or maximum number of iterations reached
#' currently only works for univariate proposals
#' @param startTab the parameter table controlling information 
#' such as bounds, initial values etc.
#' Because the different chains have different initial values, this is a list 
#' where each element is of the form of parTab in run_MCMC, but the
#' values column has different starting values.
#' startTab is actually a list of n of these where n is the number of parallel chains
#' @param data the data frame of data to be fitted
#' @param mcmcPars named vector with parameters for the MCMC procedure. 
#' mandatory: iterations, popt, opt_freq, thin, burnin, adaptive_period, save_block
#' optional: max_total_iterations, max_adaptive_period, adaptiveLeeway
#' @param filenames character vector:
#' the full filepaths at which the parallel MCMC chains should be saved. 
#' "_chain.csv" will be appended to the end of this, 
#' so filenames should have no file extensions
#' the length of filenames will be the number of parallel chains run.
#' @param CREATE_POSTERIOR_FUNC pointer to posterior function creator used to 
#' calculate a likelihood. See the main example - 
#' this should return your likelihood function (that only takes a single vector 
#' of parameters as an argument).
#' @param mvrPars a list of parameters if using a multivariate proposal. 
#' Must contain an initial covariance matrix, weighting for adapting cov matrix,
#'  and an initial scaling parameter (0-1)
#' mvrPars is acutally a list of n of these where n is the number of parallel chains
#' @param PRIOR_FUNC user function of prior for model parameters. 
#' Should take values, names and local from param_table
#' @param run_parallel logical vector of length 1. If TRUE, run in parallel on cluster
#' @return a list with: 1) convergence diagnostics; 2) the output from run_MCMC
#' (the first loop around)
#' @export
run_MCMC_loop <- function(startTab, data, mcmcPars, filenames,  
                          CREATE_POSTERIOR_FUNC, 
                          mvrPars, PRIOR_FUNC, run_parallel = FALSE){
  
  parallel_tempering_flag <- 
    ("temperature" %in% names(mcmcPars) && length(mcmcPars[["temperature"]]) > 1)

  n_replicates <- length(filenames)
  
  if(parallel_tempering_flag){
    n_pars <- nrow(startTab[[1]][[1]])
  } else {
    n_pars <- nrow(startTab[[1]])
  }

  seed <- lapply(seq_len(n_replicates), identity)
  
  diagnostics <- list(converged = FALSE)
  startTab_current <- startTab
  total_iterations <- 0
  filenames_current <- filenames
  if(!("max_total_iterations" %in% names(mcmcPars))){
    mcmcPars <- c(mcmcPars, "max_total_iterations" = mcmcPars[["iterations"]])
  }
  
  thin <- mcmcPars[["thin"]]
  max_total_iterations <- mcmcPars[["max_total_iterations"]]
  iterations <- mcmcPars[["iterations"]]
  
  mcmcPars <- rep(list(mcmcPars), n_replicates)

  timing <- system.time(
    while(!diagnostics$converged && total_iterations < max_total_iterations){

      ## run MCMC for random starting values
      output_current <- parLapply_wrapper(run_parallel,seq_len(n_replicates),
                                  function(x) run_MCMC(startTab_current[[x]], data, mcmcPars[[x]],
                                                       filenames_current[x], CREATE_POSTERIOR_FUNC,
                                                       mvrPars[[x]], PRIOR_FUNC = PRIOR_FUNC,
                                                       0.1, seed = seed[[x]]))

      # if first time running
      if(total_iterations == 0){
        output <- output_current
      } else {
        # append the new output file
        append.csv <- function(x){
          # read new output file
          temp <- data.table::fread(output_current[[x]]$file)
          temp <- temp[seq(2,nrow(temp)),]
          # renumber samples to continue from old file
          temp$sampno <- seq_len(nrow(temp))*thin + (output[[x]]$adaptive_period + total_iterations + 1)
          current_pars <- as.numeric(temp[nrow(temp),seq(2,(n_pars+1))])
          # append to old file
          write.table(temp, output[[x]]$file,
                      row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
          # delete cont file
          file.remove(output_current[[x]]$file)
          invisible(NULL)
        }
        lapply(seq_len(n_replicates),append.csv)
      }
      
      # get current parameters
      current_pars <- lapply(output_current, function(x) x$current_pars)

      # write output e.g. step size to file
      output_write <- output_current
      seed_idx <- which(names(output_write[[1]]) == "seed")
      output_write <- lapply(output_write, function(x) x[-seed_idx])
      dump("output_write", paste0(filenames[1],".output"), append = (total_iterations > 0))

      ## can't calculate convergence diagnostics if only one replicate run,
      ## so stop running here
      if(n_replicates == 1){
        diagnostics$converged <- TRUE
        diagnostics$burn_in <- iterations / 2
      } else {
        if(parallel_tempering_flag){
          fixed <- startTab[[1]][[1]]$fixed
        } else {
          fixed <- startTab[[1]]$fixed
        }
        
        ## calculate convergence diagnostics
        diagnostics <- calc_diagnostics(filenames = vapply(output, function(x) x$file, character(1)),
                                        check_freq = floor(iterations / 10 / thin),
                                        fixed = fixed,
                                        skip = vapply(output, function(x) floor(x$adaptive_period / thin), 
                                                      double(1)))
      }

      total_iterations <- total_iterations + iterations

      ## get things ready to run again if it hasn't converged
      
      make_new_startTab_wrapper <- function(startTab_single){
        f <- function(values, steps){
          if(missing(steps)){
            steps <- startTab_single$steps
          }
          cbind(data.frame("values" = values),
                startTab_single[c("names","fixed","lower_bound","upper_bound")],
                data.frame("steps" = steps))
        }
        f
      }
      
      make_new_mcmcPars <- function(old_mcmcPars, temperature){
        new_mcmcPars <- old_mcmcPars
        new_mcmcPars$temperature <- temperature
        new_mcmcPars$adaptive_period <- 0
        new_mcmcPars
      }
      
      temperatures <- lapply(output_current, function(x) x$temperatures)
      mcmcPars <- Map(make_new_mcmcPars, mcmcPars, temperatures)
    
      if(is.null(mvrPars)){

        steps <- lapply(output_current, function(x) x$steps)
        if(parallel_tempering_flag){
          startTab_single <- startTab[[1]][[1]]

          make_new_startTab <- make_new_startTab_wrapper(startTab_single)
          startTab_current <- Map(function(x,y) Map(make_new_startTab, x, y), 
                                  current_pars, steps)
        } else {
          startTab_single <- startTab[[1]]
          make_new_startTab <- make_new_startTab_wrapper(startTab_single)
          startTab_current <- Map(make_new_startTab,current_pars, steps)
        }

      } else {

        make_new_mvrPars_wrapper <- function(startTab_single, w){
          f <- function(covMat, scale){
            covMat_expand <- diag(nrow(startTab_single))
            unfixed <- which(startTab_single$fixed == 0)
            covMat_expand[unfixed,unfixed] <- covMat
            list(covMat_expand, scale, w = w)
          }
          f
        }
        
        covMat <- lapply(output_current, function(x) x$covMat)
        scale <- lapply(output_current, function(x) x$scale)

        if(parallel_tempering_flag){
          startTab_single <- startTab[[1]][[1]]
          make_new_startTab <- make_new_startTab_wrapper(startTab_single)
          make_new_mvrPars <- make_new_mvrPars_wrapper(startTab_single, mvrPars[[1]][[1]]$w)
          
          startTab_current <- lapply(current_pars,
                                     function(x) lapply(x, make_new_startTab))

          mvrPars <- Map(function(x,y) Map(make_new_mvrPars, x, y), covMat, scale)
        } else {
          startTab_single <- startTab[[1]]
          make_new_startTab <- make_new_startTab_wrapper(startTab_single)
          make_new_mvrPars <- make_new_mvrPars_wrapper(startTab_single, mvrPars[[1]]$w)

          startTab_current <- lapply(current_pars, make_new_startTab)
          mvrPars <- Map(make_new_mvrPars, covMat, scale)
        }

      }

      filenames_current <- paste0(filenames,"_new")
      seed <- lapply(output_current, function(x) x$seed)
    }
  )

  # write elapsed time to file
  write(timing,paste0(filenames[1],".time"))

  # write diagnostics to file
  dump("diagnostics",
       paste0(filenames[1],".diagnostics"),
       append = FALSE)
  list("diagnostics" = diagnostics, "output" = output)
}

#' see whether MCMC chains have converged; if so, calculate effective sample size and burn-in
#'
#' see whether MCMC chains have converged; if so, calculate effective sample size and burn-in
#'
#' @param filenames character vector of filenames, each corresponding to one csv
#' file outputted by run_MCMC
#' @param check_freq: check whether convergence has occurred every check_freq iterations
#' along the chain
#' @param fixed: logical vector: vector whose entries are TRUE if the corresponding
#' parameter is fixed, FALSE if the parameter is fitted. i.e. parTab$fixed
#' @param skip: numeric vector either of length 1 or of same length as filenames:
#' skip this many entries from the start of each csv file
#' @return if convergence has occurred, return a list with the elements
#' converged: logical = TRUE
#' burn_in: the number of burn-in iterations
#' combined_size: effective sample size combined across chains
#' if convergence has not occurred, return a list with the element
#' converged: logical = FALSE
#' max.prsf:maximum value of potential scale reduction factor across parameters
#' @export
calc_diagnostics <- function(filenames,check_freq,fixed,skip = 0){

  ## replicate skip value if more than one filename given but only one skip value given
  if(length(skip) == 1){
    skip <- rep(skip, length(filenames))
  }
  ## check if number of skip values equal to number of filenames
  if(length(skip) != length(filenames)){
    stop("input vector filenames different length to input vector skip")
  }
  data <- lapply(filenames,function(x) data.table::fread(x))
  
  # discard parameters which are fixed
  data <- Map(function(x,y) x[(y+1):nrow(x),2:(length(fixed)+1)], data, skip)
  # data <- lapply(data,function(x) x[(skip+1):nrow(x),2:(length(fixed)+1)])
  
  fitted <- !as.logical(fixed)
  data <- lapply(data,function(x) x[,fitted,drop = FALSE, with = FALSE])

  # determine length of shortest chain
  min_length <- min(vapply(data,function(x)dim(x)[1], double(1)))
  # if prsf for all parameters below threshold, converged
  thres = 1.1
  # calculate potential scale reduction factor
  # note: discards first half of chain as default

  max_psrf <- numeric()

  for (k in seq(check_freq,min_length,check_freq)){

    data_temp <- lapply(data,function(x)x[1:k,])
    data_temp <- lapply(data_temp,mcmc)
    combinedchains <- mcmc.list(data_temp)

    psrf <- gelman.diag(combinedchains)
    max_psrf <- c(max_psrf,max(psrf[[1]][,2]))
    print(max_psrf[length(max_psrf)])
    # if converged, calculate summary statistics at this point and return
    if(max_psrf[length(max_psrf)] < thres){
      burn_in <- ceiling(k/2)
      # keep second half of converged chain, plus all samples afterwards
      data <- lapply(data,function(x)x[(burn_in+1):min_length,])
      data <- lapply(data,mcmc)
      combinedchains <- mcmc.list(data)
      # calculate effective sample size
      combined_size <- effectiveSize(combinedchains)
      return(list("converged" = TRUE,
                  "max_psrf" = max_psrf,
                  "burn_in" = burn_in,
                  "combined_size" = combined_size))
    }
  }
  # else return summary statistics for second half of chain
  burn_in <- ceiling(k/2)
  # keep second half of converged chain, plus all samples afterwards
  data <- lapply(data,function(x)x[(burn_in+1):min_length,])
  data <- lapply(data,mcmc)
  combinedchains <- mcmc.list(data)
  # calculate effective sample size
  combined_size <- effectiveSize(combinedchains)
  list("converged" = FALSE,
       "max_psrf" = max_psrf,
       "burn_in" = burn_in,
       "combined_size" = combined_size)

}

#' performs parallel tempering
#' 
#' @param mcmc_list_in a list of lists: values, log likelihood etc. of parallel MCMC chains
#' @param temperatures numeric vector: temperatures of chains
#' @param offset integer: 0 or 1. 0 = swap chains 1 <-> 2, 3 <-> 4...
#' 1 = swap chains 2<->3, 4<->5...
#' 
#' @return a list of lists: values, log likelihood etc. of paralle chains after parallel tempering
#' @export
parallel_tempering <- function(mcmc_list_in, temperatures, offset){
  
  mcmc_list_out <- mcmc_list_in
  probabs <- vapply(mcmc_list_in, function(x) x$probab, double(1))
  if((offset + 1) <= (length(mcmc_list_in) - 1)){
    swap_ind <- seq(offset + 1, length(mcmc_list_in) - 1, by = 2)
    
    decide_if_swap <- function(x,y){
      delta <- (1 / temperatures[y] - 1 / temperatures[x]) *
        (probabs[x] - probabs[y])
      runif(1) <= exp(delta)
    }
    
    swaps <- vapply(swap_ind,function(x) decide_if_swap(x,x+1), logical(1))
    
    perform_swap <- function(mcmc_list_in, swap, swap_ind){
      mcmc_list <- mcmc_list_in[c(swap_ind, swap_ind + 1)]
      if(swap){
        mcmc_list_out <- mcmc_list
        mcmc_list_out[[1]][["probab"]] <- mcmc_list[[2]][["probab"]]
        mcmc_list_out[[2]][["probab"]] <- mcmc_list[[1]][["probab"]]
        mcmc_list_out[[1]][["current_pars"]] <- mcmc_list[[2]][["current_pars"]]
        mcmc_list_out[[2]][["current_pars"]] <- mcmc_list[[1]][["current_pars"]]
        mcmc_list_out
      } else {
        mcmc_list
      }
    }
    
    mcmc_list_out <- Map(function(x,y) perform_swap(mcmc_list_in,x,y),
                         swaps,swap_ind)
    mcmc_list_out <- unname(unlist(mcmc_list_out, recursive = FALSE))
  } else {
    swaps <- 0
  }
  list("swaps" = swaps, "mcmc_list" = mcmc_list_out)
}

#' wrapper for parLapply for cluster
#'
#' @param run_parallel: logical: if TRUE, use parLapply, else use lapply
#' @param x first argument of lapply
#' @param fun second argument of lapply
#' @return output arguments of lapply
parLapply_wrapper <- function(run_parallel,x,fun,...){
  if(run_parallel){
    sys_info <- Sys.info()
    if(sys_info$sysname == "Windows"){
      parLapply(cl = NULL, x, fun, ...)
    } else {
      mclapply(x, fun, ..., mc.cores = length(x))
    }
  } else {
    lapply(x, fun, ...)
  }
}

#' calibrate temperatures for parallel chains
#'
#' @param temperatures vector of length n: current temperatures of chains
#' @param swap_ratio vector of length n - 1: (proportion of accepted swaps
#' out of proposed swaps)/2
#' the factor of 2 arises because of the way the swaps are recorded, and how
#' we alternate between swapping 1<->2, 3<->4.... and 2<->3, 4<->5...
#' @return vector of length n: new temperatures of chains
#'
calibrate_temperatures <- function(temperatures,swap_ratio) {
  return(temperatures)
  diff_temp <- diff(temperatures)
  # find chains between which the swap ratio is too large
  too_large = swap_ratio > .2 # note factor of 2 from main text -- see above
  # find chains between which the swap ratio is too small
  too_small = swap_ratio < .05
  # adjust differences between temperatures accordingly
  diff_temp = diff_temp*(too_large*1.5 + too_small*.75 + (1-too_large-too_small))
  # reconstruct temperatures from their differences
  cumsum(c(temperatures[1],diff_temp))
  
}
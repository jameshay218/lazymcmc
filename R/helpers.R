#' Model parameter table check
#'
#' Takes the table of paramters used for model calculation and MCMC fitting, and returns a string error if any of the inputs are incorrectly formatted
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @return a list of strings with the errors of the parameter table
#' @export
param_table_check <- function(parTab){
    needed_names <- c("values","fixed","steps","lower_bound","upper_bound")
    errors <- FALSE

    if(class(parTab) != "data.frame"){
        class_error <- "Error in parTab - must be data frame"
        errors <- list(errors, class_error)
    }
    if(!(any(needed_names %in% colnames(parTab)))){
        colnames_error <- paste0("Error in columns - missing ",
                                 needed_names[!(needed_names %in% colnames(parTab))],
                                 ". Columns should be value, lower_bound, upper_bound, steps, fixed")
        errors <- list(errors, colnames_error)
    }
    if(any(!is.numeric(parTab$values))){
        numeric_error <- "Error in values column - non numeric values provided"
        errors <- list(errors, numeric_error)
    }
    if(any(!is.numeric(parTab$lower_bound))){
        numeric_error <- "Error in lower_bound column - non numeric values provided"
        errors <- list(errors, numeric_error)
    }
    if(any(!is.numeric(parTab$upper_bound))){
        numeric_error <- "Error in upper_bound column - non numeric values provided"
        errors <- list(errors, numeric_error)
    }
    if(any(!(parTab$fixed %in% c(0, 1)))){
        numeric_error <- "Error in steps column - must be 0 (free) or 1 (fixed)"
        errors <- list(errors, numeric_error)
    }
    if(any(parTab$steps > 1 | parTab$steps < 0)){
        numeric_error <- "Error in steps column - initial step size must be between 0 and 1"
        errors <- list(errors, numeric_error)
    }
    if(any(parTab$upper_bound < parTab$lower_bound)){
        bound_error <- "Error in bounds - upper_bound lower than lower_bound"
        errors <- list(errors, bound_error)
    }

    if(any(parTab$values < parTab$lower_bound | parTab$values > parTab$upper_bound)){
        bound_error <- "Error in starting values - values outside of bounds"
        errors <- list(errors, bound_error)
    }

    if(all(parTab$fixed == 1)){
        fixed_error <- "Error in fixed parameters - all parameters fixed!"
        errors <- list(errors, fixed_error)
    }
    if(length(errors) > 1) errors[[1]] <- TRUE

    return(errors)
    
}


#' MCMC chain parameter check
#'
#' Checks that the given settings for the MCMC chain are correct
#' @param mcmcPars the vector of MCMC pars as expected by run_MCMC
#' @param mvrPars the list of multivarate MCMC pars as expected by run_MCMC
#' @return a list of errors
#' @export
mcmc_param_check <- function(mcmcPars, mvrPars){
    needed_names <- c("iterations","popt","opt_freq","adaptive_period","save_block")
    errors <- FALSE
    if(!(any(needed_names %in% names(mcmcPars)))){
        names_error <- paste0("Error in mcmcPars - missing ",
                              needed_names[!(needed_names %in% names(mcmcPars))])
        errors <- list(errors, names_error)
    }
    if(mcmcPars[["iterations"]] < 1){
        errors <- list(errors, "Error in iterations - less than 1 iteration specified")        
    }
    if(mcmcPars[["popt"]] < 0 | mcmcPars[["popt"]] > 1){
        errors <- list(errors, "Error in popt - invalid desired acceptance rate")
    }
    if(mcmcPars[["thin"]] > mcmcPars[["iterations"]]){
        errors <- list(errors, "Error in thin value - thinning more than number of iterations")
    }

    if(!is.null(mvrPars)){
        if(length(mvrPars) != 3){
            errors <- list(errors, "Error in mvrPars - list should have 3 elements")
        }
    }
    
    if("parallel_tempering_iter" %in% names(mcmcPars) && 
       (!("temperature" %in% names(mcmcPars)) || length(mcmcPars[["temperature"]]) < 2)){
      errors <- list(errors, "insufficient number of temperatures specified for parallel tempering")
    }
    
    if(length(errors) > 1) errors[[1]] <- TRUE
    return(errors)
}

#' data.table::fread with tryCatch
#' 
#' data.table::fread sometimes crashes.  If crashing, try read.csv then converting.
#' 
#' @param filename character vector of length 1 ending in .csv
#' @return the read data frame
try_fread <- function(filename, ...) {
  tryCatch ({
    data.table::fread(filename, ...)
  }, error = function(c) {
    data_table <- read.csv(filename, ...)
    data.table::as.data.table(data_table)
  })
}

#' Generate starting parameter table
#'
#' Generates a version of \code{par_tab} with random values between \code{lower_start} and \code{upper_start}
#' @param par_tab See \code{\link{example_par_tab}}
#' @return a data frame matching par_tab
#' @family mcmc
#' @examples
#' data(example_par_tab)
#' start_tab <- generate_start_tab(example_par_tab)
#' @export
generate_start_tab <- function(par_tab){
  for(i in 1:nrow(par_tab)){
    if(par_tab[i,"fixed"] == 0){
      par_tab[i, "values"] <- runif(1,par_tab[i,"lower_start"], par_tab[i, "upper_start"])
    }
  }
  return(par_tab)        
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

#' Read in MCMC chains
#'
#' Loads all available MCMC chains from the chosen working directory, allowing the user to specify properties of the return chain
#' @param location Either the full file path to the MCMC chain containing directory, or a vector of file paths with the MCMC chainsa
#' @param parTab the parameter table that was used to solve the model (mainly used to find which were free parameters)
#' @param unfixed Boolean, if TRUE, only returns free parameters
#' @param thin thin chain by this much
#' @param burnin number of iterations to discard
#' @param multi if TRUE, looks for chains generated using the multivariate sampler. Otherwise, looks for univariate sampled chains.
#' @param chainNo if TRUE, adds the chain number to the MCMC chain as a column
#' @param PTchain DEV - if TRUE, looks for chains generated by the parallel tempering algorithm
#' @return a list containing both the MCMC chains appended to each other, and an MCMC list.
#' @export
load_mcmc_chains <- function(location="",parTab,unfixed=TRUE, thin=1,
                             burnin=100000, multi=TRUE, chainNo=FALSE,
                             PTchain=FALSE){
  if(length(location) == 1 & dir.exists(location)){
    print("Reading in chains from directory")
    if(multi){
      chains <- Sys.glob(file.path(location,"*multivariate_chain.csv"))
    } else {
      chains <- Sys.glob(file.path(location,"*univariate_chain.csv"))
    }
    if(PTchain){
      chains <- Sys.glob(file.path(location,"*_chain.csv"))
    }
  } else {
    print("Reading in chains from filepaths")
    if(length(location) == 1){
      chains <- list(location)
    } else {
      chains <- as.list(location)
    }
  }
  
  print(chains)
  if(length(chains) < 1){
    message("Error - no chains found")
    return(NULL)
  }
  
  ## Read in the MCMC chains with fread for speed
  read_chains <- lapply(chains,data.table::fread,data.table=FALSE)
  
  ## Thin and remove burn in
  read_chains <- lapply(read_chains, function(x) x[seq(1,nrow(x),by=thin),])
  read_chains <- lapply(read_chains,function(x) x[x$sampno > burnin,])
  print(lapply(read_chains, nrow))
  
  if(chainNo){
    for(i in 1:length(read_chains)) read_chains[[i]]$chain <- i
  }
  
  ## Get the estimated parameters only
  if(unfixed){
    fixed <- parTab$fixed
    read_chains <- lapply(read_chains, function(x) x[,c(which(fixed==0)+1,ncol(x))])
  }
  
  ## Try to create an MCMC list. This might not work, which is why we have a try catch
  list_chains <- tryCatch({
    tmp_list <- lapply(read_chains,coda::as.mcmc)
    tmp_list <- coda::as.mcmc.list(tmp_list)
  }, warning = function(w){
    print(w)
    NULL
  }, error = function(e){
    print(e)
    NULL
  },
  finally = {
    tmp_list
  })
  
  chain <- coda::as.mcmc(do.call("rbind",read_chains))
  return(list("list"=list_chains,"chain"=chain))
}

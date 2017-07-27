#' Model parameter table check
#'
#' Takes the table of paramters used for model calculation and MCMC fitting, and returns a string error if any of the inputs are incorrectly formatted
#' @param parTab the parameter table controlling information such as bounds, initial values etc
#' @return a list of strings with the errors of the parameter table
#' @export
#' useDynLib lazymcmc
parameter_table_check <- function(parTab){
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
    
    if(length(errors > 1)) errors[[1]] <- TRUE
    
    return(errors)
    
}


#' MCMC chain parameter check
#'
#' Checks that the given settings for the MCMC chain are correct
#' @param mcmcPars the vector of MCMC pars as expected by \link{\code{run_MCMC}}
#' @param mvrPars the list of multivarate MCMC pars as expected by \link{\code{run_MCMC}}
#' @return a list of errors
#' @export
mcmc_param_check <- function(mcmcPars, mvrPars){
    needed_names <- c("values","fixed","steps","lower_bound","upper_bound")
    errors <- FALSE
    if(!(any(needed_names %in% names(mcmcPars)))){
        names_error <- paste0("Error in mcmcPars - missing ",
                              needed_names[!(needed_names %in% names(mcmcPars))])
        errors <- list(errors, names_error)
    }
    if(mcmcPars["iterations"] < 1){
        errors <- list(errors, "Error in iterations - less than 1 iteration specified")        
    }
    if(mcmcPars["popt"] < 0 | mcmcPars["popt"] > 1){
        errors <- list(errors, "Error in popt - invalid desired acceptance rate")
    }
    if(mcmcPars["thin"] > mcmcPars["iterations"]){
        errors <- list(errors, "Error in thin value - thinning more than number of iterations")
    }

    if(!is.null(mvrPars)){
        if(length(mvrPars) != 3){
            errors <- list(errors, "Error in mvrPars - list should have 3 elements")
        }
    }
    if(length(errors > 1)) errors[[1]] <- TRUE
    return(errors)
}

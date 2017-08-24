#' BIC calculation
#' 
#' Given an MCMC chain, calculates the BIC
#' @param chain the MCMC chain to be tested
#' @param parTab the parameter table used for this chain
#' @param dat the data
#' @return a single BIC value
#' @export
calculate_BIC <- function(chain, parTab, dat){
    n <- nrow(dat)*ncol(dat)
    maxLik <- -2*max(chain$lnlike)
    B <- length(parTab[parTab$fixed==0,"values"])*log(n)
    return(maxLik + B)
}

calculate_WAIC <- function(chain, parTab, dat, f, N){
    expectation_posterior <- 0
    tmp <- matrix(nrow=chain,ncol=nrow(dat)*ncol(dat))
    for(i in 1:nrow(chain)){
        pars <- zikaProj::get_index_pars(chain,i)
        y <- f(pars)
        index <- 1
        for(j in 1:nrow(dat)){
            for(x in 1:ncol(dat)){
                wow <-norm_error(y[j,x],dat[j,x],pars["S"],pars["MAX_TITRE"])
                expectation_posterior <- expectation_posterior + log(wow)
                tmp[i, index] <- wow
                index <- index + 1                  
            }
        }
    }
    lppd <- sum(log(colMeans(tmp)))
    pwaic1 <- 2*(lppd - expectation_posterior)
}

#' @export
calc_DIC <- function(lik.fun,chain){
    D.bar <- -2*mean(chain$lnlike)
    theta.bar <- as.numeric(summary(as.mcmc(chain[,2:(ncol(chain)-1)]))$statistics[,"Mean"])
   # print(theta.bar)
    D.hat <- -2*lik.fun(theta.bar)
    pD <- D.bar - D.hat
    pV <- var(-2*chain$lnlike)/2
    list(DIC=2*pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}

#' @export
calculate_AIC <- function(chain, parTab){
    k <- nrow(parTab[parTab$fixed == 0,])
    AIC <- 2*k - 2*(max(chain$lnlike))
    return(AIC)
}

#' calculate and plot profile likelihood around a point in parameter space
#' 
#' calculate and plot profile likelihood around a point in parameter space
#' 
#' @param stuck_values a named numeric vector of parameter values (length m)
#' around which to construct the profile likelihood
#' @param f a function handle to calculate the log likelihood. 
#' f <- CREATE_POSTERIOR_FUNC(parTab, data, PRIOR_FUNC)
#' @param fixed a numeric vector of length m indicating whether a parameter is 
#' fixed (1) or fitted (0)
#' @param n_values a numeric vector of length 1 indicating the number of points 
#' to scan over in 1D parameter space for each profile. Must take an integer value.
#' @param range a numeric vector of length 1 indicating the range over which to
#' scan in 1D parameter space for each profile (so we go range/2 on each side of
#' the central value)
#' @param par_names_plot a character vector of length m to label the horizontal
#' axes of the profile plots
#' @return a list of ggplot objects containing the profile likelihood plots
#' @export 
profile_likelihood_all <- function(stuck_values, f, fixed, n_values, range, par_names_plot){
  f_lik <- function(pars){
    out <- f(pars)
    if(is.atomic(out)){
      out
    } else {
      out$lik
    }
  }
  
  par_names <- names(stuck_values)
  stuck_values <- as.numeric(stuck_values)  
  stuck_values_df <- data.frame(matrix(rep(stuck_values,n_values),nrow = n_values,byrow = TRUE))
  colnames(stuck_values_df) <- par_names
  
  profile_likelihood <- function(col){
    stuck_values_df[,col] <- seq(stuck_values_df[1,col] - range / 2,
                                 stuck_values_df[1,col] + range / 2, 
                                 length.out = n_values)
    stuck_values_df$lnlike <- apply(stuck_values_df,1,f_lik)
    g <- ggplot(stuck_values_df,aes_string(x = par_names[col], y = "lnlike"))
    g <- g + geom_point() +
      geom_vline(xintercept = stuck_values[col]) +
      theme_bw() +
      theme(text = element_text(size = 24)) +
      xlab(TeX(par_names_plot[col])) + ylab("LL")
    g
  }
  
  unfixed <- which(fixed == 0)
  lapply(unfixed, profile_likelihood)
}

#' calculate percentiles for fitted parameters
#' 
#' calculate percentiles for fitted parameters
#' 
#' @param chain the MCMC chain: a data frame with n columns, where n is the number
#' of fitted parameters
#' @param prctiles a numeric vector of length m containing the percentiles to be calculated
#' (between 0 and 1)
#' @param par_names_plot parameter names for the output data frame
#' @return a data frame with n rows and m columns containing the percentiles
#' @export 
print_prctiles <- function(chain, prctiles = c(.025,.5,.975), par_names_plot){
  prctile_table <- lapply(chain,function(x) quantile(x,prctiles))
  prctile_table <- t(as.data.frame(prctile_table))
  rownames(prctile_table) <- par_names_plot
  col_names <- as.character(prctiles*100)
  col_names <- trimws(format(col_names, digits = 3))
  col_names <- paste0(col_names,"\\%")
  colnames(prctile_table) <- col_names
  prctile_table
}
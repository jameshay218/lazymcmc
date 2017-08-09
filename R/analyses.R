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



GINAR.CLS.U <- function (eta, x, p, transform=TRUE) {

    ks <- 1:p
    n  <- length(x)

    if (transform) {

        alpha  <- GINAR.kappa.to.alpha(eta[ks])
        mu.eps <- exp(eta[p+1])
        
    } else {

        alpha  <- eta[ks]
        mu.eps <- eta[p+1]
    }

    ts <- (p+1):n
    us <- x[ts] - GINAR.cond.means(x, alpha, mu.eps)

    sum(us^2)
}



GINAR.CLS.S.sigma.eps <- function (sigma.eps, alpha, mu.eps, x, process="INAR") {

    ts <- (length(alpha)+1):length(x)
    
    us <- x[ts] - GINAR.cond.means(x, alpha, mu.eps)

    vs <- GINAR.cond.vars(x, alpha, sigma.eps, process)

    sum( (us^2 - vs)^2  )
}


GINAR.CLS.S.r <- function (r, alpha, mu.eps, x, process="INAR") {

    sigma.eps <- GINAR.r.to.sigma.eps(r, mu.eps)
    
    ts <- (length(alpha)+1):length(x)
    
    us <- x[ts] - GINAR.cond.means(x, alpha, mu.eps)

    vs <- GINAR.cond.vars(x, alpha, sigma.eps, process)

    sum( (us^2 - vs)^2  )
}



GINAR.CLS <- function (x, p, process="INAR", r.max=10, alpha0, mu.eps0) {
    ## ======================================================================
    ## Estimate the parameters of a GINAR('p') process for time series
    ## data in the vector 'x' using the conditional least square method.
    ##
    ## Currently supports the following values for 'process':
    ## 
    ## INAR      : binomial thinning with Poisson innovations;
    ## INAR-NB   : binomial thinning with negative binomial innovations;
    ## Geom-INAR : geometric thinning with Poisson innovations.
    ## ======================================================================

    if (missing(alpha0)) {
        alpha0 <- rep(1e-4, p)
    }
    if (missing(mu.eps0)) {

        mu.eps0 <- GINAR.mu.eps0(x, alpha0)
    }
    eta0 <- c(GINAR.alpha.to.kappa(alpha0), log(mu.eps0))

    if (process=="INAR" || process=="Geom-INAR") {
        ## Binomial thinning and Poisson innovations or
        ## Geometric thinning and Poission innovations

        U.est <- optim(par=eta0, fn=GINAR.CLS.U, method="BFGS",
                       x=x, p=p, transform=TRUE)       
        
        alpha.hat <- GINAR.kappa.to.alpha(U.est$par[1:p])
        
        mu.eps.hat <- exp(U.est$par[p+1])
        
        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             process = process)
        
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations
                
        U.est <- optim(par=eta0, fn=GINAR.CLS.U, method="BFGS",
                       x=x, p=p, transform=TRUE)       
        
        alpha.hat <- GINAR.kappa.to.alpha(U.est$par[1:p])
        
        mu.eps.hat <- exp(U.est$par[p+1])

        r.est <- optimize(GINAR.CLS.S.r, c(1e-5, r.max), alpha=alpha.hat,
                          mu.eps=mu.eps.hat, x=x, process=process)

        r.hat <- r.est$minimum

        sigma.eps.hat <- GINAR.r.to.sigma.eps(r.hat, mu.eps.hat)
        
        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             sigma.eps.hat = sigma.eps.hat,
             r.hat = r.hat,
             process = process)
        
    } else {

        stop(paste("The GINAR process:", process, "is not supported"))
    }
}

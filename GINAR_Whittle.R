

GINAR.Whittle.lW <- function (eta, sp, p, process, transform=TRUE) {

    fs <- sp$freq[-1]
    
    if (transform) {

        alpha  <- GINAR.kappa.to.alpha(eta[1:p])
        mu.eps <- exp(eta[p+1])
        
    } else {

        alpha  <- eta[1:p]
        mu.eps <- eta[p+1]
    }

    if (process=="INAR" || process=="Geom-INAR") {

        sigma.eps <- mu.eps

    } else if (process=="INAR-NB") {

        if (transform) {
            
            r <- exp(eta[p+2])
            
        } else {

            r <- eta[p+2]
        }

        sigma.eps <- GINAR.r.to.sigma.eps(r, mu.eps)
    }

    sdfs <- GINAR.sdf(fs, alpha, mu.eps, sigma.eps, process)
    
    2 * sum(log(sdfs) + sp$spec[-1]/sdfs)
}







GINAR.Whittle <- function (x, sp, p, process="INAR", alpha0, mu.eps0) {
    ## ======================================================================
    ## Estimate the parameters of a GINAR('p') process for time series
    ## data in the vector 'x' using the Whittle likelihood method.
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
        
        P.est <- optim(par=eta0, fn=GINAR.Whittle.lW, method="BFGS",
                       sp=sp, p=p, process=process, transform=TRUE)       
        
        alpha.hat <- GINAR.kappa.to.alpha(P.est$par[1:p])
        
        mu.eps.hat <- exp(P.est$par[p+1])
        
        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             process = process)
        
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations

        eta0   <- c(eta0, log(1))
        
        P.est <- optim(par=eta0, fn=GINAR.Whittle.lW, method="BFGS",
                       sp=sp, p=p, process=process, transform=TRUE)       
        
        alpha.hat <- GINAR.kappa.to.alpha(P.est$par[1:p])
        
        mu.eps.hat <- exp(P.est$par[p+1])

        r.hat <- exp(P.est$par[p+2])        

        sigma.eps.hat = mu.eps * (1.0 + r.hat * mu.eps.hat)        
        
        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             sigma.eps.hat = sigma.eps.hat,
             r.hat = r.hat,
             process = process)
        
    } else {

        stop(paste("The GINAR process:", process, "is not supported"))
    }
}



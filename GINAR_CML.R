

GINAR.CML.l <- function (eta, x, p, process, transform=TRUE) {

    if (transform) {

        alpha  <- GINAR.kappa.to.alpha(eta[1:p])
        mu.eps <- exp(eta[p+1])
        
    } else {

        alpha  <- eta[1:p]
        mu.eps <- eta[p+1]
    }

    n  <- length(x)
    ts <- (p+1):n

    zs   <- exp(1i * GINAR.gw$absc)
    zsm1 <- zs - 1.0
    ys   <- 1.0/zs

    if (process=="INAR") {
        ## Binomial thinning and Poisson innovations
        
        ws    <- GINAR.gw$weights * exp(mu.eps * zsm1)
        inner <- outer(alpha, zs) - alpha + 1.0
                
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations
        
        if (transform) {
            
            r <- exp(eta[p+2])
            
        } else {
            
            r <- eta[p+2]
        }
        
        ws    <- GINAR.gw$weights * (1.0 - r * mu.eps * zsm1)^(-1.0/r)
        inner <- outer(alpha, zs) - alpha + 1.0
        
    } else if (process=="Geom-INAR") {
        ## Geometric thinning and Poisson innovations
        
        ws    <- GINAR.gw$weights * exp(mu.eps * zsm1)
        inner <- 1.0/(outer(-alpha, zsm1) + 1.0)
    }

    ws0 <- ws / zsm1

    terms <- rep(NA, n - p)
    
    for (t in ts) {
        
        if (x[t]>0) {
            
            aa <- ws
            for (j in 1:p) {
                
                aa <- aa * inner[j,]^x[t-j]
            }

            terms[t-p] <- sum(Re(aa * ys^x[t]))
            
        } else {
            
            aa <- ws0
            for (j in 1:p) {
                
                aa <- aa * inner[j,]^x[t-j]
            }

            terms[t-p] <- pi/2 - sum(Re(aa))
        }
    }

    -sum(log(pmax(terms, .Machine$double.eps))) #+ (n-p)*log(pi)
}




GINAR.CML <- function (x, p, process="INAR", alpha0, mu.eps0) {
    ## ======================================================================
    ## Estimate the parameters of a GINAR('p') process for time series
    ## data in the vector 'x' using the conditional maximum likelihood
    ## method.
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
        
        P.est <- optim(par=eta0, fn=GINAR.CML.l, method="BFGS",
                       x=x, p=p, process=process, transform=TRUE)
        
        alpha.hat <- GINAR.kappa.to.alpha(P.est$par[1:p])
        
        mu.eps.hat <- exp(P.est$par[p+1])
        
        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             process = process)
        
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations

        eta0   <- c(eta0, log(1))
        
        P.est <- optim(par=eta0, fn=GINAR.CML.l, method="BFGS",
                       x=x, p=p, process=process, transform=TRUE)       
        
        alpha.hat <- GINAR.kappa.to.alpha(P.est$par[1:p])
        
        mu.eps.hat <- exp(P.est$par[p+1])

        r.hat <- exp(P.est$par[p+2])        

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

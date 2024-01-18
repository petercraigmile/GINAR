

GINAR.YW <- function (x, p, process="INAR") {
    ## ======================================================================
    ## Estimate the parameters of a GINAR('p') process for time series
    ## data in the vector 'x' using the Yule-Walker method.
    ##
    ## Currently supports the following values for 'process':
    ## 
    ## INAR      : binomial thinning with Poisson innovations;
    ## INAR-NB   : binomial thinning with negative binomial innovations;
    ## Geom-INAR : geometric thinning with Poisson innovations.
    ## ======================================================================

    yw.est <- ar.yw(x, aic = FALSE, order = p)

    alpha.hat <-  yw.est$ar

    x.bar <- mean(x)
    
    mu.eps.hat <- x.bar * ( 1.0 - sum(alpha.hat) )

    if (process=="INAR" || process=="Geom-INAR") {
        ## Binomial thinning and Poisson innovations or
        ## Geometric thinning and Poission innovations

        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             process = process)
        
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations
        
        sample.acvf   <- acf(x, lag.max=p, type="cov", plot=FALSE)$acf
        
        V.p.hat       <- sample.acvf[1] - sum(alpha.hat * sample.acvf[-1])

        beta.hat      <- alpha.hat * (1.0 - alpha.hat)

        sigma.eps.hat <- V.p.hat - x.bar * sum(beta.hat)

        r.hat         <- GINAR.sigma.eps.to.r(sigma.eps.hat, mu.eps.hat)

        list(alpha.hat  = alpha.hat,
             mu.eps.hat = mu.eps.hat,
             sigma.eps.hat = sigma.eps.hat,
             r.hat = r.hat,
             process = process)
        
    } else {

        stop(paste("The GINAR process:", process, "is not supported"))
    }
}



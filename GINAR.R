


GINAR.sim <- function (n, alpha, mu.eps, process="INAR", r=NA, n.start=1000) {

    ## Error checking to implement:
    ## 0 < alpha[i] < 1 for each i with sum(alpha) < 1
    ## if process is "INAR-NB", then there must be a value of r

    p <- length(alpha)    
    N <- n + n.start

    ## Generate the innovations
    if (process=="INAR" || process=="Geom-INAR") {

        x <- rpois(N, mu.eps)
        
    } else if (process=="INAR-NB") {

        x <- rnbinom(n=N, size=r, mu=mu.eps)
        
    } else {
        
        stop(paste("The GINAR process:", process, "is not supported"))
    }

    if (process=="INAR" || process=="INAR-NB") { ## Either we use binomial thinning

        for (t in (p+1):N) {
            for (j in 1:p) {
                
                x[t] <- x[t] + rbinom(1, x[t-j], alpha[j])
            }
        }

    } else { ## or Geometric thinning

        for (t in (p+1):N) {    
            for (j in 1:p) {
                
                if (x[t-j] > 0) {
                    
                    x[t] <- x[t] + rnbinom(1, size=x[t-j], prob=1/(1+alpha[j]))
                }
            }   
        }
    }

    x[(n.start+1):N]
}




GINAR.beta <- function (alpha, process) {

    if ((process=="INAR") || (process=="INAR-NB")) {
        ## Binomial thinning and Poisson or negative binomial innovations

        alpha * (1.0 - alpha)
        
    } else if (process=="Geom-INAR")  {
        ## Binomial thinning and negative binomial innovations

        alpha * (1.0 + alpha)
        
    } else {
        
        stop("beta function not implemented for this process")
    }
}



GINAR.cond.means <- function (x, alpha, mu.eps) {

    p  <- length(alpha)
    ks <- 1:p    
    ts <- (p+1):length(x)

    sapply(ts, function (t) {

        sum(alpha * x[t-ks]) + mu.eps
    })
}




GINAR.cond.vars <- function (x, alpha, sigma.eps, process="INAR") {

    p  <- length(alpha)
    ks <- 1:p    
    ts <- (p+1):length(x)

    beta <- GINAR.beta(alpha, process)

    sapply(ts, function (t) {

        sum(beta * x[t-ks]) + sigma.eps
    })
}




GINAR.chf <- function (u, past.x, alpha, mu.eps, r, process="INAR") {
    ## ======================================================================
    ## Calculate phi_{X_t | X_{t-1}, \ldots, X_{t-p}}(u)
    ## only works for scalar 'u'
    ## ======================================================================

    w <- exp(1i * u)
    B <- mu.eps * (w - 1.0)
    
    if (process=="INAR") {
        ## Binomial thinning and Poisson innovations

        exp(B) * prod((1.0 - alpha + alpha * w)^past.x)
        
    } else if (process=="INAR-NB") {
        ## Binomial thinning and negative binomial innovations

        (1.0 - B)^(-1/r) * prod((1.0 - alpha + alpha * w)^past.x)
        
    } else if (process=="Geom-INAR") {
        ## Geometric thinning and Poisson innovations

        exp(B) * prod( (1.0 - alpha * (w - 1.0) )^(-past.x) )
    }
}




GINAR.sdf <- function (freqs, alpha, mu.eps, sigma.eps,
                       process = "INAR", delta.t = 1) {
    
    ws <- -2 * pi * delta.t * freqs
    
    js <- seq(length(alpha))
    
    reals <- sapply(ws, function(w, js, alpha)
        1 - sum(alpha * cos(js * w)), js = js, alpha = alpha)
    
    imags <- sapply(ws, function(w, js, alpha)
        sum(alpha * sin(js * w)), js = js, alpha = alpha)

    mu.X   <- mu.eps / (1.0 - sum(alpha))

    beta   <- GINAR.beta(alpha, process)

    sigma2 <- sigma.eps + mu.X * sum(beta)

    (sigma2 * delta.t)/(reals * reals + imags * imags)
}



## Define the Gauss weights used in the CML method

GINAR.gw <- Gauss.weights(0, pi, 300)



## Define some transformations

GINAR.alpha.to.kappa <- function (alpha) {

    log(alpha / (1.0 - sum(alpha)))
}



GINAR.kappa.to.alpha <- function (kappa) {

    exp.kappa <- exp(kappa)

    exp.kappa / (1.0 + sum(exp.kappa))
}




GINAR.sigma.eps.to.r <- function (sigma.eps, mu.eps) {
    
    ( sigma.eps - mu.eps ) / mu.eps^2
}



GINAR.r.to.sigma.eps <- function (r, mu.eps) {
    
    mu.eps + r * mu.eps^2
}




GINAR.mu.eps0 <- function (x, alpha0, delta=1e-10) {

    alpha.star <- ifelse(alpha0<delta, delta,
                         ifelse(alpha0>(1-delta), 1-delta, alpha0))

    mean(x) * (1 - sum(alpha.star))
}

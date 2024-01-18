
Gauss.weights <- function (lower, upper, n)
  ## ======================================================================
  ## Calculates the Gaussian weights for Legendre quadrature
  ## i.e. \int_{lower}^{upper} f(x) dx is approximated by 
  ##      sum_{i=1}^n weights_i f(abscissa_i)
  ##
  ## Requires: statmod R library.
  ## =======================================================================
{
  require(statmod)
  mid  <- 0.5 * (lower + upper)
  hrng <- 0.5 * (upper - lower)
  gq   <- gauss.quad(n)

  list(abscissa = gq$nodes * hrng + mid,
       weights  = hrng * gq$weights,
       n        = n)
}

## Last updated 2007-05-09, pfc@stat.osu.edu
## Created 2002-06-01, pfc@stat.osu.edu


Fourier.frequencies <- function (N, delta.t=1) {
  ## ======================================================================
  ## Purpose : Calculates all the Fourier frequencies obtained when
  ##           performing spectral analysis on a time series of length 'N',
  ##           with sampling interval 'delta.t'.
  ## Updated : pfc@stat.osu.edu, April 2007.
  ## ======================================================================
  
  (seq(floor(N/2+1))-1) / (N * delta.t)
}



spect <- function (spec, N, delta.t=1)
  ## ======================================================================
  ## Purpose : Create a class of type "spect" from a spectrum 'spec' of a 
  ##           time series of length 'N' regularly sampled at interval 'delta.t'.
  ## ======================================================================
{
  freq <- Fourier.frequencies(N, delta.t)
  n.Fourier <- length(freq)

  if (is.null(dim(spec)))
    actual.spec <- spec[1:n.Fourier]
  else
    actual.spec <- spec[1:n.Fourier,] 
  
  structure(list(freq      = freq,
                 spec      = actual.spec,
                 delta.t   = delta.t,
                 N         = N,
                 n.Fourier = n.Fourier),
            class = "spect")
}




pgram <- function (x, delta.t=deltat(x))
{
  N <- length(x)
  spect(Mod(fft(x))^2/N, N, delta.t)
}


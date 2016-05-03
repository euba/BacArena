% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Substance.R
\docType{methods}
\name{diffusePDE}
\alias{diffusePDE}
\alias{diffusePDE,Substance-method}
\title{Function for diffusion of the Substance matrix}
\usage{
diffusePDE(object, init_mat, gridgeometry, lrw = NULL, tstep)

\S4method{diffusePDE}{Substance}(object, init_mat, gridgeometry, lrw = NULL,
  tstep)
}
\arguments{
\item{object}{An object of class Substance.}

\item{init_mat}{A matrix with values to be used by the diffusion.}

\item{gridgeometry}{A list specifying the geometry of the Arena}

\item{lrw}{A numeric value needed by solver to estimate array size (by default lwr is estimated in simEnv() by the function estimate_lrw())}

\item{tstep}{A numeric value giving the time step of integration}
}
\description{
The generic function \code{diffusePDE} implements the diffusion by the solving diffusion equation.
}
\details{
Partial differential equation is solved to model 2d diffusion process in the arena.
}
\examples{
arena <- Arena(n=100, m=100, stir=FALSE, Lx=0.025, Ly=0.025)
sub <- Substance(n=100,m=100,smax=0,name='test', difspeed=0.1, 
                 gridgeometry=arena@gridgeometry) #initialize test substance
sub@diffmat[ceiling(100/2),ceiling(100/2)] <- 40
diffusePDE(sub, init_mat=as.matrix(sub@diffmat),
           gridgeometry=arena@gridgeometry, tstep=arena@tstep)
}
\seealso{
\code{\link{Substance-class}} and \code{\link{diffuseR}}
}

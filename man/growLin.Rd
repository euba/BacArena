% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/Organism.R
\name{growLin}
\alias{growLin}
\title{Function for letting organisms grow linearly}
\usage{
growLin(object, growth)
}
\arguments{
\item{object}{An object of class Organisms.}

\item{growth}{A number indicating the current biomass, which has to be updated.}
}
\value{
Returns the updated biomass of the organisms of interest.
}
\description{
The generic function \code{growLin} implements a growth model of organisms in their environment.
}
\details{
Linear growth of organisms is implemented by adding the calculated growthrate by \code{optimizeLP} to the already present growth value.
}
\examples{
\dontrun{
ecore <- model #get Escherichia coli core metabolic model
org <- Organism(ecore,deathrate=0.05,duplirate=0.5,
           growthlimit=0.05,growtype="exponential") #initialize a organism
growLin(org,1)
}
}
\seealso{
\code{\link{Organism-class}} and \code{\link{optimizeLP}}
}

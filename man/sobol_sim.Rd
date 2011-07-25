\name{sobol_sim}
\alias{sobol_sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulates an SBML model given the model and sobol set matrix}
\description{
  This function converts an SBML level 2 file into a corresponding R model structure of class \code{SBML}.
}
\usage{
sobol_sim(model, sobolset, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  An SBML level 2 model input file.}
  \item{sobolset}{ Sobol Matrix Set.}
  \item{time}{  Time series for simulation.}
}
\details{
A limited subset of SBML level 2 models is currently supported, e.g. events and 
function definitions are not covered.
}
\value{A corresponding SBML model object in R. }
\author{ Vishak Venkateswaran}
\note{This function replaces \code{read.SBML} of older versions. }

\seealso{ \code{\link{readSBMLR}} }
\examples{
library(SBMLR)
library(odesolve)
sobolmatrix <- sobolset(path_to_sobol)
output <- sobol_sim(model, sobolmatrix, time)
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

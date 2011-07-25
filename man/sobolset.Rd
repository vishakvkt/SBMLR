\name{sobolset}
\alias{sobolset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulates an SBML model given the model and sobol set matrix}
\description{
  This function creates the Sobol Matrix, each row is a set and each column is a parameter \code{SBML}.
}
\usage{
sobolset(param_path, parameter_count=0, no_of_sets=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param_path}{  An SBML level 2 model input file.}
  \item{parameter_count}{ Optional count of number of parameters.}
  \item{no_of_sets}{Number of SobolSets to simulate.}
}
\details{
Function does not depend on the SBML parser. This is an independent entity.
}
\value{A data frame matrix containing the set of model parameters. }
\author{ Vishak Venkateswaran}
\note{This function creates the sobol matrix for the next step of the simulation. }

\seealso{ \code{\link{readSBMLR}} }
\examples{
library(SBMLR)
library(odesolve)
sobolmatrix <- sobolset(path_to_sobol)
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

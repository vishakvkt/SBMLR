\name{area_matrix}
\alias{area_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates area under the curve given the simulation results and model.}
\description{
Calculates area under the curve for the species in the model.
}
\usage{
area_matrix(model, model_result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{  An SBML level 2 model input file.}
  \item{model_result}{The set of simulation results.}
}
\details{
This is also independent of the parser. IT simply calculates area under the curve for every species in the model.
}
\value{matrix containing values. }
\author{ Vishak Venkateswaran }
\note{This function replaces \code{read.SBML} of older versions. }

\seealso{ \code{\link{readSBMLR}} }
\examples{
library(SBMLR)
library(odesolve)
areamatrix <- area_matrix(model, model_result)
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

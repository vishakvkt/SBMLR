\name{readSBML}
\alias{readSBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert an SBML file into an R model object of class SBML}
\description{
  This function converts an SBML level 2 file into a corresponding R model structure of class \code{SBML}.
}
\usage{
readSBML(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{  An SBML level 2 model input file.}
}
\details{
A limited subset of SBML level 2 models is currently supported, e.g. events and 
function definitions are not covered.
}
\value{A corresponding SBML model object in R. }
\author{ Tom Radivoyevitch }
\note{This function replaces \code{read.SBML} of older versions. }

\seealso{ \code{\link{readSBMLR}} }
\examples{
library(SBMLR)
library(odesolve)
curtoX=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
curtoR=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
curtoX==curtoR
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

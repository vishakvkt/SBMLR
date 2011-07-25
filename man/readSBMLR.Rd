\name{readSBMLR}
\alias{readSBMLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert an SBMLR file into an R model object of class SBML}
\description{
  This function converts an SBMLR model definition in filename into a corresponding returned SBML model structure.
}
\usage{
readSBMLR(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{  An SBMLR model definition file. }
}
\details{
A limited subset of SBML level 2 models is currently supported, e.g. events and 
function definitions are not covered.
}
\value{A corresponding \code{SBML} model object in R. }
\author{ Tom Radivoyevitch  }
\note{This function replaces the use of \code{source} in older versions of SBMLR. 
It includes rate law and rule string to function, expression and MathML mappings. }

\seealso{ \code{\link{readSBML}} }
\examples{
library(SBMLR)
library(odesolve)
curtoX=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
curtoR=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
curtoX==curtoR
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

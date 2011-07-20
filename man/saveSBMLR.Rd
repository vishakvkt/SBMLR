\name{saveSBMLR}
\alias{saveSBMLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Save an R model object of class SBML as an SBMLR file}
\description{
  This function converts \code{SBML} model object in R into an SBMLR model definition file. }
\usage{
saveSBMLR(model,filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ The \code{SBML} model object to be mapped into the SBMLR model definition file. }
  \item{filename}{ The file name of the destination SBMLR model definition file.  }
}
\value{
  No value returned. 
}
\references{
Radivoyevitch, T. A two-way interface between limited Systems Biology Markup Language and R. BMC Bioinformatics 5, 190 (2004).
}
\author{ Tom Radivoyevitch  }
\note{Similar to saveSBML, the file is written incrementally. }

\section{Warning }{SBML events and function definitions are NOT implemented.}

\seealso{ \code{\link{saveSBML}} }
\examples{
library(SBMLR)
library(odesolve)
curto=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))
saveSBMLR(curto,"curtoR.r")  
curtoR=readSBMLR("curtoR.r")  
curto==curtoR
summary(curtoR)
unlink("curtoR.r")
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

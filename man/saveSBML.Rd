\name{saveSBML}
\alias{saveSBML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Saves an R model object as an SBML file}
\description{
  This function converts a class \code{SBML} model object in R into an SBML level 2 file.
}
\usage{
saveSBML(model,filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ The model object in R.  }
  \item{filename}{ The name of the SBML file}
}
\details{
  The output file is SBML level 2.  
}
\value{
  No value returned. 
}
\references{
Radivoyevitch, T. A two-way interface between limited Systems Biology Markup Language and R. BMC Bioinformatics 5, 190 (2004).
}
\author{ Tom Radivoyevitch }
\note{For speed, the SBML file is written incrementally, rather than first built as a DOM in R and then saved using \code{xmlSave}.
}

\section{Warning }{SBML events and 
function definitions are NOT implemented.}

\seealso{ \code{\link{saveSBMLR}} }
\examples{
library(SBMLR)
library(odesolve)
curtoR=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
saveSBML(curtoR,"curtoR.xml") 
curtoX=readSBML("curtoR.xml")  
curtoX==curtoR
summary(curtoR)
unlink("curtoR.xml")
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line

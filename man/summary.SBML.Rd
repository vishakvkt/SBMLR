\name{summary.SBML}
\alias{summary.SBML}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get summary information from an SBML model}
\description{
  This function extracts information from a model of class \code{SBML} and returns it as a list. The list includes
  species and reaction information tabularized as data frames. 
}
\usage{
\method{summary}{SBML}(object,\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A model object of class \code{SBML} from which information is to be extracted.}
  \item{...}{ For compatibility with \code{summary} of the base package.}
 }
\details{no details}
\value{ A list containing the following elements
 \item{BC}{ A logical vector indicating which species are not state variables, 
  i.e. which species are boundary conditions or auxillary variables.}
  \item{y0}{ The initial state (boundary conditions excluded!). }
  \item{nStates}{ The length of the state vector, i.e. the number of system states.}
 \item{S0}{ The full set of species initial values.}
  \item{nReactions}{ The number of reactions.}
  \item{nSpecies}{ The number of species, including 
  states, boundary conditions and possibly auxillary variables such as the total concentration of 
  dihydofolate reductase in the morrison.r  model. }
  \item{incid}{
  The incidence/stoichiometry matrix. This usually contains ones and minus ones corresponding to fluxes either synthesizing or degrading (respectively) 
  a state variable
  chemical species. This matrix multiplied by the flux vector on its right yields 
  the corresponding concentration state variable time derivatives. }
 \item{species}{ Species information (i.e. names, ICs, BCs, and compartments) as a data frame.}
  \item{reactions}{ Reaction information tabularized as a dataframe, including string laws and initial fluxes.}
}

\note{The list output can be attached to immediately define many model variables of interest. }


\author{ Tom Radivoyevitch }

%\seealso{ \code{\link{==.SBML}}  }
%\seealso{ \code{\link{equateModels}}  }
\examples{
library(SBMLR)  
curto=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
summary(curto)

}
\keyword{ math }% at least one, from doc/KEYWORDS

\name{load_simulate_sensitivity}
\alias{load_simulate_sensitivity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulates an SBML model given the model path and sobol set path and time series and calculates Prcc using . This is a wrapper.}
\description{
	Wrapper function to load, simulate and retrieve the results.
}
\usage{
load_simulate_sensitivity(set_path, model_path, time)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{set_path}{ Path to the file containing the table of sobol sets.}
  \item{model_path}{ Path to the model you want to simulate.}
  \item{time}{  Time series for simulation.}
}
\details{
Level 2 parser seems to work accurately. Level 3 parser switch is also in place. Certain SBML entries have not been implemented. They will be skipped.
}
\value{A list of simulation result objects. }
\author{ Vishak Venkateswaran}
\note{This function replaces \code{read.SBML} of older versions. }

\seealso{ \code{\link{readSBMLR}} }
\examples{
library(senstivity)
library(SBMLR)
library(odesolve)
outputs <- load_simulate_sensitivity(path_to_sobol, path_to_model, minutes_to_run_simulation)
}
\keyword{ arith }% at least one, from doc/KEYWORDS
\keyword{ math }% __ONLY ONE__ keyword per line
